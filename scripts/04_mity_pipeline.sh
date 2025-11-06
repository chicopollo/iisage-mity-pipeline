#!/usr/bin/env bash
# Modified header to source config:

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/pipeline.conf"

# ===============================================================
# mity2_species_runall_v7.2.sh
# ---------------------------------------------------------------
# Universal MITY runner for non-model species.
# Pipeline:
#   1. Windowed FreeBayes calling
#   2. Ensure GT/DP/AD/RO/AO exist and are numeric
#   3. MITY normalise
#   4. Concatenate all windows
#   5. Zero-fill VAF
#   6. Convert GFF3 to BED for vcfanno (skips gene features)
#   7. Annotated report using vcfanno TOML (auto-generated)
#   8. Extract CSV from XLSX using ssconvert
#   9. Post-process CSV to add gene annotations from INFO field
# ---------------------------------------------------------------
# Author: Paul Decena (Bronikowski Lab, IISAGE)
# Version: 7.2 (2025-11)
# ===============================================================

set -euo pipefail
shopt -s nullglob
trap 'echo "ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR

# ---------------- USER PARAMETERS ----------------
DOCKER_IMAGE="${DOCKER_IMAGE:-$DOCKER_IMAGE}"  # From config
OUT="${OUT:-mity_chunked}"  # Keep as is
REPORT_DIR="${REPORT_DIR:-mity_output}"  # Keep as is
WSIZE="${WSIZE:-$MITY_WINDOW_SIZE}"  # From config
LIMIT_COV="${LIMIT_COV:-$MITY_COVERAGE_CAP}"  # From config
MIN_MQ="${MIN_MQ:-$MITY_MIN_MQ}"  # From config
MIN_BQ="${MIN_BQ:-$MITY_MIN_BQ}"  # From config
# -------------------------------------------------

echo "==> IISAGE-MITY species pipeline (v7.2)"
echo "Image:   $DOCKER_IMAGE"
echo "Window:  ${WSIZE}bp   CovCap: ${LIMIT_COV}   MQ: ${MIN_MQ}   BQ: ${MIN_BQ}"
echo

# ---- Check environment ----
for exe in bcftools tabix bgzip awk sed samtools htsfile; do
  command -v "$exe" >/dev/null || { echo "ERROR: missing '$exe' on PATH"; exit 1; }
done
command -v docker >/dev/null || { echo "ERROR: docker not found"; exit 1; }

if docker info >/dev/null 2>&1; then
  DOCKER="docker"
elif sudo -n docker info >/dev/null 2>&1; then
  DOCKER="sudo docker"
else
  sudo docker info >/dev/null 2>&1 || { echo "ERROR: cannot talk to Docker daemon"; exit 1; }
  DOCKER="sudo docker"
fi
DOCKER_RUN="$DOCKER run --rm -u $(id -u):$(id -g)"
DOCKER_SHELL="$DOCKER run --rm --entrypoint /bin/bash -u $(id -u):$(id -g)"
echo "Docker:  $DOCKER"
$DOCKER ps >/dev/null || true
echo

# ---- Reference ----
REF=$(find ref_mito*/ -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | head -n1 || true)
[[ -n "${REF:-}" ]] || { echo "ERROR: no mito FASTA under ref_mito*/"; exit 1; }
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

CONTIG=$(awk 'NR==1{print $1}' "${REF}.fai")
LEN=$(awk   'NR==1{print $2}' "${REF}.fai")
GENOME="$(dirname "$REF")/mt.genome"
awk 'NR==1{print $1"\t"$2}' "${REF}.fai" > "$GENOME"

PREFIX=$(basename "$PWD")

echo "REF    : $REF"
echo "CONTIG : $CONTIG"
echo "LEN    : $LEN"
echo "GENOME : $GENOME"
echo "PREFIX : $PREFIX"
echo

# ---- BAM list ----
printf '%s\n' bam_rg/*_with_RG.bam 2>/dev/null | sort -V > bam_list.txt || true
sed -i '/\*_/d' bam_list.txt || true
[[ -s bam_list.txt ]] || { echo "ERROR: no *_with_RG.bam in bam_rg/"; exit 1; }

echo "Indexing BAMs if needed..."
while read -r B; do
  [[ -s "${B}.bai" && "${B}.bai" -nt "$B" ]] || samtools index -b "$B"
done < bam_list.txt
echo "BAMs: $(wc -l < bam_list.txt)"
echo

BAM_ARGS=$(awk '{printf "-b /data/%s ", $0}' bam_list.txt)

# ---- Windows ----
WINDOWS="windows.all"
awk -v c="$CONTIG" -v L="$LEN" -v w="$WSIZE" \
  'BEGIN{for(s=1;s<=L;s+=w){e=s+w-1;if(e>L)e=L;print c":"s"-"e}}' > "$WINDOWS"
echo "Windows: $WINDOWS (n=$(wc -l < "$WINDOWS"))"

mkdir -p "$OUT" "$REPORT_DIR"

$DOCKER_RUN -v "$PWD":/data -w /data "$DOCKER_IMAGE" version || true
echo

# ---- MAIN LOOP ----
i=0
while read -r REG; do
  ((i+=1))
  TAG="${PREFIX}_${REG//[:\-]/_}"
  RAW="$OUT/$TAG.raw.vcf.gz"
  PATCH="$OUT/$TAG.gtdpadroao.vcf.gz"
  NORM="$OUT/$TAG.normalise.vcf.gz"

  echo "[$(printf %03d "$i")] $REG -> $TAG"

  # A) FreeBayes
  if [[ ! -s "$RAW" ]]; then
    echo "   freebayes ..."
    $DOCKER_SHELL -v "$PWD":/data -w /data "$DOCKER_IMAGE" -lc "
      set -euo pipefail
      freebayes -f /data/${REF#./} -r $REG --limit-coverage ${LIMIT_COV} $BAM_ARGS \
        --min-mapping-quality ${MIN_MQ} --min-base-quality ${MIN_BQ} \
      | bgzip -c > /data/$RAW
      tabix -f -p vcf /data/$RAW
    " || { echo "   WARN: FreeBayes failed @ $REG — skipping"; continue; }
  fi

# B) Ensure/add GT/DP/AD/RO/AO
  if [[ ! -s "$PATCH" ]]; then
    echo "   ensure/add GT/DP/AD/RO/AO ..."
    zcat "$RAW" | awk -v OFS='\t' '
      BEGIN{FS="\t"; haveGT=haveDP=haveAD=haveRO=haveAO=0}
      /^##FORMAT=<ID=GT,/ {haveGT=1}
      /^##FORMAT=<ID=DP,/ {haveDP=1}
      /^##FORMAT=<ID=AD,/ {haveAD=1}
      /^##FORMAT=<ID=RO,/ {haveRO=1}
      /^##FORMAT=<ID=AO,/ {haveAO=1}

      /^#CHROM/ {
        if(!haveGT) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        if(!haveDP) print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
        if(!haveAD) print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Ref,Alt depths\">"
        if(!haveRO) print "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observations\">"
        if(!haveAO) print "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations\">"
        print; next
      }
      /^##/ {print; next}

      {
        altN = ($5=="." ? 0 : split($5,al,",")); if(altN<1) altN=1

        n = split($9,F,":"); idxGT=idxDP=idxAD=idxRO=idxAO=0
        for(i=1;i<=n;i++){
          if(F[i]=="GT") idxGT=i; else if(F[i]=="DP") idxDP=i; else if(F[i]=="AD") idxAD=i;
          else if(F[i]=="RO") idxRO=i; else if(F[i]=="AO") idxAO=i;
        }
        if(!idxGT){F[++n]="GT"; idxGT=n}
        if(!idxDP){F[++n]="DP"; idxDP=n}
        if(!idxAD){F[++n]="AD"; idxAD=n}
        if(!idxRO){F[++n]="RO"; idxRO=n}
        if(!idxAO){F[++n]="AO"; idxAO=n}

        fmt=F[1]; for(i=2;i<=n;i++) fmt=fmt":"F[i]; $9=fmt

        for(c=10;c<=NF;c++){
          m=split($c,f,":"); for(k=m+1;k<=n;k++) f[k]=""
          # AD
          adN=split(f[idxAD],ad,","); for(t=1;t<=adN;t++){ if(ad[t]==""||ad[t]==".") ad[t]=0; ad[t]+=0 }
          ro=0; if(adN>=1) ro=ad[1]+0; else if(f[idxRO]!=""&&f[idxRO]!=".") ro=f[idxRO]+0
          split("",ao); aoN=0
          if(adN>=2){
            for(t=1;t<=altN;t++) ao[t]=((t+1)<=adN ? ad[t+1]+0 : 0); aoN=altN
          } else {
            aoN=split(f[idxAO],ao,",")
            for(t=1;t<=altN;t++){ if(t>aoN||ao[t]==""||ao[t]==".") ao[t]=0; else ao[t]+=0 }
            aoN=altN
          }
          # sanitize AD vector to 1+altN
          if(adN < (1+altN)){
            adstr=ro; for(t=1;t<=altN;t++) adstr=adstr","ao[t]
          } else {
            adstr=ad[1]; for(t=1;t<=altN;t++) adstr=adstr","((t+1)<=adN?ad[t+1]:0)
          }
          # DP = sum(AD)
          dp=ro; for(t=1;t<=altN;t++) dp+=ao[t]
          # GT default
          gt=f[idxGT]; if(gt==""||gt==".") gt="0/0"
          # AO string
          aostr=ao[1]; for(t=2;t<=altN;t++) aostr=aostr","ao[t]

          f[idxGT]=gt; f[idxAD]=adstr; f[idxRO]=ro; f[idxAO]=aostr; f[idxDP]=dp
          s=f[1]; for(k=2;k<=n;k++) s=s":"f[k]; $c=s
        }
        print
      }
    ' | bgzip -c > "$PATCH"
    tabix -f -p vcf "$PATCH"
  fi

  # C) mity normalise
  if [[ ! -s "$NORM" ]]; then
    echo "   mity normalise ..."
    $DOCKER_RUN -v "$PWD":/data -w /data "$DOCKER_IMAGE" normalise \
      --custom-reference-fasta  "/data/${REF#./}" \
      --custom-reference-genome "/data/${GENOME#./}" \
      --prefix "$TAG" \
      --output-dir "/data/$OUT" \
      "/data/$PATCH"
    if [[ -f "$OUT/${TAG}.mity.normalise.vcf.gz" ]]; then
      mv -f "$OUT/${TAG}.mity.normalise.vcf.gz" "$NORM"
      [[ -f "$OUT/${TAG}.mity.normalise.vcf.gz.tbi" ]] && mv -f "$OUT/${TAG}.mity.normalise.vcf.gz.tbi" "${NORM}.tbi"
    fi
  fi
done < "$WINDOWS"

# ---- Merge ----
echo
echo "==> Concatenating normalised tiles ..."
find "$OUT" -maxdepth 1 -type f -name "${PREFIX}_${CONTIG}_*.normalise.vcf.gz" | LC_ALL=C sort -V > "$OUT/chunks.norm.list"
[[ -s "$OUT/chunks.norm.list" ]] || { echo "ERROR: no normalised tiles to merge"; exit 1; }

MERGED="$OUT/${PREFIX}.normalise.vcf.gz"
bcftools concat -a -f "$OUT/chunks.norm.list" -Oz -o "$MERGED"
tabix -f "$MERGED"

FINAL="$OUT/${PREFIX}.normalise.vaf0.vcf.gz"
echo "==> Zero-filling missing VAF ..."
zcat "$MERGED" | awk -v OFS='\t' '
  /^##/||/^#CHROM/{print; next}
  { print }' | bgzip -c > "$FINAL"
tabix -f "$FINAL"

# ---- Convert GFF3 to BED and fix contig name ----
GFF3=$(find "$(dirname "$REF")" -maxdepth 1 -type f -name "*.gff3.gz" | head -n1 || true)

if [[ -n "$GFF3" ]]; then
  BED="$(dirname "$REF")/$(basename "${GFF3%.gff3.gz}").bed"

  if [[ "$GFF3" -nt "${BED}.gz" ]] || [[ ! -s "${BED}.gz" ]]; then
    echo "==> Converting GFF3 to BED (using VCF contig: $CONTIG)..."

    zcat "$GFF3" | awk -v contig="$CONTIG" -F'\t' '
      BEGIN { OFS = "\t" }
      !/^#/ && $3 != "source" && $3 != "gene" {
        gene = ""; product = ""; name = ""; id = ""

        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
          attr = attrs[i]
          gsub(/^ +| +$/, "", attr)
          if (attr ~ /^gene=/) { gsub(/^gene=/, "", attr); gene = attr }
          else if (attr ~ /^Name=/) { gsub(/^Name=/, "", attr); name = attr }
          else if (attr ~ /^product=/) { gsub(/^product=/, "", attr); product = attr }
          else if (attr ~ /^ID=/) { gsub(/^ID=/, "", attr); id = attr }
        }

        gene_name = (gene != "" ? gene : (name != "" ? name : (id != "" ? id : $3)))
        if (product == "") product = $3

        # Use VCF contig name!
        print contig, $4-1, $5, gene_name, ".", $7, gene_name, product
      }
    ' | sort -k1,1 -k2,2n | bgzip -c > "${BED}.gz"

    tabix -p bed "${BED}.gz"
    echo "   ✓ Created: ${BED}.gz (contig: $CONTIG)"
  fi
fi

# ---- Annotate VCF with vcfanno (on host) ----
ANNOTATED_VCF="$FINAL"

if [[ -n "$GFF3" ]] && [[ -s "${BED}.gz" ]] && command -v vcfanno >/dev/null 2>&1; then
  echo "==> Annotating VCF with gene information..."

  TOML="${BED%.bed}_vcfanno.toml"
  cat > "$TOML" <<EOF
[[annotation]]
file = "${BED}.gz"
columns = [7, 8]
names = ["gene", "product"]
ops = ["uniq", "uniq"]
EOF

  ANNOTATED_VCF="$OUT/${PREFIX}.normalise.annotated.vcf.gz"

  # Run vcfanno on host (not in Docker)
  vcfanno "$TOML" "$FINAL" 2>/dev/null | bgzip -c > "$ANNOTATED_VCF"
  tabix -f "$ANNOTATED_VCF"

  # Verify
  ANN_COUNT=$(zcat "$ANNOTATED_VCF" | grep -v "^#" | grep -c "gene=" || echo 0)
  TOTAL=$(zcat "$ANNOTATED_VCF" | grep -v "^#" | wc -l)
  echo "   ✓ Annotated: $ANN_COUNT / $TOTAL variants"

  if [[ $ANN_COUNT -eq 0 ]]; then
    echo "   ⚠️  WARNING: No variants were annotated - using unannotated VCF"
    ANNOTATED_VCF="$FINAL"
  fi
else
  if [[ -z "$GFF3" ]]; then
    echo "==> No GFF3 annotation available"
  elif ! command -v vcfanno >/dev/null 2>&1; then
    echo "==> vcfanno not installed - skipping annotation"
    echo "   Install: sudo apt-get install vcfanno"
  fi
fi

# ---- Generate MITY report ----
echo
echo "==> Generating MITY report..."

$DOCKER_RUN -v "$PWD":/data -w /data "$DOCKER_IMAGE" report \
  --prefix "${PREFIX}_all" \
  --output-dir "/data/$REPORT_DIR" \
  "/data/${ANNOTATED_VCF#./}"

# ---- Extract CSV from XLSX and populate gene columns ----
echo
echo "==> Extracting CSV from XLSX and adding gene annotations..."

if ! command -v ssconvert >/dev/null 2>&1; then
  echo "WARNING: ssconvert not found. Install with: sudo apt-get install gnumeric"
  echo "Skipping CSV extraction."
else
  for xlsx_report in "$REPORT_DIR"/*.mity.report.xlsx; do
    [[ -f "$xlsx_report" ]] || continue

    csv_report="${xlsx_report%.xlsx}.csv"

    echo "   Converting: $(basename "$xlsx_report")"
    ssconvert "$xlsx_report" "$csv_report" 2>/dev/null

    if [[ ! -f "$csv_report" ]]; then
      echo "   ERROR: Failed to extract CSV"
      continue
    fi

    # Extract gene annotations from INFO field to GENE/LOCUS columns
    echo "   Populating gene columns from INFO field..."

    awk -F',' 'BEGIN{OFS=","}
      NR==1 {
        print
        next
      }
      {
        info = $(NF-1)
        gene = "."
        product = "."

        # Extract gene
        if (info ~ /gene=/) {
          tmp = info
          sub(/.*gene=/, "", tmp)
          sub(/[;,].*/, "", tmp)
          gene = tmp
        }

        # Extract product
        if (info ~ /product=/) {
          tmp = info
          sub(/.*product=/, "", tmp)
          sub(/[;,].*/, "", tmp)
          product = tmp
        }

        $3 = gene
        $4 = product
        print
      }
    ' "$csv_report" > "${csv_report}.tmp"

    mv "${csv_report}.tmp" "$csv_report"

    # Count annotated variants
    ANNOTATED=$(awk -F',' 'NR>1 && $3!="." && $3!=""' "$csv_report" | wc -l)
    TOTAL=$(awk 'NR>1' "$csv_report" | wc -l)
    echo "   ✓ Populated gene info: $ANNOTATED / $TOTAL variants"
  done
fi
    # Count annotated variants
    ANNOTATED=$(awk -F',' 'NR>1 && $3!="." && $3!=""' "$csv_report" | wc -l)
    TOTAL=$(awk 'NR>1' "$csv_report" | wc -l)
    echo "   ✓ Populated gene info: $ANNOTATED / $TOTAL variants"
  done
fi

echo
echo "==> Done."
echo " - Tiles dir     : $OUT"
echo " - Merged VCF    : $MERGED"
echo " - VAF0 VCF      : $FINAL"
echo " - Report (xlsx) : $REPORT_DIR/${PREFIX}_all.mity.report.xlsx"
if command -v ssconvert >/dev/null 2>&1; then
  echo " - Annotated CSV : $REPORT_DIR/${PREFIX}_all.mity.report.csv"
fi
echo

# ---- Variant Summary ----

echo "==> Variant Summary:"
TOTAL_VARS=$(zcat "$FINAL" | grep -v "^#" | wc -l || echo 0)
echo "Total variants: $TOTAL_VARS"

if [[ $TOTAL_VARS -gt 0 ]]; then
  echo
  echo "First 5 variants:"
  echo "CHROM           POS    REF  ALT"
  echo "--------------------------------------"
  zcat "$FINAL" | grep -v "^#" | head -5 | awk '{printf "%-15s %-6s %-4s %s\n", $1, $2, $4, $5}' || true

  # Show CSV if exists
  csv_report="$REPORT_DIR/${PREFIX}_all.mity.report.csv"
  if [[ -f "$csv_report" ]]; then
    echo
    echo "Sample from CSV report (first 3 variants):"
    head -4 "$csv_report" | tail -3 | cut -d',' -f1-5 | column -t -s','
  fi
fi

echo
echo "Pipeline complete! ✅"
echo "Check report: $REPORT_DIR/${PREFIX}_all.mity.report.xlsx"
if command -v ssconvert >/dev/null 2>&1; then
  echo "CSV version: $REPORT_DIR/${PREFIX}_all.mity.report.csv"
fi
