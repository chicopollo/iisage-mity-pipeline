#!/usr/bin/env bash
# Modified header to source config:

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/pipeline.conf"

# ===============================================================
# mity2_species_runall_v7.1.sh
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
# Author: Paul Ten (Bronikowski Lab, IISAGE)
# Version: 7.1 (2025-10)
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

echo "==> MITY species pipeline (v7.1)"
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

  # B) Sanitize FORMAT fields
  if [[ ! -s "$PATCH" ]]; then
    echo "   ensure/add GT/DP/AD/RO/AO (sanitize NoneType) ..."
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
        if(NF < 10) {print; next}

        n = split($9,F,":"); idxRO=idxDP=idxAO=0
        for(i=1;i<=n;i++){
          if(F[i]=="RO") idxRO=i
          else if(F[i]=="DP") idxDP=i
          else if(F[i]=="AO") idxAO=i
        }

        for(c=10;c<=NF;c++){
          m=split($c,f,":")
          if(m<n){ for(j=m+1;j<=n;j++) f[j]="." }
          if(idxRO && (f[idxRO]==""||f[idxRO]=="."||f[idxRO]=="None")) f[idxRO]=0
          if(idxDP && (f[idxDP]==""||f[idxDP]=="."||f[idxDP]=="None")) f[idxDP]=0
          if(idxAO && (f[idxAO]==""||f[idxAO]=="."||f[idxAO]=="None")) f[idxAO]="0"
          s=f[1]; for(j=2;j<=n;j++) s=s":"f[j]; $c=s
        }
        print
      }' | bgzip -c > "$PATCH"
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

# ---- Convert GFF3 to BED for vcfanno ----
GFF3=$(find "$(dirname "$REF")" -maxdepth 1 -type f -name "*.gff3.gz" | head -n1 || true)

if [[ -n "$GFF3" ]]; then
  BED="$(dirname "$REF")/$(basename "${GFF3%.gff3.gz}").bed"
  
  # Force regeneration if GFF3 is newer
  if [[ "$GFF3" -nt "${BED}.gz" ]] || [[ ! -s "${BED}.gz" ]]; then
    echo "==> Converting GFF3 to BED for vcfanno..."
    zcat "$GFF3" | awk -F'\t' '
      BEGIN {
        OFS = "\t"
      }
      
      # Skip comments, source, and gene features (gene features typically lack useful annotations)
      !/^#/ && $3 != "source" && $3 != "gene" {
        gene = ""
        product = ""
        name = ""
        id = ""
        
        # Parse attributes
        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
          attr = attrs[i]
          gsub(/^ +| +$/, "", attr)
          
          if (attr ~ /^gene=/) {
            gsub(/^gene=/, "", attr)
            gene = attr
          }
          else if (attr ~ /^Name=/) {
            gsub(/^Name=/, "", attr)
            name = attr
          }
          else if (attr ~ /^product=/) {
            gsub(/^product=/, "", attr)
            product = attr
          }
          else if (attr ~ /^ID=/) {
            gsub(/^ID=/, "", attr)
            id = attr
          }
        }
        
        # Prioritize: gene > Name > ID > feature_type
        gene_name = gene
        if (gene_name == "") gene_name = name
        if (gene_name == "") gene_name = id
        if (gene_name == "") gene_name = $3
        
        # Use product if available, otherwise feature type
        if (product == "") product = $3
        
        # BED format: chr, start(0-based), end, name, score, strand, gene, product
        print $1, $4-1, $5, gene_name, ".", $7, gene_name, product
      }
    ' | sort -k1,1 -k2,2n | bgzip -c > "${BED}.gz"
    
    tabix -p bed "${BED}.gz"
    echo "   Created: ${BED}.gz"
    
    # Show statistics
    FEAT_COUNT=$(zcat "${BED}.gz" | wc -l)
    echo "   Features: $FEAT_COUNT"
    echo "   Sample annotations:"
    zcat "${BED}.gz" | head -5 | awk '{printf "     %s:%d-%d  %s (%s)\n", $1, $2, $3, $7, $8}'
  else
    echo "==> BED annotation exists: ${BED}.gz"
  fi
fi

# ---- Annotation-aware report ----
echo
echo "==> Generating mity report ..."

if [[ -n "$GFF3" ]] && [[ -s "${BED}.gz" ]]; then
  echo "Using BED annotation: ${BED}.gz"
  
  # Create TOML config for vcfanno
  TOML="${BED%.bed}_vcfanno.toml"
  cat > "$TOML" <<EOF
[[annotation]]
file = "/data/${BED#./}.gz"
columns = [7, 8]
names = ["gene", "product"]
ops = ["uniq", "uniq"]
EOF
  
  echo "TOML config: $TOML"
  
  # Run mity report with annotations
  $DOCKER_RUN -v "$PWD":/data -w /data "$DOCKER_IMAGE" report \
    --prefix "${PREFIX}_all" \
    --vcfanno-base-path "/data" \
    --custom-vcfanno-config "/data/${TOML#./}" \
    --output-dir "/data/$REPORT_DIR" \
    "/data/${FINAL#./}" || {
      echo "⚠️  Report with annotations failed, trying without..."
      $DOCKER_RUN -v "$PWD":/data -w /data "$DOCKER_IMAGE" report \
        --prefix "${PREFIX}_all" \
        --output-dir "/data/$REPORT_DIR" \
        "/data/${FINAL#./}"
    }
else
  echo "⚠️  No valid annotation; running report without annotations."
  $DOCKER_RUN -v "$PWD":/data -w /data "$DOCKER_IMAGE" report \
    --prefix "${PREFIX}_all" \
    --output-dir "/data/$REPORT_DIR" \
    "/data/${FINAL#./}"
fi

# ---- Extract CSV from XLSX and add annotations ----
echo
echo "==> Extracting CSV from XLSX reports..."

# Check if ssconvert is available
if ! command -v ssconvert >/dev/null 2>&1; then
  echo "WARNING: ssconvert not found. Install with: sudo apt-get install gnumeric"
  echo "Skipping CSV extraction and annotation."
else
  for xlsx_report in "$REPORT_DIR"/*.mity.report.xlsx; do
    [[ -f "$xlsx_report" ]] || continue
    
    csv_report="${xlsx_report%.xlsx}.csv"
    annotated_csv="${xlsx_report%.xlsx}_annotated.csv"
    
    echo "   Converting: $(basename "$xlsx_report")"
    
    # Convert XLSX to CSV
    ssconvert "$xlsx_report" "$csv_report" 2>/dev/null
    
    if [[ ! -f "$csv_report" ]]; then
      echo "   ERROR: Failed to extract CSV"
      continue
    fi
    
    echo "   Annotating: $(basename "$csv_report")"
    
    # Add annotations from INFO field
    awk -F',' 'BEGIN {OFS=","}
      NR==1 {
        print
        next
      }
      {
        info = $(NF-1)
        gene = "."
        product = "."
        
        if (info ~ /gene=/) {
          tmp = info
          sub(/.*gene=/, "", tmp)
          sub(/;.*/, "", tmp)
          gene = tmp
        }
        
        if (info ~ /product=/) {
          tmp = info
          sub(/.*product=/, "", tmp)
          sub(/;.*/, "", tmp)
          product = tmp
        }
        
        $3 = gene
        $4 = product
        
        print
      }
    ' "$csv_report" > "$annotated_csv"
    
    echo "   Created: $(basename "$annotated_csv")"
    
    ANNOTATED_COUNT=$(awk -F',' 'NR>1 && $3!="." && $3!=""' "$annotated_csv" | wc -l)
    TOTAL_COUNT=$(awk 'NR>1' "$annotated_csv" | wc -l)
    echo "   Annotated: $ANNOTATED_COUNT / $TOTAL_COUNT variants"
    
  done
fi

echo
echo "==> Done."
echo " - Tiles dir     : $OUT"
echo " - Merged VCF    : $MERGED"
echo " - VAF0 VCF      : $FINAL"
echo " - Report (xlsx) : $REPORT_DIR/${PREFIX}_all.mity.report.xlsx"
if command -v ssconvert >/dev/null 2>&1; then
  echo " - Annotated CSV : $REPORT_DIR/${PREFIX}_all.mity.report_annotated.csv"
fi
echo

# ---- Variant Summary ----
echo "==> Variant Summary:"
TOTAL_VARS=$(zcat "$FINAL" | grep -v "^#" | wc -l)
echo "Total variants: $TOTAL_VARS"

if [[ $TOTAL_VARS -gt 0 ]]; then
  echo
  echo "First 5 variants:"
  echo "CHROM           POS    REF  ALT"
  echo "--------------------------------------"
  zcat "$FINAL" | grep -v "^#" | head -5 | awk '{printf "%-15s %-6s %-4s %s\n", $1, $2, $4, $5}'
  
  # Show sample from annotated CSV
  annotated_csv="$REPORT_DIR/${PREFIX}_all.mity.report_annotated.csv"
  if [[ -f "$annotated_csv" ]]; then
    echo
    echo "Sample annotations from CSV (first 3 variants):"
    echo "SAMPLE                         HGVS       GENE       PRODUCT"
    echo "------------------------------------------------------------------------"
    awk -F',' 'NR>1 && NR<=4 {printf "%-30s %-10s %-10s %s\n", substr($1,1,30), substr($2,1,10), substr($3,1,10), substr($4,1,40)}' \
      "$annotated_csv"
  fi
fi

echo
echo "Pipeline complete! ✅"
if command -v ssconvert >/dev/null 2>&1; then
  echo "Check annotated CSV for gene names and products: $REPORT_DIR/${PREFIX}_all.mity.report_annotated.csv"
else
  echo "Install ssconvert (gnumeric) to enable CSV extraction and annotation."
fi
