#!/usr/bin/env bash
# gb_to_gff3_v4.4.sh - Final working version
set -euo pipefail

DOCKER_IMAGE="biopython/biopython:latest"

echo "==> Starting GenBank to GFF3 conversion"
echo

for gb in $(find . -type f \( -name "*.gb" -o -name "*.gbk" -o -name "*.gbff" \) | grep "ref_mito_" | sort); do
  
  # Skip backup files
  [[ "$(basename "$gb")" == *"_fixed"* ]] && continue
  
  gff3="${gb%.*}.gff3"
  species=$(basename "$(dirname "$gb")")
  
  echo "[$species]"
  echo "  Converting: $(basename "$gb") → $(basename "$gff3")"
  
  # Remove old files
  rm -f "$gff3" "${gff3}.gz" "${gff3}.gz.tbi"
  
  # Get absolute paths
  GB_ABS=$(cd "$(dirname "$gb")" && pwd)/$(basename "$gb")
  GFF3_ABS=$(cd "$(dirname "$gb")" && pwd)/$(basename "$gff3")
  WORKDIR=$(dirname "$GB_ABS")
  
  # Run conversion with proper mounting
  sudo docker run --rm -v "$WORKDIR":/work -w /work "$DOCKER_IMAGE" python3 -c '
from Bio import SeqIO
import sys

try:
    gb_file = "'"$(basename "$GB_ABS")"'"
    gff_file = "'"$(basename "$GFF3_ABS")"'"
    
    with open(gff_file, "w") as out:
        out.write("##gff-version 3\n")
        
        total_features = 0
        for record in SeqIO.parse(gb_file, "genbank"):
            out.write("##sequence-region {} 1 {}\n".format(record.id, len(record.seq)))
            
            for feature in record.features:
                # Only relevant mitochondrial features
                if feature.type not in ["source", "gene", "CDS", "rRNA", "tRNA", 
                                       "misc_feature", "D-loop", "repeat_region", 
                                       "ncRNA", "misc_RNA", "D_loop"]:
                    continue
                
                # Get coordinates (GFF3 is 1-based, inclusive)
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                
                # Strand
                if feature.location.strand == 1:
                    strand = "+"
                elif feature.location.strand == -1:
                    strand = "-"
                else:
                    strand = "."
                
                # Build attributes
                attrs = []
                
                # ID (mandatory)
                if "locus_tag" in feature.qualifiers:
                    feat_id = feature.qualifiers["locus_tag"][0]
                elif "gene" in feature.qualifiers:
                    feat_id = feature.qualifiers["gene"][0]
                else:
                    feat_id = "{}_{}".format(feature.type, start)
                
                # Clean ID
                feat_id = feat_id.replace(" ", "_").replace(";", "_").replace(",", "_")
                attrs.append("ID=" + feat_id)
                
                # Name and gene
                if "gene" in feature.qualifiers:
                    gene = feature.qualifiers["gene"][0].replace(";", ",")
                    attrs.append("Name=" + gene)
                    attrs.append("gene=" + gene)
                
                # Product
                if "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0].replace(";", ",")
                    attrs.append("product=" + product)
                
                # Note
                if "note" in feature.qualifiers:
                    note = feature.qualifiers["note"][0].replace(";", ",")
                    attrs.append("Note=" + note)
                
                # Write GFF3 line
                out.write("{}\tGenBank\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(
                    record.id, feature.type, start, end, strand, ";".join(attrs)))
                total_features += 1
        
        print("SUCCESS: Extracted {} features from {} ({}bp)".format(
            total_features, record.id, len(record.seq)), file=sys.stderr)
        
        if total_features == 0:
            sys.exit(1)

except Exception as e:
    print("ERROR: " + str(e), file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)
'
  
  # Check result
  if [[ ! -s "$gff3" ]]; then
    echo "  ❌ FAILED - no output created"
    continue
  fi
  
  # Count features
  feature_count=$(grep -v "^#" "$gff3" | wc -l)
  
  if [[ $feature_count -eq 0 ]]; then
    echo "  ❌ FAILED - GFF3 is empty"
    rm -f "$gff3"
    continue
  fi
  
  # Compress and index
  bgzip -f "$gff3"
  tabix -f -p gff "${gff3}.gz"
  
  echo "  ✅ SUCCESS: $feature_count features"
  
  # Show sample
  echo "  Sample features:"
  zcat "${gff3}.gz" | grep -v "^#" | head -3 | awk -F'\t' '{
    split($9, a, ";")
    for(i in a) {
      if(a[i] ~ /^product=/) {
        sub(/^product=/, "", a[i])
        product = a[i]
        break
      }
    }
    printf "    %s (%d-%d): %s\n", $3, $4, $5, product
  }'
  echo
  
done

echo
echo "=== CONVERSION SUMMARY ==="
printf "%-40s %10s\n" "Species" "Features"
printf "%s\n" "------------------------------------------------------------"
for gff3 in $(find . -type f -name "*.gff3.gz" | grep "ref_mito_" | sort); do
  count=$(zcat "$gff3" | grep -v "^#" | wc -l)
  species=$(basename "$(dirname "$gff3")")
  printf "%-40s %10d\n" "$species" "$count"
done
echo
echo "All GFF3 files created and indexed!"
echo
echo "Next: Run your MITY pipeline to generate annotated reports."
