#!/usr/bin/env bash
# Modified to work from parent directory
set -euo pipefail

SPECIES=$1
[[ -z "$SPECIES" ]] && { echo "Usage: $0 <species_directory>"; exit 1; }

# Source config
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/config/pipeline.conf"

cd "$SPECIES"
mkdir -p bam_rg

for bam_file in *.sorted.bam; do
  [[ -e $bam_file ]] || continue
  
  base_name=$(basename "$bam_file" .sorted.bam)
  output_bam="bam_rg/${base_name}_with_RG.bam"
  
  [[ -f $output_bam ]] && continue
  
  echo "  Adding RG: $bam_file"
  java -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$bam_file" \
    O="$output_bam" \
    RGID="$base_name" \
    RGLB="$RG_LIBRARY" \
    RGPL="$RG_PLATFORM" \
    RGPU="${RG_PLATFORM_UNIT}_${base_name}" \
    RGSM="$base_name" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT
done

cd ..
