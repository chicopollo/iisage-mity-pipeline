#!/usr/bin/env bash
# Modified to accept species as argument
set -euo pipefail

SPECIES=$1
[[ -z "$SPECIES" ]] && { echo "Usage: $0 <species_directory>"; exit 1; }

cd "$SPECIES"

# Find reference
reference_fa=$(find ref_mito*/ -type f \( -name "*.fa" -o -name "*.fasta" \) -print -quit)

if [[ -z $reference_fa ]]; then
  reference_gb=$(find ref_mito*/ -type f -name "*.gb" -print -quit)
  if [[ ! -f $reference_gb ]]; then
    echo "ERROR: No reference found for $SPECIES"
    exit 1
  fi
  reference_fa=${reference_gb%.gb}.fasta
  seqtk seq -A "$reference_gb" > "$reference_fa"
fi

# Index if needed
[[ ! -f ${reference_fa}.bwt ]] && bwa index "$reference_fa"

# Process FASTQ pairs
for fq1 in */*_R1_*fastq.gz; do
  [[ -e $fq1 ]] || continue
  
  fq2=${fq1/_R1_/_R2_}
  [[ ! -f $fq2 ]] && continue
  
  sample_base=$(basename "$fq1" _R1_001.fastq.gz)
  bam_sorted="${SPECIES}_${sample_base}.sorted.bam"
  
  [[ -f $bam_sorted ]] && continue
  
  echo "  Aligning: $sample_base"
  bwa mem "$reference_fa" <(gunzip -c "$fq1") <(gunzip -c "$fq2") \
    | samtools view -bS \
    | samtools sort -o "$bam_sorted"
  
  samtools index "$bam_sorted"
done

cd ..
