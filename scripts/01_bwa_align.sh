#!/usr/bin/env bash
# Modified to accept species as argument
# Fixed to manage .gb files
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
  
  echo "Converting GenBank to FASTA..."
  # Use AWK instead of seqtk for reliable conversion
  awk '
    /^LOCUS/ {acc=$2}
    /^ORIGIN/ {inseq=1; next}
    /^\/\// {exit}
    inseq {
      gsub(/[^acgtACGT]/, "")
      seq = seq $0
    }
    END {
      # Use accession as contig name
      print ">" acc
      for(i=1; i<=length(seq); i+=80) {
        print substr(seq, i, 80)
      }
    }
  ' "$reference_gb" > "$reference_fa"
  
  # Verify FASTA was created
  if [[ ! -s "$reference_fa" ]] || ! grep -q "^>" "$reference_fa"; then
    echo "ERROR: Failed to create valid FASTA from $reference_gb"
    exit 1
  fi
  
  echo "Created: $reference_fa"
fi

# Index if needed
if [[ ! -f ${reference_fa}.bwt ]]; then
  echo "Indexing reference with BWA..."
  bwa index "$reference_fa" || { echo "ERROR: BWA indexing failed"; exit 1; }
fi

# Index with samtools
if [[ ! -f ${reference_fa}.fai ]]; then
  samtools faidx "$reference_fa" || { echo "ERROR: samtools faidx failed"; exit 1; }
fi

# Process FASTQ pairs
for fq1 in */*_R1_*.fastq.gz; do
  [[ -e $fq1 ]] || continue
  
  fq2=${fq1/_R1_/_R2_}
  [[ ! -f $fq2 ]] && { echo "WARNING: Missing R2 for $fq1"; continue; }
  
  sample_base=$(basename "$fq1" _R1_001.fastq.gz)
  bam_sorted="${SPECIES}_${sample_base}.sorted.bam"
  
  if [[ -f $bam_sorted ]]; then
    echo "  Skipping $sample_base (BAM exists)"
    continue
  fi
  
  echo "  Aligning: $sample_base"
  
  # Use temp files instead of process substitution for better error handling
  bwa mem -t 4 "$reference_fa" "$fq1" "$fq2" 2>/dev/null \
    | samtools view -bS - 2>/dev/null \
    | samtools sort -o "$bam_sorted" - || {
      echo "ERROR: Alignment failed for $sample_base"
      rm -f "$bam_sorted"
      exit 1
    }
  
  samtools index "$bam_sorted" || {
    echo "ERROR: Indexing failed for $bam_sorted"
    exit 1
  }
  
  echo "  ✓ Created: $bam_sorted"
done

echo "✓ Alignment complete for $SPECIES"
cd ..
