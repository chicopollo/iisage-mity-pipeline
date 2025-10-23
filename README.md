# Mitochondrial Variant Calling Pipeline

Automated mitochondrial genome variant calling pipeline for non-model species using BWA, MITY, FreeBayes, and custom annotations.
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

This pipeline performs end-to-end analysis from FASTQ files to annotated variant reports:

1. **Alignment**: BWA-MEM alignment to mitochondrial reference
2. **Read Groups**: Add sample metadata with Picard
3. **Annotation Prep**: Convert GenBank to GFF3 format
4. **Variant Calling**: FreeBayes variant calling with MITY normalization, filtering and reporting
5. **Annotation**: Gene-level variant annotation with product descriptions

## Quick Start
```bash
# 1. Clone repository
git clone https://github.com/chicopollo/iisage-mity-pipeline.git
cd iisage-mity-pipeline

# 2. Setup (downloads Picard, checks dependencies)
bash setup.sh

# 3. Organize your data (see Directory Structure below)

# 4. Run pipeline
./run_mity_pipeline.sh
```

## Requirements

### Software Dependencies
- **BWA** (≥0.7.17)
- **samtools** (≥1.10)
- **bcftools** (≥1.10)
- **tabix/bgzip** (htslib)
- **Docker** (for MITY container)
- **Java** (≥8, for Picard)
- **seqtk** (for GenBank conversion)
- **ssconvert** (gnumeric, for Excel→CSV)
- **Python 3** with Biopython (in Docker)

### Installation

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y bwa samtools bcftools tabix docker.io openjdk-11-jre seqtk gnumeric
```

**Conda:**
```bash
conda install -c bioconda bwa samtools bcftools htslib seqtk
```

**Check installation:**
```bash
bash check_dependencies.sh
```

## Directory Structure

### Input Structure
```
mito-variant-pipeline/
├── run_mity_pipeline.sh
├── config/
│   └── pipeline.conf
├── scripts/
│   ├── 01_bwa_align.sh
│   ├── 02_add_readgroups.sh
│   ├── 03_genbank_to_gff3.sh
│   └── 04_mity_pipeline.sh
└── Species_name/                    # One or more species
    ├── SampleA/
    │   ├── SampleA_R1_001.fastq.gz
    │   └── SampleA_R2_001.fastq.gz
    ├── SampleB/
    │   ├── SampleB_R1_001.fastq.gz
    │   └── SampleB_R2_001.fastq.gz
    └── ref_mito_speciesname/
        ├── genome.fa               # Mitochondrial reference fasta
        └── genome.gb               # GenBank annotation
```

### Output Structure
```
Species_name/
├── Species_name_SampleA.sorted.bam           # Aligned BAMs
├── Species_name_SampleB.sorted.bam
├── bam_rg/                                   # BAMs with read groups
│   ├── Species_name_SampleA_with_RG.bam
│   └── Species_name_SampleB_with_RG.bam
├── ref_mito_speciesname/
│   ├── genome.gff3.gz                        # Converted annotations
│   └── genome.bed.gz                         # BED format for vcfanno
├── mity_chunked/                             # Intermediate VCFs
│   └── Species_name.normalise.vaf0.vcf.gz    # Final merged VCF
└── mity_output/                              # Final reports
    ├── Species_name_all.mity.report.xlsx     # Excel report
    ├── Species_name_all.mity.report.csv      # CSV report
    └── Species_name_all.mity.report_annotated.csv  # ⭐ Annotated CSV
```

## Usage

### Run Entire Pipeline
```bash
# All species, all steps
./run_mity_pipeline.sh

# Specific species
./run_mity_pipeline.sh --species Chrysemys_picta

# Multiple species
./run_mity_pipeline.sh Chrysemys_picta Desmodus_rotundus
```

### Run Specific Steps
```bash
# Only alignment (step 1)
./run_mity_pipeline.sh --step 1

# Only variant calling (step 4)
./run_mity_pipeline.sh --step 4

# Steps 3-4 (annotation prep and variant calling)
./run_mity_pipeline.sh --from-step 3
```

### Other Options
```bash
# List available species
./run_mity_pipeline.sh --list

# Force re-run (skip completion checks)
./run_mity_pipeline.sh --force

# Show help
./run_mity_pipeline.sh --help
```

## Configuration

Edit `config/pipeline.conf` to customize parameters:

```bash
# MITY parameters
MITY_WINDOW_SIZE=300        # Variant calling window size (bp)
MITY_COVERAGE_CAP=8000      # Maximum coverage for variant calling
MITY_MIN_MQ=30              # Minimum mapping quality
MITY_MIN_BQ=24              # Minimum base quality

# BWA parameters
BWA_THREADS=8               # Number of threads for alignment

# Read group defaults
RG_LIBRARY="lib1"
RG_PLATFORM="illumina"
```

## Output Files

### Key Outputs

**Annotated CSV Report** (`*_annotated.csv`):
- Columns 3-4 contain gene names and product descriptions
- All MITY quality metrics included
- Ready for downstream analysis in R

**VCF File** (`*.normalise.vaf0.vcf.gz`):
- MITY-normalized variants
- Gene annotations in INFO field
- Compatible with standard VCF tools

## Pipeline Steps Explained

### Step 1: BWA Alignment
- Converts GenBank to FASTA if needed
- Indexes reference with BWA
- Aligns paired-end reads
- Sorts and indexes BAM files

### Step 2: Add Read Groups
- Uses Picard AddOrReplaceReadGroups
- Adds sample metadata required by MITY
- Creates indexed BAMs in `bam_rg/`

### Step 3: Annotation Preparation
- Converts GenBank (.gb) to GFF3 format
- Extracts gene, CDS, tRNA, rRNA features
- Converts GFF3 to BED for vcfanno
- Compresses and indexes annotations

### Step 4: Variant Calling
- Windowed FreeBayes calling (300bp windows)
- MITY normalization (handles circular genomes)
- Variant quality filtering
- Gene annotation using vcfanno
- Generates Excel and annotated CSV reports

## Troubleshooting

### Common Issues

**1. "No species directories found"**
```bash
# Check directory structure
./run_mity_pipeline.sh --list

# Ensure directory naming: Species_name/ref_mito_*/
```

**2. "Docker permission denied"**
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Log out and back in
```

**3. "ssconvert not found"**
```bash
# Install gnumeric
sudo apt-get install gnumeric
```

**4. Empty GFF3 files**
```bash
# Regenerate from GenBank
bash scripts/03_genbank_to_gff3.sh
```

**5. MITY container issues**
```bash
# Pull latest image
docker pull drmjc/mity:2.0.0
```

## Advanced Usage

### Resume Failed Pipeline
The pipeline automatically detects completed steps:
```bash
# Will skip completed steps for each species
./run_mity_pipeline.sh
```

### Custom Configuration
```bash
# Use custom config file
CONFIG_FILE=/path/to/custom.conf ./run_mity_pipeline.sh
```

### Manual Step Execution
```bash
# Run individual scripts (from parent directory)
bash scripts/01_bwa_align.sh Chrysemys_picta
bash scripts/02_add_readgroups.sh Chrysemys_picta

# Or from species directory
cd Chrysemys_picta
bash ../scripts/04_mity_pipeline.sh
```

## Citation

If you use this pipeline, please cite:

- **MITY**: Puttick C, et al. (2024) "mity: A highly sensitive mitochondrial variant analysis pipeline for whole genome sequencing data." Journal of bionformatics and systems biology.
- **FreeBayes**: Garrison E & Marth G. (2012) "Haplotype-based variant detection from short-read sequencing." arXiv
- **BWA**: Li H, Durbin R. (2009) "Fast and accurate short read alignment with Burrows-Wheeler transform." Bioinformatics

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Support

- **Issues**: https://github.com/chicopllo/iisage-mity-pipeline/issues
- **Contact**: decenalo@msu.edu
- **Lab**: Bronikowski Lab, IISAGE

## License

This project is licensed under the Creative Commons Attribution 4.0 International License - see the [LICENSE](LICENSE) file for details

## Authors

**L. Paul Decena-Segarra** 

**Eric Randolph**

Bronikowski Lab - 2025

Michigan State University - IISAGE

## Acknowledgments

- MITY development team
- FreeBayes development team
- Bronikowski Lab members
- IISAGE members
