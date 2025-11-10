# Future Directions: Expanding Heteroplasmy Analysis

## Overview

This pipeline, including the accompanying R analysis scripts, provides a flexible and scalable framework for investigating mitochondrial heteroplasmy across diverse taxa. While the current dataset may have limited statistical power due to sample size, the infrastructure established here opens numerous avenues for expanded research into age-related and sex-specific patterns of mtDNA variation in non-model organisms.

## Expanding Statistical Power

### Adding More Samples

The current analysis can be substantially strengthened by:

- **Increasing sample sizes per species**: More biological replicates within each age/sex category will improve statistical power and effect size estimation
- **Adding new species**: Comparative analyses across broader taxonomic ranges can reveal general patterns vs. lineage-specific effects
- **Balancing experimental design**: Ensuring equal representation of age classes, sexes, and species improves the robustness of statistical comparisons

Simply add new samples to your `SampleMetaData.csv` file with appropriate `Species`, `Sample_ID`, `Sex`, and `Status` information, run the pipeline, and the R scripts will automatically incorporate them into the analysis.

## Leveraging Publicly Available Data

### Mining NCBI Sequence Archives

A major advantage of this pipeline is its compatibility with publicly available whole genome sequencing (WGS) data. Many datasets already exist that include age and sex metadata, even if they were originally collected for other purposes:

**Finding suitable datasets:**

1. **NCBI SRA (Sequence Read Archive)**: Search for whole genome sequencing projects with metadata filters:
   - Keywords: "whole genome sequencing", species names of interest
   - Filters: age, sex, life stage information in metadata
   - Look for BioProjects with multiple individuals

2. **Example search strategy**:
   ```
   ("whole genome sequencing"[Strategy] OR "WGS"[Strategy])
   AND "age"[Attribute]
   AND "sex"[Attribute]
   AND [your_species_name]
   ```

3. **Repurposing existing data**: Population genomics studies, conservation genomics projects, and phylogeographic studies often include:
   - Age class information (juvenile, adult, aged)
   - Sex information from genetic or phenotypic data
   - Multiple individuals with sufficient mtDNA coverage

### Integrating GWS Data

**Key advantages of GWS data with mity:**

- **Higher coverage**: WGS typically provides more uniform and deeper coverage of the mitochondrial genome compared to PCR amplicons
- **Reduced amplification bias**: Avoids preferential amplification of certain mtDNA haplotypes that can occur with PCR
- **Genome-wide context**: Allows integration with nuclear genome analyses if desired

**Pipeline compatibility:**

This pipeline is **already compatible with WGS data**. The `mity` variant caller was designed to work with both:
- PCR-amplified mitochondrial DNA (as in the current dataset)
- Whole genome sequencing data (where mtDNA reads are extracted from total genomic DNA)

Simply provide WGS FASTQ files in the same directory structure and the pipeline will process them identically.

## Hypothesis Testing Opportunities

### Expanded Research Questions

With larger and more diverse datasets, this analytical framework enables testing of:

1. **Age-related accumulation of heteroplasmy**
   - Do heteroplasmy levels increase with age across taxa?
   - Are there species-specific rates of mtDNA mutation accumulation?
   - Do long-lived vs. short-lived species differ in age-related patterns?

2. **Sex-specific differences**
   - Do males and females show different heteroplasmy levels?
   - Are there sex × age interactions in mtDNA variation?
   - How does maternal inheritance affect patterns in males vs. females?

3. **Life history correlates**
   - Do metabolic rate, body size, or longevity predict heteroplasmy levels?
   - Are there differences between endotherms and ectotherms?
   - Do reproductive strategies correlate with mtDNA variation?

4. **Functional consequences**
   - Which genes/regions accumulate more variation?
   - Are nonsynonymous variants more common in specific proteins?
   - Do variants in certain genes (e.g., OXPHOS complexes) show stronger age effects?

5. **Comparative evolutionary patterns**
   - How do mutation spectra differ across vertebrate clades?
   - Are transition/transversion ratios consistent across taxa?
   - Do mutation hotspots exist at conserved positions?

### Modifying the Analysis Scripts

The R scripts are designed to be easily modified:

- **Filter criteria**: Adjust `TIER_THRESHOLD`, `MAX_HETEROPLASMY`, and `MIN_VARIANTS_PER_SAMPLE` in the configuration section
- **Statistical tests**: Switch between parametric (t-test) and non-parametric (Wilcoxon) tests via `STAT_TEST` parameter
- **Additional covariates**: Add new metadata columns (e.g., body mass, geographic location, sequencing depth) and incorporate them into linear models
- **Gene-specific analyses**: Filter by specific mitochondrial genes using the `GENE` column from mity output

## Practical Steps to Get Started

### 1. Identify and Download Data

```bash
# Example: Download SRA data
prefetch SRR1234567
fastq-dump --split-files --gzip SRR1234567
```

### 2. Update Metadata

Add entries to `SampleMetaData.csv`:

```csv
Species,Sample_ID,Sex,Status
Mus_musculus,SRR1234567,M,young
Mus_musculus,SRR1234568,F,old
```

### 3. Run the Pipeline

```bash
# Add new species data directory
mkdir -p Mus_musculus/raw_reads/
mv SRR*.fastq.gz Mus_musculus/raw_reads/

# Run pipeline
./run_mity_pipeline.sh Mus_musculus
```

### 4. Analyze Results

The R scripts will automatically include new species in all analyses and visualizations.

## Extending to Non-Model Organisms

### Why This Matters

Most heteroplasmy research focuses on traditional model organisms (humans, mice, *Drosophila*, *C. elegans*). However:

- **Biodiversity insights**: Non-model organisms represent >99% of animal diversity
- **Life history variation**: Non-model systems span extreme ranges of longevity, metabolic rates, and reproductive strategies
- **Evolutionary context**: Comparative approaches require taxonomic breadth
- **Conservation applications**: Understanding mtDNA variation can inform population health assessments

### Advantages of This Pipeline for Non-Model Systems

1. **No reference bias**: Works with any species with a mitochondrial reference genome
2. **De novo capability**: `mity` can work with assembled mitochondrial genomes from the same dataset
3. **Low computational requirements**: Can run on standard workstations
4. **Flexible metadata**: Accommodates diverse biological information beyond age/sex

## Example Future Projects

### Project 1: Comparative Aging Across Vertebrates

- **Hypothesis**: Long-lived species show slower accumulation of mtDNA heteroplasmy
- **Data**: Collect WGS data from NCBI for species with known maximum lifespans
- **Analysis**: Correlate heteroplasmy levels with species maximum lifespan and age class

### Project 2: Sex-Specific Selection on mtDNA

- **Hypothesis**: Males show higher heteroplasmy due to relaxed selection (maternal inheritance)
- **Data**: Population genomics datasets with sex information
- **Analysis**: Compare heteroplasmy levels and variant effects between sexes

### Project 3: Thermal Physiology and mtDNA Stability

- **Hypothesis**: Ectotherms show different heteroplasmy patterns than endotherms
- **Data**: WGS from diverse vertebrates across thermal niches
- **Analysis**: Test for correlations between body temperature, metabolic rate, and mtDNA variation

### Project 4: Repurposing Conservation Genomics Data

- **Hypothesis**: Population bottlenecks affect mtDNA heteroplasmy levels
- **Data**: Conservation genomics datasets often include age and demographic information
- **Analysis**: Compare heteroplasmy in stable vs. declining populations

## Conclusion

The analytical framework developed here is intentionally flexible and expandable. While the current analysis may not show statistically significant patterns due to sample size limitations, the infrastructure is in place to:

1. Seamlessly integrate new samples and species
2. Leverage publicly available WGS datasets from NCBI and other repositories
3. Repurpose existing genomic data collected for population genetics, phylogenomics, or conservation
4. Test diverse hypotheses about age-related and sex-specific mtDNA variation
5. Extend heteroplasmy research into the vast diversity of non-model organisms

By combining this pipeline with the rich genomic resources now available for diverse taxa, researchers can address fundamental questions about mitochondrial evolution, aging biology, and the functional consequences of mtDNA variation across the tree of life.

---

**Ready to expand your analysis?** Start by exploring NCBI SRA for relevant datasets, or reach out to collaborators who may have existing WGS data with age and sex information. The pipeline is ready when you are.
