#!/usr/bin/env Rscript
################################################################################
# MITY Pipeline Statistical Analysis and Figure Generation
################################################################################
#
# This script performs comprehensive mitochondrial heteroplasmy analysis using
# outputs from the mity pipeline. It generates two complementary types of
# analyses:
#   1. Sample-level heteroplasmy metrics (from CSV outputs)
#   2. Gene-level variant analysis (from VCF files and GFF3 annotations)
#
# USAGE:
#   Rscript mity_analysis.R
#   or run interactively in RStudio
#
# REQUIREMENTS:
#   - mity pipeline outputs (CSV files and VCF files)
#   - GFF3 annotation files for each species
#   - SampleMetaData.csv with sample information
#
# CUSTOMIZATION:
#   See the CONFIGURATION section below to:
#   - Select which species to analyze
#   - Choose which plots to generate
#   - Adjust filtering thresholds
#
################################################################################

# ==============================================================================
# CONFIGURATION - Modify these settings as needed
# ==============================================================================

# Which analyses to run?
RUN_HETEROPLASMY_ANALYSIS <- TRUE   # Sample-level heteroplasmy from CSV
RUN_GENE_ANALYSIS         <- TRUE   # Gene-level analysis from VCF

# Which plots to generate for heteroplasmy analysis?
PLOT_QUALITY_DEPTH        <- TRUE   # Alt depth vs variant quality
PLOT_QUALITY_HETERO       <- TRUE   # Variant quality vs heteroplasmy
PLOT_DEPTH_HETERO         <- TRUE   # Alt depth vs heteroplasmy
PLOT_POSITION             <- TRUE   # Cohort frequency by position
PLOT_SEX                  <- TRUE   # Heteroplasmy metrics by sex
PLOT_AGE_CHRON            <- TRUE   # Heteroplasmy metrics by chronological age
PLOT_AGE_PERCENT          <- TRUE   # Heteroplasmy metrics by age percent
PLOT_POS_HETERO           <- TRUE   # Variant heteroplasmy by position

# Which plots to generate for gene analysis?
PLOT_GENE_BOXPLOTS        <- TRUE   # Mean heteroplasmy by gene
PLOT_FUNCTIONAL_REGIONS   <- TRUE   # Heteroplasmy by functional region
PLOT_GENE_HEATMAP         <- TRUE   # Gene x Sample heatmap
PLOT_VARIABLE_GENES       <- TRUE   # Most variable genes

# Filter settings
TIER_THRESHOLD            <- 3      # Only include variants with TIER < this
MAX_HETEROPLASMY          <- 1      # Exclude variants with heteroplasmy >= this
MIN_VAF                   <- 0      # Minimum variant allele frequency
MAX_VAF                   <- 1      # Maximum variant allele frequency

# Specific species to analyze (leave empty to analyze all)
SPECIES_FILTER            <- c()    # e.g., c("Musca_domestica", "Desmodus_rotundus")

# Output settings
PAUSE_BETWEEN_PLOTS       <- TRUE   # Wait for user input between plots?

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

cat("Loading required libraries...\n")

# Data manipulation and visualization
suppressPackageStartupMessages({
  library(readr)      # Read CSV/TSV files
  library(dplyr)      # Data manipulation
  library(purrr)      # Functional programming tools
  library(tidyr)      # Data tidying
  library(ggplot2)    # Plotting
  library(stringr)    # String manipulation
  library(broom)      # Tidy statistical outputs
  library(scales)     # Scale functions for plots (pvalue formatting)

  # Only load vcfR if we're running gene analysis
  if (RUN_GENE_ANALYSIS) {
    library(vcfR)     # Read and process VCF files
  }
})

cat("Libraries loaded successfully.\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# pause_plot: Pause between plots for interactive viewing
# ------------------------------------------------------------------------------
pause_plot <- function() {
  if (PAUSE_BETWEEN_PLOTS && interactive()) {
    cat("\nPress [Enter] to continue to next plot...")
    readline()
  }
}

# ------------------------------------------------------------------------------
# mk_age_labeller: Create facet labels with statistical test results
#
# For 2 age groups: performs a t-test
# For 3+ age groups: performs linear regression
#
# Parameters:
#   df_sum: Long-format data frame with Metric, Value, and age variable
#   age_var: Name of age variable (e.g., "Age_chron" or "Age_percent")
#
# Returns:
#   Named vector for use with as_labeller() in facet_wrap()
# ------------------------------------------------------------------------------
mk_age_labeller <- function(df_sum, age_var) {
  av <- rlang::sym(age_var)

  # Count distinct ages per metric
  n_by_metric <- df_sum %>%
    group_by(Metric) %>%
    summarise(n_ages = n_distinct(!!av, na.rm = TRUE), .groups = "drop")

  # Build labels with appropriate statistical test
  out <- n_by_metric %>%
    rowwise() %>%
    do({
      met <- .$Metric
      nA <- .$n_ages
      dat <- df_sum %>% filter(Metric == met, !is.na(!!av), !is.na(Value))

      if (nrow(dat) == 0) {
        # No data available
        tibble(Metric = met, strip = sprintf("%s\nno data", met))
      } else if (nA == 2) {
        # Two groups: use t-test
        tt <- tidy(t.test(Value ~ factor(!!av), data = dat))
        ages <- sort(unique(pull(dat, !!av)))
        tibble(Metric = met,
               strip = sprintf("%s\nt = %.2f, p = %s (%g vs %g)",
                               met, tt$statistic, pvalue(tt$p.value, accuracy = 1e-3),
                               ages[1], ages[2]))
      } else {
        # Three or more groups: use linear regression
        lm_t <- tidy(lm(Value ~ !!av, data = dat)) %>%
          filter(term == rlang::as_label(av))
        tibble(Metric = met,
               strip = sprintf("%s\nslope t = %.2f, p = %s",
                               met, lm_t$statistic, pvalue(lm_t$p.value, accuracy = 1e-3)))
      }
    }) %>%
    ungroup()

  setNames(out$strip, out$Metric)
}

# ------------------------------------------------------------------------------
# read_gene_annotations: Parse GFF3 file to extract gene annotations
#
# Parameters:
#   gff_path: Path to GFF3 file (can be gzipped)
#
# Returns:
#   Data frame with start, end, gene_type, and name columns
# ------------------------------------------------------------------------------
read_gene_annotations <- function(gff_path) {
  gff <- read_tsv(gff_path, comment = "#", col_names = FALSE,
                  show_col_types = FALSE)
  colnames(gff) <- c("seqname", "source", "type", "start", "end",
                     "score", "strand", "phase", "attributes")

  gff %>%
    filter(type %in% c("CDS", "tRNA", "rRNA")) %>%
    mutate(
      gene = str_extract(attributes, "gene=([^;]+)") %>% str_remove("gene="),
      product = str_extract(attributes, "product=([^;]+)") %>% str_remove("product="),
      name = coalesce(gene, product, type),
      gene_type = type
    ) %>%
    select(start, end, gene_type, name)
}

# ------------------------------------------------------------------------------
# assign_to_genes: Assign variants to genes based on position
#
# Parameters:
#   df_vaf: Data frame with variant positions
#   gff_data: Gene annotation data from read_gene_annotations()
#
# Returns:
#   df_vaf with added gene and gene_type columns
# ------------------------------------------------------------------------------
assign_to_genes <- function(df_vaf, gff_data) {
  df_vaf %>%
    rowwise() %>%
    mutate(
      gene = {
        hits <- gff_data %>% filter(pos >= start & pos <= end)
        if (nrow(hits) > 0) hits$name[1] else "Intergenic"
      },
      gene_type = {
        hits <- gff_data %>% filter(pos >= start & pos <= end)
        if (nrow(hits) > 0) hits$gene_type[1] else "Intergenic"
      }
    ) %>%
    ungroup()
}

# ==============================================================================
# HETEROPLASMY ANALYSIS (CSV-based)
# ==============================================================================

if (RUN_HETEROPLASMY_ANALYSIS) {
  cat(strrep("=", 78), "\n", sep = "")
  cat("PART 1: SAMPLE-LEVEL HETEROPLASMY ANALYSIS\n")
  cat(strrep("=", 78), "\n\n", sep = "")

  # ----------------------------------------------------------------------------
  # Load data
  # ----------------------------------------------------------------------------

  cat("Loading mity CSV output files...\n")

  # Find all CSV files in mity_output directories
  csv_files <- list.files(path = ".", pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
  csv_files <- grep("/mity_output/", csv_files, value = TRUE)

  if (length(csv_files) == 0) {
    cat("WARNING: No CSV files found in mity_output directories.\n")
    cat("         Skipping heteroplasmy analysis.\n\n")
  } else {
    # Extract species names from file paths
    species <- basename(dirname(dirname(csv_files)))

    # Read all CSV files and combine with species labels
    all_mity <- map2_dfr(csv_files, species, ~ {
      read_csv(.x, show_col_types = FALSE) %>%
        mutate(Species = .y)
    })

    cat(sprintf("  Loaded %d variants from %d files\n",
                nrow(all_mity), length(csv_files)))

    # Load sample metadata
    cat("Loading sample metadata...\n")
    metadata <- read_csv("SampleMetaData.csv", show_col_types = FALSE)

    # ----------------------------------------------------------------------------
    # Join data and prepare numeric columns
    # ----------------------------------------------------------------------------

    cat("Joining variant data with metadata...\n")

    # Extract Sample_ID from the SAMPLE column
    # Pattern: Species_SampleID_... -> extract SampleID
    all_mity_meta <- all_mity %>%
      mutate(Sample_ID = sub("^([^_]+_[^_]+)_([^_]+)_.*$", "\\2", SAMPLE)) %>%
      inner_join(metadata, by = c("Sample_ID", "Species")) %>%
      mutate(
        # Convert string columns to numeric
        `VARIANT HETEROPLASMY` = as.numeric(`VARIANT HETEROPLASMY`),
        `VARIANT QUALITY`      = as.numeric(`VARIANT QUALITY`),
        `ALT DEPTH`            = as.numeric(`ALT DEPTH`),
        Age_chron              = suppressWarnings(as.numeric(`Age_(chron)`)),
        Age_percent            = suppressWarnings(as.numeric(`Age_(%)`))
      )

    cat(sprintf("  Final dataset: %d variants across %d samples\n\n",
                nrow(all_mity_meta), n_distinct(all_mity_meta$Sample_ID)))

    # ----------------------------------------------------------------------------
    # Analyze each species
    # ----------------------------------------------------------------------------

    # Split data by species
    sp_list <- split(all_mity_meta, all_mity_meta$Species)
    sp_list <- sp_list[vapply(sp_list, nrow, integer(1)) > 0]

    # Filter to specific species if requested
    if (length(SPECIES_FILTER) > 0) {
      sp_list <- sp_list[names(sp_list) %in% SPECIES_FILTER]
    }

    cat(sprintf("Analyzing %d species...\n\n", length(sp_list)))

    # Process each species
    for (sp_name in names(sp_list)) {
      cat("\n")
      cat(strrep("-", 78), "\n", sep = "")
      cat(sprintf("Species: %s\n", sp_name))
      cat(strrep("-", 78), "\n\n", sep = "")

      df <- sp_list[[sp_name]]

      # ------------------------------------------------------------------------
      # Prepare summary statistics
      # ------------------------------------------------------------------------

      # Calculate per-sample heteroplasmy metrics
      df_sum <- df %>%
        filter(TIER < TIER_THRESHOLD, `VARIANT HETEROPLASMY` < MAX_HETEROPLASMY) %>%
        group_by(Sample_ID, Sex, Age_chron, Age_percent) %>%
        summarise(
          meanVariant = mean(`VARIANT HETEROPLASMY`, na.rm = TRUE),
          minVariant  = suppressWarnings(min(`VARIANT HETEROPLASMY`, na.rm = TRUE)),
          varVariant  = var(`VARIANT HETEROPLASMY`, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        pivot_longer(
          c(meanVariant, minVariant, varVariant),
          names_to = "Metric",
          values_to = "Value"
        ) %>%
        mutate(Metric = recode(Metric,
                               meanVariant = "Mean heteroplasmy",
                               minVariant  = "Min heteroplasmy",
                               varVariant  = "Variance of heteroplasmy"))

      # Create facet labels with statistical tests
      sex_lab <- df_sum %>%
        group_by(Metric) %>%
        do(tidy(t.test(Value ~ Sex, data = .))) %>%
        transmute(Metric,
                  strip = sprintf("%s\nt = %.2f, p = %s",
                                  Metric, statistic, pvalue(p.value, accuracy = 1e-3))) %>%
        { setNames(.$strip, .$Metric) }

      agec_lab <- mk_age_labeller(df_sum, "Age_chron")
      agep_lab <- mk_age_labeller(df_sum, "Age_percent")

      # ------------------------------------------------------------------------
      # Generate plots
      # ------------------------------------------------------------------------

      # Plot 1: Alt depth vs Variant quality (colored by TIER)
      if (PLOT_QUALITY_DEPTH) {
        cat("  Generating: Alt depth vs Variant quality plot...\n")
        p <- ggplot(df, aes(log(`VARIANT QUALITY`), log(`ALT DEPTH`), colour = TIER)) +
          geom_point(alpha = 0.7) +
          labs(title = paste(sp_name, "- Alt depth vs Variant quality by Tier"),
               x = "log(Variant quality)",
               y = "log(Alt depth)") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 2: Variant quality vs Heteroplasmy (colored by TIER)
      if (PLOT_QUALITY_HETERO) {
        cat("  Generating: Variant quality vs Heteroplasmy plot...\n")
        p <- ggplot(df, aes(log(`VARIANT QUALITY`), `VARIANT HETEROPLASMY`, colour = TIER)) +
          geom_point(alpha = 0.7) +
          labs(title = paste(sp_name, "- Variant quality vs Heteroplasmy"),
               x = "log(Variant quality)",
               y = "Variant heteroplasmy") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 3: Alt depth vs Heteroplasmy (colored by TIER)
      if (PLOT_DEPTH_HETERO) {
        cat("  Generating: Alt depth vs Heteroplasmy plot...\n")
        p <- ggplot(df, aes(log10(`ALT DEPTH`), `VARIANT HETEROPLASMY`, colour = TIER)) +
          geom_point(alpha = 0.7) +
          labs(title = paste(sp_name, "- Alt depth vs Heteroplasmy"),
               x = "log10(Alt depth)",
               y = "Variant heteroplasmy") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 4: Cohort frequency by position
      if (PLOT_POSITION) {
        cat("  Generating: Cohort frequency by position plot...\n")
        p <- df %>%
          filter(TIER < TIER_THRESHOLD) %>%
          ggplot(aes(POS, `COHORT FREQUENCY`, colour = `VARIANT HETEROPLASMY`)) +
          geom_point(size = 2) +
          labs(title = paste(sp_name, "- Cohort frequency by position"),
               x = "Position",
               y = "Cohort frequency") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 5: Heteroplasmy metrics by Sex (with t-test results)
      if (PLOT_SEX) {
        cat("  Generating: Heteroplasmy metrics by Sex plot...\n")
        p <- ggplot(df_sum, aes(x = Sex, y = Value, colour = Sex)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.85,
                      show.legend = FALSE) +
          facet_wrap(~ Metric, ncol = 1, scales = "free_y",
                     labeller = as_labeller(sex_lab)) +
          labs(title = paste(sp_name, "- Heteroplasmy metrics by Sex"),
               x = "Sex",
               y = "Value") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 6: Heteroplasmy metrics by Age (chronological)
      if (PLOT_AGE_CHRON && any(!is.na(df_sum$Age_chron))) {
        cat("  Generating: Heteroplasmy metrics by Age (chronological) plot...\n")
        p <- ggplot(df_sum, aes(x = factor(Age_chron), y = Value)) +
          geom_boxplot(outlier.shape = NA, colour = "grey30") +
          geom_jitter(aes(colour = factor(Age_chron)), width = 0.15, height = 0,
                      size = 1.8, alpha = 0.85, show.legend = FALSE) +
          facet_wrap(~ Metric, ncol = 1, scales = "free_y",
                     labeller = as_labeller(agec_lab)) +
          labs(title = paste(sp_name, "- Heteroplasmy metrics by Age (chron)"),
               x = "Age (chron)",
               y = "Value") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 7: Heteroplasmy metrics by Age (percent)
      if (PLOT_AGE_PERCENT && any(!is.na(df_sum$Age_percent))) {
        cat("  Generating: Heteroplasmy metrics by Age (%) plot...\n")
        p <- ggplot(df_sum, aes(x = factor(Age_percent), y = Value)) +
          geom_boxplot(outlier.shape = NA, colour = "grey30") +
          geom_jitter(aes(colour = factor(Age_percent)), width = 0.15, height = 0,
                      size = 1.8, alpha = 0.85, show.legend = FALSE) +
          facet_wrap(~ Metric, ncol = 1, scales = "free_y",
                     labeller = as_labeller(agep_lab)) +
          labs(title = paste(sp_name, "- Heteroplasmy metrics by Age (%)"),
               x = "Age (%)",
               y = "Value") +
          theme_minimal()
        print(p)
        pause_plot()
      }

      # Plot 8: Variant heteroplasmy by position (colored by Age)
      if (PLOT_POS_HETERO) {
        cat("  Generating: Variant heteroplasmy by position plot...\n")
        p <- df %>%
          filter(TIER <= TIER_THRESHOLD, `COHORT FREQUENCY` < 1) %>%
          ggplot(aes(factor(POS), `VARIANT HETEROPLASMY`)) +
          geom_boxplot(aes(colour = factor(Age_chron)), outlier.shape = NA) +
          labs(title = paste(sp_name, "- Variant heteroplasmy by position"),
               x = "Position",
               y = "Variant heteroplasmy",
               colour = "Age (chron)") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        pause_plot()
      }

      cat(sprintf("  Completed analysis for %s\n", sp_name))
    }

    cat("\n")
    cat("Heteroplasmy analysis complete.\n\n")
  }
}

# ==============================================================================
# GENE-LEVEL ANALYSIS (VCF-based)
# ==============================================================================

if (RUN_GENE_ANALYSIS) {
  cat(strrep("=", 78), "\n", sep = "")
  cat("PART 2: GENE-LEVEL VARIANT ANALYSIS\n")
  cat(strrep("=", 78), "\n\n", sep = "")

  # ----------------------------------------------------------------------------
  # Find VCF and GFF files
  # ----------------------------------------------------------------------------

  cat("Locating VCF and GFF3 files...\n")

  vcf_files <- list.files(".", pattern = "normalise\\.annotated\\.vcf\\.gz$",
                          recursive = TRUE, full.names = TRUE)
  vcf_files <- grep("/mity_chunked/", vcf_files, value = TRUE)

  gff_files <- list.files(".", pattern = "\\.gff3\\.gz$",
                          recursive = TRUE, full.names = TRUE)
  gff_files <- grep("/ref_mito", gff_files, value = TRUE)

  if (length(vcf_files) == 0 || length(gff_files) == 0) {
    cat("WARNING: VCF or GFF3 files not found.\n")
    cat("         VCF files found:", length(vcf_files), "\n")
    cat("         GFF3 files found:", length(gff_files), "\n")
    cat("         Skipping gene-level analysis.\n\n")
  } else {
    # Extract species names
    species_vcf <- sub("^\\./([^/]+)/.*$", "\\1", vcf_files)
    species_gff <- sub("^\\./([^/]+)/.*$", "\\1", gff_files)

    cat(sprintf("  Found %d VCF files and %d GFF3 files\n\n",
                length(vcf_files), length(gff_files)))

    # Load metadata
    cat("Loading sample metadata...\n")
    metadata <- read_csv("SampleMetaData.csv", show_col_types = FALSE)

    # ----------------------------------------------------------------------------
    # Process each species
    # ----------------------------------------------------------------------------

    # Get unique species with both VCF and GFF
    available_species <- intersect(species_vcf, species_gff)

    # Filter to specific species if requested
    if (length(SPECIES_FILTER) > 0) {
      available_species <- intersect(available_species, SPECIES_FILTER)
    }

    cat(sprintf("Analyzing %d species...\n\n", length(available_species)))

    for (sp_name in available_species) {
      cat("\n")
      cat(strrep("-", 78), "\n", sep = "")
      cat(sprintf("Species: %s\n", sp_name))
      cat(strrep("-", 78), "\n\n", sep = "")

      # Find files for this species
      vcf_idx <- which(species_vcf == sp_name)[1]
      gff_idx <- which(species_gff == sp_name)[1]

      vcf_path <- vcf_files[vcf_idx]
      gff_path <- gff_files[gff_idx]

      # ------------------------------------------------------------------------
      # Load and process data
      # ------------------------------------------------------------------------

      cat("  Reading VCF file...\n")
      vcf <- read.vcfR(vcf_path, verbose = FALSE)

      cat("  Reading GFF3 annotations...\n")
      gff_data <- read_gene_annotations(gff_path)

      cat("  Extracting allele depth information...\n")
      ad_mat <- extract.gt(vcf, element = "AD", as.numeric = FALSE)
      positions <- as.integer(getPOS(vcf))

      # Convert allele depth matrix to long format and calculate VAF
      df_vaf <- as_tibble(ad_mat, rownames = "variant") %>%
        mutate(pos = positions) %>%
        pivot_longer(cols = -c(variant, pos),
                     names_to = "sample",
                     values_to = "ad") %>%
        separate(ad, into = c("ref_reads", "alt_reads"),
                 sep = ",", convert = TRUE) %>%
        mutate(
          vaf = alt_reads / (ref_reads + alt_reads),
          Sample_ID = str_extract(sample, "[A-Za-z]+\\d+")
        ) %>%
        filter(vaf > MIN_VAF, vaf < MAX_VAF)

      # Join with metadata
      cat("  Joining with metadata...\n")
      df_vaf <- df_vaf %>%
        left_join(metadata %>% filter(Species == sp_name), by = "Sample_ID")

      # Check if we have required columns
      if (!("Sex" %in% colnames(df_vaf))) {
        cat("  WARNING: No Sex information in metadata for", sp_name, "\n")
        cat("           Skipping this species.\n")
        next
      }

      # Assign variants to genes
      cat("  Assigning variants to genes...\n")
      df_vaf <- assign_to_genes(df_vaf, gff_data)

      cat(sprintf("  Processed %d variants across %d samples\n",
                  nrow(df_vaf), n_distinct(df_vaf$Sample_ID)))
      cat(sprintf("  Identified %d genes\n\n", n_distinct(df_vaf$gene)))

      # ------------------------------------------------------------------------
      # Generate plots
      # ------------------------------------------------------------------------

      # Plot 1: Mean heteroplasmy by gene
      if (PLOT_GENE_BOXPLOTS) {
        cat("  Generating: Mean heteroplasmy by gene plot...\n")

        gene_summary <- df_vaf %>%
          group_by(gene, gene_type, Sex, Status) %>%
          summarise(mean_vaf = mean(vaf, na.rm = TRUE), .groups = "drop") %>%
          filter(gene != "Intergenic")

        if (nrow(gene_summary) > 0) {
          p <- ggplot(gene_summary, aes(x = reorder(gene, mean_vaf),
                                        y = mean_vaf, fill = gene_type)) +
            geom_boxplot() +
            facet_grid(Sex ~ Status, scales = "free_x") +
            coord_flip() +
            scale_fill_manual(values = c(
              "CDS" = "steelblue",
              "tRNA" = "orange",
              "rRNA" = "green3"
            )) +
            labs(
              x = "Gene",
              y = "Mean VAF",
              title = paste(sp_name, "- Mean heteroplasmy by gene"),
              fill = "Gene type"
            ) +
            theme_minimal() +
            theme(legend.position = "bottom")
          print(p)
          pause_plot()
        }
      }

      # Plot 2: Heteroplasmy by functional region
      if (PLOT_FUNCTIONAL_REGIONS) {
        cat("  Generating: Heteroplasmy by functional region plot...\n")

        func_summary <- df_vaf %>%
          group_by(Sample_ID, Sex, Status, gene_type) %>%
          summarise(mean_vaf = mean(vaf, na.rm = TRUE), .groups = "drop")

        if (nrow(func_summary) > 0) {
          p <- ggplot(func_summary, aes(x = gene_type, y = mean_vaf,
                                        fill = gene_type)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            facet_grid(Sex ~ Status) +
            scale_fill_manual(values = c(
              "CDS" = "steelblue",
              "tRNA" = "orange",
              "rRNA" = "green3",
              "Intergenic" = "gray70"
            )) +
            labs(
              x = "Functional Region",
              y = "Mean VAF per Sample",
              title = paste(sp_name, "- Heteroplasmy by functional region"),
              fill = "Region type"
            ) +
            theme_minimal() +
            theme(legend.position = "none")
          print(p)
          pause_plot()
        }
      }

      # Plot 3: Gene x Sample heatmap
      if (PLOT_GENE_HEATMAP) {
        cat("  Generating: Gene x Sample heatmap...\n")

        gene_sample <- df_vaf %>%
          filter(gene != "Intergenic") %>%
          group_by(gene, Sample_ID, Sex, Status) %>%
          summarise(mean_vaf = mean(vaf, na.rm = TRUE), .groups = "drop") %>%
          mutate(sample_group = paste(Sex, Status, Sample_ID, sep = "_"))

        if (nrow(gene_sample) > 0) {
          p <- ggplot(gene_sample, aes(x = sample_group, y = gene,
                                       fill = mean_vaf)) +
            geom_tile(color = "grey", linewidth = 0.5) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                 midpoint = 0.5, limits = c(0, 1)) +
            facet_grid(. ~ Sex + Status, scales = "free_x", space = "free_x") +
            labs(
              x = "Sample",
              y = "Gene",
              title = paste(sp_name, "- Heteroplasmy heatmap"),
              fill = "Mean VAF"
            ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              panel.grid = element_blank()
            )
          print(p)
          pause_plot()
        }
      }

      # Plot 4: Most variable genes
      if (PLOT_VARIABLE_GENES) {
        cat("  Generating: Most variable genes plot...\n")

        gene_var <- df_vaf %>%
          filter(gene != "Intergenic") %>%
          group_by(gene, gene_type) %>%
          summarise(
            mean_vaf = mean(vaf, na.rm = TRUE),
            sd_vaf = sd(vaf, na.rm = TRUE),
            cv = sd_vaf / mean_vaf,
            .groups = "drop"
          ) %>%
          filter(!is.na(cv), is.finite(cv)) %>%
          slice_max(order_by = cv, n = 15)

        if (nrow(gene_var) > 0) {
          p <- ggplot(gene_var, aes(x = reorder(gene, cv), y = cv,
                                    fill = gene_type)) +
            geom_col() +
            coord_flip() +
            scale_fill_manual(values = c(
              "CDS" = "steelblue",
              "tRNA" = "orange",
              "rRNA" = "green3"
            )) +
            labs(
              x = "Gene",
              y = "Coefficient of Variation",
              title = paste(sp_name, "- Most variable genes"),
              fill = "Gene type"
            ) +
            theme_minimal()
          print(p)
          pause_plot()
        }
      }

      cat(sprintf("  Completed analysis for %s\n", sp_name))
    }

    cat("\n")
    cat("Gene-level analysis complete.\n\n")
  }
}

# ==============================================================================
# ANALYSIS COMPLETE
# ==============================================================================

cat(strrep("=", 78), "\n", sep = "")
cat("ALL ANALYSES COMPLETE\n")
cat(strrep("=", 78), "\n\n", sep = "")

cat("To customize this analysis, edit the CONFIGURATION section at the top\n")
cat("of this script to:\n")
cat("  - Select which analyses to run\n")
cat("  - Choose which plots to generate\n")
cat("  - Filter to specific species\n")
cat("  - Adjust filtering thresholds\n\n")
