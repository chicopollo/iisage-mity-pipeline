#!/usr/bin/env Rscript
################################################################################
# Status-Based Heteroplasmy Analysis (Old vs Young)
################################################################################
#
# This script analyzes mitochondrial heteroplasmy levels comparing old vs young
# individuals across species using the Status column from metadata.
#
# USAGE:
#   Rscript status_heteroplasmy_analysis.R
#   or run interactively in RStudio
#
# REQUIREMENTS:
#   - mity pipeline CSV outputs in species/mity_output/ directories
#   - SampleMetaData.csv with Species, Sample_ID, Sex, Status columns
#
# OUTPUTS:
#   - Interactive plots displayed in console/RStudio
#   - Summary statistics printed to console
#   - Statistical test results for old vs young comparisons
#
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Filter settings
TIER_THRESHOLD        <- 3      # Only include variants with TIER < this value
MAX_HETEROPLASMY      <- 1      # Exclude homoplasmic variants (heteroplasmy = 1)
MIN_VARIANTS_PER_SAMPLE <- 1    # Minimum variants per sample to include

# Which plots to generate?
PLOT_STATUS_OVERVIEW   <- TRUE   # Overview of old vs young across species
PLOT_STATUS_BOXPLOT    <- TRUE   # Boxplots comparing old vs young
PLOT_STATUS_BY_SPECIES <- TRUE   # Status comparison within each species
PLOT_STATUS_BY_SEX     <- TRUE   # Status x Sex interaction
PLOT_VARIANT_COUNTS    <- TRUE   # Number of variants by status
PLOT_EFFECT_SIZES      <- TRUE   # Effect sizes across species

# Species to exclude from analysis (if any)
EXCLUDE_SPECIES       <- c()    # e.g., c("Species1", "Species2")

# Pause between plots for interactive viewing?
PAUSE_BETWEEN_PLOTS   <- TRUE

# Statistical test to use: "t.test" or "wilcox" (Wilcoxon/Mann-Whitney)
STAT_TEST             <- "wilcox"  # "wilcox" is more robust for non-normal data

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

cat("Loading required libraries...\n")

suppressPackageStartupMessages({
  library(readr)      # Read CSV files
  library(dplyr)      # Data manipulation
  library(tidyr)      # Data tidying
  library(purrr)      # Functional programming
  library(ggplot2)    # Plotting
  library(broom)      # Statistical modeling
  library(scales)     # Scale functions
  library(stringr)    # String manipulation
})

cat("Libraries loaded.\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Pause between plots if in interactive mode
pause_plot <- function() {
  if (PAUSE_BETWEEN_PLOTS && interactive()) {
    cat("\nPress [Enter] to continue...")
    readline()
  }
}

# Calculate Cohen's d effect size
cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * sd(x)^2 + (ny - 1) * sd(y)^2) / (nx + ny - 2))
  d <- (mean(x) - mean(y)) / pooled_sd
  return(d)
}

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

cat(strrep("=", 78), "\n", sep = "")
cat("LOADING DATA\n")
cat(strrep("=", 78), "\n\n", sep = "")

# Find all mity CSV output files
cat("Searching for mity CSV output files...\n")
csv_files <- list.files(path = ".", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
csv_files <- grep("/mity_output/", csv_files, value = TRUE)

if (length(csv_files) == 0) {
  stop("ERROR: No CSV files found in mity_output directories.\n",
       "       Please check your directory structure.")
}

cat(sprintf("  Found %d CSV files\n", length(csv_files)))

# Extract species names from directory structure
species_names <- basename(dirname(dirname(csv_files)))

# Read all CSV files and combine
cat("Reading CSV files...\n")
all_mity <- map2_dfr(csv_files, species_names, function(file, sp) {
  read_csv(file, show_col_types = FALSE) %>%
    mutate(Species = sp)
})

cat(sprintf("  Loaded %d total variants\n", nrow(all_mity)))

# Load sample metadata
cat("Loading sample metadata...\n")
if (!file.exists("SampleMetaData.csv")) {
  stop("ERROR: SampleMetaData.csv not found in current directory.")
}

metadata <- read_csv("SampleMetaData.csv", show_col_types = FALSE)
cat(sprintf("  Loaded metadata for %d samples\n", nrow(metadata)))

# ==============================================================================
# DATA PROCESSING
# ==============================================================================

cat("\nProcessing and filtering data...\n")

# Extract Sample_ID from SAMPLE column and join with metadata
all_data <- all_mity %>%
  mutate(
    Sample_ID = sub("^([^_]+_[^_]+)_([^_]+)_.*$", "\\2", SAMPLE),
    # Convert to numeric
    `VARIANT HETEROPLASMY` = as.numeric(`VARIANT HETEROPLASMY`),
    `VARIANT QUALITY` = as.numeric(`VARIANT QUALITY`),
    `ALT DEPTH` = as.numeric(`ALT DEPTH`),
    TIER = as.numeric(TIER)
  ) %>%
  inner_join(metadata, by = c("Sample_ID", "Species"))

# Apply filters
filtered_data <- all_data %>%
  filter(
    TIER < TIER_THRESHOLD,
    `VARIANT HETEROPLASMY` < MAX_HETEROPLASMY,
    !is.na(`VARIANT HETEROPLASMY`),
    !is.na(Status),
    Status %in% c("old", "young")  # Only keep old and young
  )

# Exclude specific species if requested
if (length(EXCLUDE_SPECIES) > 0) {
  filtered_data <- filtered_data %>%
    filter(!Species %in% EXCLUDE_SPECIES)
}

cat(sprintf("  After filtering: %d variants\n", nrow(filtered_data)))
cat(sprintf("  Across %d species\n", n_distinct(filtered_data$Species)))
cat(sprintf("  Across %d samples\n", n_distinct(filtered_data$Sample_ID)))
cat(sprintf("  Old samples: %d, Young samples: %d\n\n",
            n_distinct(filtered_data$Sample_ID[filtered_data$Status == "old"]),
            n_distinct(filtered_data$Sample_ID[filtered_data$Status == "young"])))

# ==============================================================================
# CALCULATE SUMMARY STATISTICS
# ==============================================================================

cat(strrep("=", 78), "\n", sep = "")
cat("SUMMARY STATISTICS BY STATUS\n")
cat(strrep("=", 78), "\n\n", sep = "")

# Per-sample summary
sample_summary <- filtered_data %>%
  group_by(Species, Sample_ID, Sex, Status) %>%
  summarise(
    n_variants = n(),
    mean_heteroplasmy = mean(`VARIANT HETEROPLASMY`, na.rm = TRUE),
    median_heteroplasmy = median(`VARIANT HETEROPLASMY`, na.rm = TRUE),
    sd_heteroplasmy = sd(`VARIANT HETEROPLASMY`, na.rm = TRUE),
    min_heteroplasmy = min(`VARIANT HETEROPLASMY`, na.rm = TRUE),
    max_heteroplasmy = max(`VARIANT HETEROPLASMY`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_variants >= MIN_VARIANTS_PER_SAMPLE)

# Overall status summary
status_summary <- sample_summary %>%
  group_by(Status) %>%
  summarise(
    n_samples = n(),
    n_variants_total = sum(n_variants),
    mean_variants_per_sample = mean(n_variants),
    mean_heteroplasmy = mean(mean_heteroplasmy, na.rm = TRUE),
    sd_heteroplasmy = sd(mean_heteroplasmy, na.rm = TRUE),
    median_heteroplasmy = median(mean_heteroplasmy, na.rm = TRUE),
    .groups = "drop"
  )

# Print status summary
cat("Overall status comparison:\n\n")
print(status_summary, n = Inf)

# Species x Status summary
species_status_summary <- sample_summary %>%
  group_by(Species, Status) %>%
  summarise(
    n_samples = n(),
    mean_heteroplasmy = mean(mean_heteroplasmy, na.rm = TRUE),
    sd_heteroplasmy = sd(mean_heteroplasmy, na.rm = TRUE),
    median_heteroplasmy = median(mean_heteroplasmy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Status,
    values_from = c(n_samples, mean_heteroplasmy, sd_heteroplasmy, median_heteroplasmy),
    names_glue = "{Status}_{.value}"
  )

cat("\n\nSpecies x Status summary:\n\n")
print(species_status_summary, n = Inf)

# ==============================================================================
# VISUALIZATION
# ==============================================================================

cat("\n")
cat(strrep("=", 78), "\n", sep = "")
cat("GENERATING PLOTS\n")
cat(strrep("=", 78), "\n\n", sep = "")

# ------------------------------------------------------------------------------
# Plot 1: Overall status comparison across all species
# ------------------------------------------------------------------------------

if (PLOT_STATUS_OVERVIEW) {
  cat("Plot 1: Overall heteroplasmy comparison - Old vs Young\n")

  p <- ggplot(sample_summary,
              aes(x = Status, y = mean_heteroplasmy, fill = Status)) +
    geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
                 fill = "red", color = "black") +
    scale_fill_manual(values = c("young" = "#3498db", "old" = "#e74c3c")) +
    labs(
      title = "Mitochondrial Heteroplasmy: Old vs Young",
      subtitle = "All species combined (red diamonds = means, lines = quartiles)",
      x = "Status",
      y = "Mean Heteroplasmy per Sample"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "none"
    )

  print(p)
  pause_plot()
}

# ------------------------------------------------------------------------------
# Plot 2: Status boxplot comparison
# ------------------------------------------------------------------------------

if (PLOT_STATUS_BOXPLOT) {
  cat("Plot 2: Boxplot comparison of old vs young\n")

  p <- ggplot(sample_summary,
              aes(x = Status, y = mean_heteroplasmy, fill = Status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
                 fill = "yellow", color = "black") +
    scale_fill_manual(values = c("young" = "#3498db", "old" = "#e74c3c")) +
    labs(
      title = "Status Comparison of Heteroplasmy",
      subtitle = "Yellow diamonds indicate means; boxes show IQR",
      x = "Status",
      y = "Mean Heteroplasmy per Sample"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none"
    )

  print(p)
  pause_plot()
}

# ------------------------------------------------------------------------------
# Plot 3: Status comparison within each species
# ------------------------------------------------------------------------------

if (PLOT_STATUS_BY_SPECIES) {
  cat("Plot 3: Old vs Young comparison within each species\n")

  p <- ggplot(sample_summary,
              aes(x = Status, y = mean_heteroplasmy, fill = Status)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~ Species, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("young" = "#3498db", "old" = "#e74c3c")) +
    labs(
      title = "Old vs Young Heteroplasmy Comparison by Species",
      subtitle = "Each panel shows one species",
      x = "Status",
      y = "Mean Heteroplasmy per Sample",
      fill = "Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  print(p)
  pause_plot()
}


# ------------------------------------------------------------------------------
# Plot 4: Status x Sex interaction
# ------------------------------------------------------------------------------

if (PLOT_STATUS_BY_SEX) {
  cat("Plot 4: Status comparison by sex\n")

  # Filter to samples with known sex
  sex_data <- sample_summary %>%
    filter(!is.na(Sex), Sex %in% c("M", "F"))

  if (nrow(sex_data) > 0) {
    p <- ggplot(sex_data,
                aes(x = Status, y = mean_heteroplasmy, fill = Status)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
      facet_grid(Species ~ Sex, scales = "free_y",
                 labeller = labeller(Sex = c("M" = "Male", "F" = "Female"))) +
      scale_fill_manual(values = c("young" = "#3498db", "old" = "#e74c3c")) +
      labs(
        title = "Heteroplasmy: Status x Sex Interaction",
        subtitle = "Comparing old vs young within each sex and species",
        x = "Status",
        y = "Mean Heteroplasmy per Sample",
        fill = "Status"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom"
      )

    print(p)
    pause_plot()
  } else {
    cat("  Skipping: No sex data available\n")
  }
}

# ------------------------------------------------------------------------------
# Plot 5: Variant counts by status
# ------------------------------------------------------------------------------

if (PLOT_VARIANT_COUNTS) {
  cat("Plot 5: Number of variants by status\n")

  variant_counts <- sample_summary %>%
    group_by(Status) %>%
    summarise(
      total_variants = sum(n_variants),
      mean_per_sample = mean(n_variants),
      .groups = "drop"
    )

  p <- ggplot(variant_counts,
              aes(x = Status, y = total_variants, fill = Status)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = total_variants), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c("young" = "#3498db", "old" = "#e74c3c")) +
    labs(
      title = "Total Heteroplasmic Variants by Status",
      subtitle = "Number of variants passing quality filters",
      x = "Status",
      y = "Total Number of Variants"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14)
    )

  print(p)
  pause_plot()
}

# ------------------------------------------------------------------------------
# Plot 6: Effect sizes across species
# ------------------------------------------------------------------------------

if (PLOT_EFFECT_SIZES) {
  cat("Plot 6: Effect sizes (old - young) across species\n")

  # Calculate effect sizes for each species
  effect_sizes <- sample_summary %>%
    group_by(Species) %>%
    filter(all(c("old", "young") %in% Status)) %>%  # Only species with both groups
    summarise(
      n_old = sum(Status == "old"),
      n_young = sum(Status == "young"),
      mean_old = mean(mean_heteroplasmy[Status == "old"]),
      mean_young = mean(mean_heteroplasmy[Status == "young"]),
      diff = mean_old - mean_young,
      cohens_d = cohens_d(
        mean_heteroplasmy[Status == "old"],
        mean_heteroplasmy[Status == "young"]
      ),
      .groups = "drop"
    ) %>%
    mutate(
      effect_direction = ifelse(diff > 0, "Old > Young", "Young > Old"),
      abs_d = abs(cohens_d)
    )

  p <- ggplot(effect_sizes,
              aes(x = reorder(Species, cohens_d), y = cohens_d,
                  fill = effect_direction)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray50") +
    coord_flip() +
    scale_fill_manual(values = c("Old > Young" = "#e74c3c",
                                  "Young > Old" = "#3498db")) +
    labs(
      title = "Effect Sizes: Old vs Young Heteroplasmy",
      subtitle = "Cohen's d (positive = old higher, negative = young higher)",
      x = "Species",
      y = "Cohen's d",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  print(p)
  pause_plot()
}

# ==============================================================================
# STATISTICAL COMPARISONS
# ==============================================================================

cat("\n")
cat(strrep("=", 78), "\n", sep = "")
cat("STATISTICAL COMPARISONS\n")
cat(strrep("=", 78), "\n\n", sep = "")

# Overall test: Old vs Young across all species
cat("Overall comparison of heteroplasmy: Old vs Young (all species combined)\n")
cat(sprintf("Statistical test: %s\n\n", STAT_TEST))

old_data <- sample_summary %>% filter(Status == "old") %>% pull(mean_heteroplasmy)
young_data <- sample_summary %>% filter(Status == "young") %>% pull(mean_heteroplasmy)

if (STAT_TEST == "wilcox") {
  overall_test <- wilcox.test(old_data, young_data)
  cat("Wilcoxon rank-sum test (Mann-Whitney U test):\n")
} else {
  overall_test <- t.test(old_data, young_data)
  cat("Welch's t-test:\n")
}

print(overall_test)

# Calculate effect size
overall_d <- cohens_d(old_data, young_data)
cat(sprintf("\nCohen's d effect size: %.3f\n", overall_d))
cat(sprintf("Interpretation: %s\n",
            ifelse(abs(overall_d) < 0.2, "negligible",
            ifelse(abs(overall_d) < 0.5, "small",
            ifelse(abs(overall_d) < 0.8, "medium", "large")))))

# Species-specific tests
cat("\n\nSpecies-specific comparisons (Old vs Young):\n")
cat(strrep("-", 78), "\n", sep = "")

# Get list of species with both old and young samples
species_with_both <- sample_summary %>%
  group_by(Species) %>%
  filter(all(c("old", "young") %in% Status)) %>%
  pull(Species) %>%
  unique()

# Function to run test for a single species
run_species_test <- function(species_name) {
  old_vals <- sample_summary %>%
    filter(Species == species_name, Status == "old") %>%
    pull(mean_heteroplasmy)

  young_vals <- sample_summary %>%
    filter(Species == species_name, Status == "young") %>%
    pull(mean_heteroplasmy)

  # Run statistical test
  if (STAT_TEST == "wilcox") {
    test <- wilcox.test(old_vals, young_vals)
  } else {
    test <- t.test(old_vals, young_vals)
  }

  # Return results as a tibble
  tibble(
    Species = species_name,
    n_old = length(old_vals),
    n_young = length(young_vals),
    mean_old = mean(old_vals),
    mean_young = mean(young_vals),
    p_value = test$p.value,
    statistic = as.numeric(test$statistic),
    cohens_d = cohens_d(old_vals, young_vals)
  )
}

# Run tests for all species
species_tests <- map_dfr(species_with_both, run_species_test) %>%
  arrange(p_value)

print(species_tests, n = Inf)

# Apply multiple testing correction
species_tests <- species_tests %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "fdr"),
    significant = p_adjusted < 0.05
  )

cat("\n\nWith FDR correction (Benjamini-Hochberg):\n")
print(species_tests %>% select(Species, p_value, p_adjusted, significant, cohens_d),
      n = Inf)

sig_species <- species_tests %>% filter(significant)
if (nrow(sig_species) > 0) {
  cat("\n\nSpecies with significant differences (FDR-adjusted p < 0.05):\n")
  print(sig_species %>% select(Species, n_old, n_young, mean_old, mean_young,
                                p_adjusted, cohens_d), n = Inf)
} else {
  cat("\n\nNo species show significant differences after FDR correction.\n")
}

# Test for Sex x Status interaction (if sex data available)
sex_data <- sample_summary %>%
  filter(!is.na(Sex), Sex %in% c("M", "F"))

if (nrow(sex_data) > 10) {
  cat("\n\nTesting Sex x Status interaction:\n")
  cat(strrep("-", 78), "\n", sep = "")

  # Two-way ANOVA
  aov_result <- aov(mean_heteroplasmy ~ Status * Sex + Species, data = sex_data)
  cat("\nTwo-way ANOVA (Status x Sex, controlling for Species):\n")
  print(summary(aov_result))
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat(strrep("=", 78), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 78), "\n\n", sep = "")

cat("Key findings:\n")
cat(sprintf("  - Analyzed %d species\n", n_distinct(sample_summary$Species)))
cat(sprintf("  - Old samples: %d (mean heteroplasmy: %.4f)\n",
            nrow(sample_summary %>% filter(Status == "old")),
            mean(old_data)))
cat(sprintf("  - Young samples: %d (mean heteroplasmy: %.4f)\n",
            nrow(sample_summary %>% filter(Status == "young")),
            mean(young_data)))
cat(sprintf("  - Overall effect size (Cohen's d): %.3f\n", overall_d))
cat(sprintf("  - Overall p-value: %.4f\n", overall_test$p.value))

if (nrow(sig_species) > 0) {
  cat(sprintf("  - Species with significant old/young differences: %d\n",
              nrow(sig_species)))
  cat(sprintf("    %s\n",
              paste(sig_species$Species, collapse = ", ")))
} else {
  cat("  - No individual species showed significant differences after correction\n")
}

cat("\nTo customize this analysis, edit the CONFIGURATION section at the top.\n")
cat("You can adjust filters, select specific plots, and change statistical tests.\n\n")
