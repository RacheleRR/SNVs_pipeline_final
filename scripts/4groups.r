# ============================================================
# SNV ANALYSIS - EXECUTION SCRIPT
# Uses functions from all_together.r
# Covers: 4-group (SCZ/BD/Converter/Non_Converter)
# Both private and non-private variants
# ============================================================

# ============================================================
# 1. LIBRARIES
# ============================================================

library(readr)
library(data.table)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)
library(purrr)
library(broom)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(openxlsx)
library(gridExtra)
library(ggsignif)
library(FSA)
library(dunn.test)
library(RColorBrewer)

# Source all refactored functions
setwd("/home/rachele/SNVs/try/new/SNVs_pipeline_final/scripts")
source("Prep_data.r")
source("report.r")
source("statistics functions.r")
source("PLOTS.r")


# ============================================================
# 2. PARAMETERS — edit these to change paths or outputs
# ============================================================

# Input paths
VARIANT_DIR        <- "/media/rachele/DATA/SNVs/05_tsv"
GENE_SETS_DIR      <- "/home/rachele/SNVs/results_final_2026_march"   # folder with csv/tsv gene set files
# MANIFEST_PATH      <- "/home/rachele/SNVs/results_pasteur_tsv_With_FEATURE/manifest_correct.tsv"
MANIFEST_PATH      <- "/home/rachele/SNVs/results_final_2026_march/manifest_clean_4Groups.tsv"
# QC_UHR_PATH        <- "~/UHR_NA_SAMPLE_IDS.csv"
# QC_LOWQUAL_PATH    <- "~/Low_Quality_SAMPLES_AND_UHR_NA.tsv"

# Output base
OUTPUT_BASE        <- "/home/rachele/SNVs/results_final_2026_march"

# Samples to eliminate (hardcoded from original executed files)
# ELIMINATE_SAMPLES  <- c("S36913", "S36916", "S36917", "S36919", "S36959")


# ============================================================
# 3. LOAD RAW DATA
# ============================================================

cat("=== Loading variant data ===\n")
raw_variant_data <- load_variant_data(VARIANT_DIR)

# Rename from ugly VCF filenames to short names
# (matches what was done manually in the executed files)
Missense_canonical <- raw_variant_data[["missense_MPC_AM_canonical.vcf.gz"]]
PTV_canonic       <- raw_variant_data[["PTV_HC_canonical.vcf.gz"]]
rm(raw_variant_data)

Missense_canonical$SAMPLES <- gsub("_pool", "", Missense_canonical$SAMPLES)
PTV_canonic$SAMPLES <- gsub("_pool", "", PTV_canonic$SAMPLES)

cat("Variant data loaded.\n\n")

# ============================================================
# 4. APPLY CONSTRAINT FILTER (lof.oe_ci.upper < 0.6)
# ============================================================

cat("=== Applying constraint filter ===\n")
# gnomad.v4.1.constraint_metrics must already be loaded in your environment
gnomad.v4.1.constraint_metrics <- read.delim("~/SNVs/results_final_2026_march/gnomad.v4.1.constraint_metrics.tsv")
Missense_canonical <- apply_constraint_single(Missense_canonical, gnomad.v4.1.constraint_metrics)
PTV_canonic        <- apply_constraint_single(PTV_canonic,        gnomad.v4.1.constraint_metrics)
cat("Constraint filter applied.\n\n")


PTV_canonic <- PTV_canonic%>% filter(AF < 0.05)
# ============================================================
# 5. LOAD GENE SETS FROM FOLDER
# ============================================================

cat("=== Loading gene sets ===\n")
# load_gene_sets reads all csv/tsv/txt from the folder
# and assigns genes_<filename> to global env automatically
gene_sets <- load_gene_sets(GENE_SETS_DIR, gene_column = "Gene")

# Build the gene_lists named list used by all downstream functions
# Adjust names to match what is in your gene sets folder
# gene_lists <- list(
#   brain      = unique(gene_sets[["brain_gene_consensus_filtered_consensus_no_pitular"]]),
#   brain_ntpm = unique(gene_sets[["brain_gene_consensus_ntm_consensus_no_pitular"]]),
#   schema_pval = unique(gene_sets[["SCHEMA_pval"]]),      # adjust key to match your filename
#   schema_qval = unique(gene_sets[["SCHEMA_qval"]]),
#   bipolar     = unique(gene_sets[["BipEx_Bipolar"]]),
#   gwas        = unique(gene_sets[["GWAS_120"]])
# )

gene_lists <- list(
  brain      = unique(gene_sets[["brain_gene_consensus_filtered_consensus_no_pitular"]]),
  brain_ntpm = unique(gene_sets[["brain_gene_consensus_ntm_consensus_no_pitular"]])
)
cat("Gene sets loaded:", paste(names(gene_lists), collapse = ", "), "\n\n")

# ============================================================
# 6. LOAD AND CLEAN MANIFEST
# ============================================================

cat("=== Loading manifest ===\n")
manifest_raw <- read.delim(MANIFEST_PATH)

# Rename columns to match what refactored functions expect:
# sample_id and Status
# Check if 'Sequencing_number' exists before renaming
if ("Sequencing_number" %in% colnames(manifest_raw)) {
  manifest_clean <- manifest_raw %>%
    rename(sample_id = Sequencing_number) %>%      # rename if exists
    filter(!grepl("UHR_NA", Status)) %>%           # remove UHR_NA group entries
    filter(!grepl("UHR-NA", Status)) %>%           # remove UHR_NA group entries
    filter(sample_id != "")                         # remove empty sample_ids
} else {
  manifest_clean <- manifest_raw %>%
    filter(!grepl("UHR_NA", Status)) %>%           # remove UHR_NA group entries
    filter(!grepl("UHR-NA", Status)) %>%           # remove UHR_NA group entries
    filter(sample_id != "")                         # remove empty sample_ids
}


cat("Manifest cleaned. Sample counts per group:\n")
print(table(manifest_clean$Status))
cat("\n")

# ============================================================
# 7. VARIANT DATASETS LIST (used across both analyses)
# ============================================================

all_variant_datasets <- list(
  Missense_canonical = Missense_canonical,
  PTV_canonic        = PTV_canonic
)


# ============================================================
# ============================================================
# SECTION B: 4-GROUP ANALYSIS (SCZ / BD / Converter / Non_Converter)
# ============================================================
# ============================================================

cat("\n\n========================================\n")
cat("SECTION B: 4-GROUP ANALYSIS\n")
cat("========================================\n\n")

# Recode manifest for 4 groups
manifest_4group <- manifest_clean %>%
  mutate(Status = case_when(
    Status == "FEP-SCZ"       ~ "SCZ",
    Status == "FEP-BD"        ~ "BD",
    Status == "Converter"     ~ "Converter",
    Status == "Non_Converter" ~ "Non_Converter",
    TRUE                      ~ Status
  )) %>%
  filter(Status %in% c("SCZ", "BD", "Converter", "Non_Converter"))

cat("4-group manifest counts:\n")
print(table(manifest_4group$Status))
cat("\n")

# --------------------------------------------------
# B1. Process variant data (adds sample labels)
# Only PTV_canonic used for 4-group (matches executed file)
# --------------------------------------------------
cat("--- Processing variant data for 4-group ---\n")
datasets_4group <- list(PTV_canonic = PTV_canonic)

processed_4group <- process_variant_data(
  variant_data = datasets_4group,
  manifest     = manifest_4group
)


# --------------------------------------------------
# B2. NOT PRIVATE — Pairwise Fisher + Kruskal
# --------------------------------------------------
cat("\n--- B2: NOT PRIVATE ---\n")
output_4group_notprivate <- file.path(OUTPUT_BASE, "group4", "not_private")
dir.create(output_4group_notprivate, recursive = TRUE, showWarnings = FALSE)

# Fisher count tables (4 groups → pairwise Fisher automatically)
cat("Creating Fisher count tables (not private)...\n")
fisher_counts_4group_notprivate <- list()
for (ds_name in names(processed_4group)) {
  fisher_counts_4group_notprivate[[ds_name]] <- create_count_dataframes(
    df          = processed_4group[[ds_name]],
    name        = ds_name,
    gene_lists  = list(
      brain      = gene_lists$brain,
      brain_ntpm = gene_lists$brain_ntpm
    ),  # 4-group executed file only used brain gene sets for Fisher
    manifest_df = manifest_4group,
    private     = FALSE,
    row_selection = "unique_individuals"
  )
}

# Pairwise Fisher (auto-detected because >2 groups in manifest)
cat("Running pairwise Fisher tests (not private)...\n")
fisher_results_4group_notprivate <- run_fisher_with_manifest_enhanced(
  count_df    = fisher_counts_4group_notprivate,
  manifest    = manifest_4group,
  fdr_scope   = "both",
  clean_names = TRUE
)

# Rates
cat("Calculating rates (not private)...\n")
rates_4group_notprivate <- calculate_rates(
  data_df    = processed_4group,
  manifest   = manifest_4group,
  gene_lists = list(brain = gene_lists$brain, brain_ntpm = gene_lists$brain_ntpm),
  rows = "all",
  private = FALSE,
  row_selection = "variants"
)

# Kruskal count tables
cat("Creating Kruskal count tables (not private)...\n")
kruskal_counts_4group_notprivate <- list()
kruskal_counts_4group_notprivate <- create_wilcoxon_kruskal_tables(
     processed_data    = processed_4group,
    gene_lists  = list(
      brain      = gene_lists$brain,
      brain_ntpm = gene_lists$brain_ntpm
    ),
    manifest_data= manifest_4group,
    private     = FALSE,
    pure        = FALSE
  )



# Kruskal-Wallis + Dunn post-hoc (auto-detected because >2 groups)
cat("Running Kruskal-Wallis tests (not private)...\n")
kruskal_results_4group_notprivate <- run_wilcoxon_kruskal_analysis(
  data_df          = kruskal_counts_4group_notprivate$complete_datasets,
  fdr_scope        = "both",
  run_dunn_posthoc = TRUE
)
kruskal_df_notprivate <- results_to_dataframe(kruskal_results_4group_notprivate)

# Plots
cat("Creating plots (not private)...\n")
fisher_plots_4group_notprivate <- create_comprehensive_visualization(
  results_df        = fisher_results_4group_notprivate,
  plot_types        = c("volcano", "parallel", "dot", "manhattan", "forest"),
  significant_only  = FALSE,
  save_plots        = TRUE,
  output_dir        = output_4group_notprivate,
  prefix            = "fisher_4group_notprivate"
)

kruskal_plots_notprivate <- create_variant_plots(
  plot_data        = kruskal_counts_4group_notprivate$complete_datasets,
  results_table    = kruskal_df_notprivate,
  significant_only = TRUE
)
save_variant_plots(
  kruskal_plots_notprivate,
  prefix     = "kruskal_4group_notprivate",
  output_dir = output_4group_notprivate
)

# Save tables
cat("Saving tables (not private)...\n")
setwd(output_4group_notprivate)
write.csv(fisher_results_4group_notprivate, "fisher_results_4group_notprivate.csv",   row.names = FALSE)
write.csv(rates_4group_notprivate,          "rates_4group_notprivate.csv",            row.names = FALSE)
write.csv(kruskal_df_notprivate,            "kruskal_results_4group_notprivate.csv",  row.names = FALSE)
for (ds_name in names(fisher_counts_4group_notprivate)) {
  write.csv(fisher_counts_4group_notprivate[[ds_name]],
            paste0("count_table_", ds_name, "_notprivate.csv"), row.names = FALSE)
}
cat("Not-private 4-group done.\n\n")


# --------------------------------------------------
# B3. Gene lists Excel (4-group)
# --------------------------------------------------
cat("--- B3: Gene lists Excel (4-group) ---\n")
output_4group_gene <- file.path(OUTPUT_BASE, "group4", "gene")
dir.create(output_4group_gene, recursive = TRUE, showWarnings = FALSE)

for (ds_name in names(processed_4group)) {
  create_gene_lists_excel(
    df           = processed_4group[[ds_name]],
    gene_lists   = gene_lists,
    output_dir   = output_4group_gene,
    dataset_name = ds_name,
    group_column = "sample_label",
    gene_column  = "SYMBOL"
  )
  create_gene_list_summary(
    df            = processed_4group[[ds_name]],
    gene_lists    = gene_lists,
    output_dir    = output_4group_gene,
    dataset_name  = ds_name,
    group_column  = "sample_label",
    gene_column   = "SYMBOL",
    sample_column = "SAMPLES"
  )
  create_enrichment_excel(
    datasets = processed_4group,
    output_dir   = output_4group_gene, 
    gene_lists   = gene_lists
  )
}
cat("Gene lists Excel done.\n\n")


# ============================================================
# DONE
# ============================================================
cat("\n\n========================================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("========================================\n")
cat("Outputs saved under:", OUTPUT_BASE, "\n")
cat("  group4/not_private/  — 4-group pairwise Fisher + Kruskal\n")
cat("  group4/gene/         — Gene list Excel files\n")