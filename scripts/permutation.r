# ============================================================
# PERMUTATION TESTS FOR GENE-SET BURDEN
# permutation_tests.r
# ============================================================
#
# PURPOSE:
#   Build empirical null distributions by shuffling sample Status
#   labels, for both Fisher and Wilcoxon/Kruskal gene-set burden tests.
#
# IMPORTANT — NO DUPLICATED LOGIC:
#   This script only contains the shuffling loop.
#   All statistical functions are the ones already in statistics functions.r:
#     - perform_fisher_core()          → Fisher 2x2 table
#     - normalize_manifest()           → manifest cleaning
#     - extract_groups_from_manifest() → group names + totals
#     - wilcox.test() / kruskal.test() → base R, called directly
#       (the existing wrappers like run_wilcoxon_kruskal_analysis are too
#        heavy for a tight permutation loop — they apply FDR, reshape, log)
#
# INPUT (both functions):
#   complete_datasets — from create_wilcoxon_kruskal_tables$complete_datasets
#   This already has per-individual carrier counts per gene set, so no
#   QC files or variant data need to be re-read inside the loop.
#
# USAGE (drop-in after your existing analysis blocks in 2group.r / 4groups.r):
#
#   perm_results <- run_permutation_burden_tests(
#     complete_datasets = wilcoxon_counts_2group_notprivate$complete_datasets,
#     manifest          = manifest_2group,
#     gene_lists        = gene_lists,
#     n_perm            = 1000   # use 10000 for publication
#   )
#   save_permutation_results(perm_results,
#                             output_dir = output_2group_notprivate,
#                             prefix     = "permutation_notprivate")
#
# ============================================================

library(dplyr)
library(purrr)
library(tibble)


# ============================================================
# INTERNAL HELPERS
# ============================================================

# Return ALL Number_of_* columns from a complete_dataset df.
# This includes total burden (Number_of_pure_Variants) AND all gene-set
# columns (Number_of_pure_brain_Variants etc.) so both are permuted.
.all_number_cols <- function(df) {
  grep("^Number_of", names(df), value = TRUE)
}

# When gene_lists is provided, find the specific gene-set columns.
# Handles both "Number_of_<gs>_Variants" and "Number_of_pure_<gs>_Variants".
.find_gene_set_cols <- function(df, gene_set_names) {
  vapply(gene_set_names, function(gs) {
    pattern <- paste0("Number_of(?:_pure)?_", gs, "_Variants")
    matches <- grep(pattern, names(df), value = TRUE, perl = TRUE)
    if (length(matches) == 0) NA_character_ else matches[1]
  }, character(1))
}


# ============================================================
# SECTION 1: FISHER PERMUTATION
# ============================================================

#' Permutation test for Fisher gene-set burden
#'
#' Derives binary carrier status from complete_datasets
#' (carrier = Number_of_X_Variants > 0, mirroring the "unique_individuals"
#' row in create_count_dataframes), shuffles Status labels n_perm times,
#' and calls perform_fisher_core() from statistics functions.r each time.
#'
#' @param complete_datasets Named list from create_wilcoxon_kruskal_tables$complete_datasets.
#' @param manifest Manifest with sample_id and Status.
#' @param n_perm Number of permutations (1000 = exploration, 10000 = publication).
#' @param seed Random seed.
#' @param verbose Print per-test progress.
#'
#' @return Dataframe: Dataset, metric, Group1, Group2,
#'   n_carriers_g1, n_carriers_g2, n_total_g1, n_total_g2,
#'   observed_OR, observed_p, perm_p, perm_p_adj, n_perm
#' @export
run_permutation_fisher <- function(complete_datasets,
                                   manifest,
                                   n_perm  = 1000,
                                   seed    = 42,
                                   verbose = TRUE) {

  set.seed(seed)

  # Reuse existing manifest utilities from statistics functions.r
  manifest   <- normalize_manifest(manifest)
  group_info <- extract_groups_from_manifest(manifest)
  groups  <- group_info$groups
  totals  <- group_info$counts   # named vector: group -> total N
  pairs   <- combn(groups, 2, simplify = FALSE)

  all_results <- list()

  for (ds_name in names(complete_datasets)) {
    if (verbose) cat(sprintf("\n[Fisher Permutation] Dataset: %s\n", ds_name))

    df       <- complete_datasets[[ds_name]]
    all_cols <- .all_number_cols(df)   # ALL Number_of_* cols: total + gene sets

    if (length(all_cols) == 0) {
      warning(sprintf("Dataset %s: no Number_of_ columns found, skipping.", ds_name))
      next
    }

    # Build binary carrier matrix once — reused across all permutations.
    # carrier = 1 if the individual has >= 1 variant in that column.
    # Applies to total burden (Number_of_pure_Variants) AND each gene set.
    carrier_mat <- df %>%
      select(sample_id, Status) %>%
      bind_cols(
        map_dfc(all_cols, function(col) {
          tibble(!!col := as.integer(df[[col]] > 0))
        })
      )

    status_vec <- carrier_mat$Status

    for (col in all_cols) {
      carrier_vec <- carrier_mat[[col]]

      for (pair in pairs) {
        g1 <- pair[1]; g2 <- pair[2]

        mask         <- status_vec %in% c(g1, g2)
        pair_status  <- status_vec[mask]
        pair_carrier <- carrier_vec[mask]

        # Observed — call existing perform_fisher_core() from statistics functions.r
        a_obs <- sum(pair_carrier[pair_status == g1])
        c_obs <- sum(pair_carrier[pair_status == g2])
        obs   <- perform_fisher_core(a_obs, c_obs, totals[[g1]], totals[[g2]])

        if (is.null(obs)) {
          if (verbose) cat(sprintf("  Skipping %s [%s vs %s]: degenerate table\n", col, g1, g2))
          next
        }

        obs_stat <- -log10(obs$p_value + 1e-300)

        # Permutation loop — only shuffling, reuse perform_fisher_core()
        perm_stats <- vapply(seq_len(n_perm), function(i) {
          shuffled <- sample(pair_status)
          a_perm   <- sum(pair_carrier[shuffled == g1])
          c_perm   <- sum(pair_carrier[shuffled == g2])
          res      <- perform_fisher_core(a_perm, c_perm, totals[[g1]], totals[[g2]])
          if (is.null(res)) 0 else -log10(res$p_value + 1e-300)
        }, numeric(1))

        # Empirical p (Phipson-Smyth +1 correction for finite permutations)
        perm_p <- (sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)

        if (verbose) {
          cat(sprintf("  %s [%s vs %s]: OR=%.2f, obs_p=%.4f, perm_p=%.4f\n",
                      col, g1, g2, obs$OR, obs$p_value, perm_p))
        }

        all_results[[length(all_results) + 1]] <- tibble(
          Dataset       = ds_name,
          metric        = col,
          Group1        = g1,
          Group2        = g2,
          n_carriers_g1 = a_obs,
          n_carriers_g2 = c_obs,
          n_total_g1    = totals[[g1]],
          n_total_g2    = totals[[g2]],
          observed_OR   = obs$OR,
          observed_p    = obs$p_value,
          perm_p        = perm_p,
          n_perm        = n_perm
        )
      }
    }
  }

  if (length(all_results) == 0) { warning("No Fisher permutation results produced."); return(NULL) }

  bind_rows(all_results) %>%
    mutate(perm_p_adj = p.adjust(perm_p, method = "fdr")) %>%
    arrange(perm_p)
}


# ============================================================
# SECTION 2: WILCOXON / KRUSKAL PERMUTATION
# ============================================================

#' Permutation test for Wilcoxon/Kruskal gene-set burden
#'
#' Shuffles Status labels in complete_datasets and calls base R
#' wilcox.test() or kruskal.test() directly. The existing wrappers
#' (run_wilcoxon_kruskal_analysis etc.) are intentionally not used here
#' because they apply FDR, reshape, and log — all wasteful in a loop.
#'
#' @param complete_datasets Named list from create_wilcoxon_kruskal_tables$complete_datasets.
#' @param n_perm Number of permutations.
#' @param seed Random seed.
#' @param verbose Print per-test progress.
#'
#' @return Dataframe: Dataset, metric, test_type,
#'   Group1, Group2 (if 2-group), observed_stat, observed_p,
#'   perm_p, perm_p_adj, n_perm
#' @export
run_permutation_wilcoxon <- function(complete_datasets,
                                     n_perm  = 1000,
                                     seed    = 42,
                                     verbose = TRUE) {

  set.seed(seed)

  all_results <- list()

  for (ds_name in names(complete_datasets)) {
    if (verbose) cat(sprintf("\n[Wilcoxon/Kruskal Permutation] Dataset: %s\n", ds_name))

    df       <- complete_datasets[[ds_name]]
    groups   <- unique(df$Status[!is.na(df$Status)])
    n_groups <- length(groups)

    if (n_groups < 2) { warning(sprintf("Dataset %s: <2 groups, skipping.", ds_name)); next }

    # Test ALL Number_of_* columns: total burden + every gene set
    test_cols <- .all_number_cols(df)

    if (length(test_cols) == 0) { warning(sprintf("Dataset %s: no columns found, skipping.", ds_name)); next }

    status_vec <- df$Status

    for (col in test_cols) {
      values <- df[[col]]

      if (n_groups == 2) {
        # ---- Wilcoxon ----
        obs <- tryCatch(
          wilcox.test(values[status_vec == groups[1]],
                      values[status_vec == groups[2]], exact = FALSE),
          error = function(e) NULL
        )
        if (is.null(obs)) next

        obs_stat <- -log10(obs$p.value + 1e-300)

        perm_stats <- vapply(seq_len(n_perm), function(i) {
          sh  <- sample(status_vec)
          res <- tryCatch(
            wilcox.test(values[sh == groups[1]], values[sh == groups[2]], exact = FALSE),
            error = function(e) NULL
          )
          if (is.null(res)) 0 else -log10(res$p.value + 1e-300)
        }, numeric(1))

        perm_p <- (sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)

        if (verbose) cat(sprintf("  %s [%s vs %s]: obs_p=%.4f, perm_p=%.4f\n",
                                  col, groups[1], groups[2], obs$p.value, perm_p))

        all_results[[length(all_results) + 1]] <- tibble(
          Dataset       = ds_name,
          metric        = col,
          test_type     = "Wilcoxon",
          Group1        = groups[1],
          Group2        = groups[2],
          observed_stat = obs$statistic,
          observed_p    = obs$p.value,
          perm_p        = perm_p,
          n_perm        = n_perm
        )

      } else {
        # ---- Kruskal-Wallis ----
        obs <- tryCatch(kruskal.test(values ~ status_vec), error = function(e) NULL)
        if (is.null(obs)) next

        obs_stat <- -log10(obs$p.value + 1e-300)

        perm_stats <- vapply(seq_len(n_perm), function(i) {
          sh  <- sample(status_vec)
          res <- tryCatch(kruskal.test(values ~ sh), error = function(e) NULL)
          if (is.null(res)) 0 else -log10(res$p.value + 1e-300)
        }, numeric(1))

        perm_p <- (sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)

        if (verbose) cat(sprintf("  %s [Kruskal]: obs_p=%.4f, perm_p=%.4f\n",
                                  col, obs$p.value, perm_p))

        all_results[[length(all_results) + 1]] <- tibble(
          Dataset       = ds_name,
          metric        = col,
          test_type     = "Kruskal-Wallis",
          Group1        = NA_character_,
          Group2        = NA_character_,
          observed_stat = obs$statistic,
          observed_p    = obs$p.value,
          perm_p        = perm_p,
          n_perm        = n_perm
        )
      }
    }
  }

  if (length(all_results) == 0) { warning("No Wilcoxon permutation results produced."); return(NULL) }

  bind_rows(all_results) %>%
    mutate(perm_p_adj = p.adjust(perm_p, method = "fdr")) %>%
    arrange(perm_p)
}


# ============================================================
# SECTION 3: COMBINED WRAPPER
# ============================================================

#' Run both Fisher and Wilcoxon/Kruskal permutation tests
#'
#' @param complete_datasets From create_wilcoxon_kruskal_tables$complete_datasets.
#' @param manifest Manifest with sample_id and Status.
#' @param gene_lists Named list of gene set vectors.
#' @param n_perm Number of permutations.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return Named list: $fisher_perm and $wilcoxon_perm dataframes.
#' @export
run_permutation_burden_tests <- function(complete_datasets,
                                          manifest,
                                          gene_lists,
                                          n_perm  = 1000,
                                          seed    = 42,
                                          verbose = TRUE) {

  cat("============================================================\n")
  cat("PERMUTATION BURDEN TESTS\n")
  cat(sprintf("Permutations: %d | Datasets: %s | Gene sets: %s\n",
              n_perm,
              paste(names(complete_datasets), collapse = ", "),
              paste(names(gene_lists), collapse = ", ")))
  cat("============================================================\n\n")

  fisher_perm <- tryCatch(
    run_permutation_fisher(complete_datasets, manifest, n_perm, seed, verbose),
    error = function(e) { warning("Fisher permutation failed: ", e$message); NULL }
  )

  wilcoxon_perm <- tryCatch(
    run_permutation_wilcoxon(complete_datasets, n_perm, seed, verbose),
    error = function(e) { warning("Wilcoxon permutation failed: ", e$message); NULL }
  )

  cat("\n============================================================\n")
  if (!is.null(fisher_perm))   cat(sprintf("Fisher results   : %d rows\n", nrow(fisher_perm)))
  if (!is.null(wilcoxon_perm)) cat(sprintf("Wilcoxon results : %d rows\n", nrow(wilcoxon_perm)))
  cat("============================================================\n")

  list(fisher_perm = fisher_perm, wilcoxon_perm = wilcoxon_perm)
}


# ============================================================
# SECTION 4: SAVE
# ============================================================

#' Save permutation results to CSV
#' @export
save_permutation_results <- function(perm_results, output_dir, prefix = "permutation") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!is.null(perm_results$fisher_perm)) {
    path <- file.path(output_dir, paste0(prefix, "_fisher.csv"))
    write.csv(perm_results$fisher_perm, path, row.names = FALSE)
    cat(sprintf("Saved: %s\n", path))
  }
  if (!is.null(perm_results$wilcoxon_perm)) {
    path <- file.path(output_dir, paste0(prefix, "_wilcoxon.csv"))
    write.csv(perm_results$wilcoxon_perm, path, row.names = FALSE)
    cat(sprintf("Saved: %s\n", path))
  }
}


# ============================================================
# USAGE EXAMPLES
# ============================================================
#
# ---- 2-GROUP (slot in after your existing Wilcoxon block) ----
#
  perm_notprivate <- run_permutation_burden_tests(
    complete_datasets = wilcoxon_counts_2group_notprivate$complete_datasets,
    manifest          = manifest_2group,
    gene_lists        = gene_lists,
    n_perm            = 10000
  )
  save_permutation_results(perm_notprivate, output_2group_notprivate, "permutation_notprivate")

#   perm_private <- run_permutation_burden_tests(
#     complete_datasets = wilcoxon_counts_2group_private$complete_datasets,
#     manifest          = manifest_2group,
#     gene_lists        = gene_lists,
#     n_perm            = 1000
#   )
#   save_permutation_results(perm_private, output_2group_private, "permutation_private")
#
#
# ---- 4-GROUP (slot in after your existing Kruskal block) ----
#
#   perm_4group <- run_permutation_burden_tests(
#     complete_datasets = kruskal_counts_4group_notprivate$complete_datasets,
#     manifest          = manifest_4group,
#     gene_lists        = list(brain = gene_lists$brain, brain_ntpm = gene_lists$brain_ntpm),
#     n_perm            = 1000
#   )
#   save_permutation_results(perm_4group, output_4group_notprivate, "permutation_4group_notprivate")
#
#
# ---- READING RESULTS ----
#
#   perm_notprivate$fisher_perm %>%
#     filter(perm_p_adj < 0.05) %>%
#     select(Dataset, Gene_Set, Group1, Group2,
#            n_carriers_g1, n_carriers_g2, observed_OR, observed_p, perm_p, perm_p_adj)
#
#   perm_notprivate$wilcoxon_perm %>%
#     filter(perm_p_adj < 0.05) %>%
#     select(Dataset, Gene_Set, test_type, observed_p, perm_p, perm_p_adj)
#
# ============================================================