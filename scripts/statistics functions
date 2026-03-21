#STatitsics

library(purrr)
library(dplyr)
library(tidyr)
library(broom)
library(dunn.test)



# FISHER TEST 


#' Run Fisher tests with manifest - enhanced version
#'
#' @param count_df Count dataframe or list of dataframes
#' @param manifest Manifest dataframe
#' @param rows Which rows to test
#' @param fdr_scope "global", "within", or "both" for FDR correction
#' @param clean_names Clean row names (TRUE/FALSE)
#'
#' @return Results with metadata
#' @export
run_fisher_with_manifest_enhanced <- function(count_df, manifest, rows = "all", 
                                              fdr_scope = "both", clean_names = TRUE) {
    
    # Handle both single dataframe and list
    if (is.list(count_df) && !is.data.frame(count_df)) {
        # Multiple datasets
        results <- run_fisher_on_multiple_datasets(
            count_list = count_df,
            manifest = manifest,
            rows = rows,
            fdr_scope = fdr_scope,
            clean_names = clean_names
        )
    } else {
        # Single dataset
        # Extract groups and totals from manifest
        manifest_info <- extract_groups_from_manifest(manifest)
        groups <- manifest_info$groups
        totals <- manifest_info$counts
        
        # Check that all groups in count_df exist in manifest
        count_groups <- detect_groups_from_counts(count_df)
        missing_in_manifest <- setdiff(count_groups, names(totals))
        if (length(missing_in_manifest) > 0) {
            stop(sprintf("Groups in count_df not found in manifest: %s\nAvailable in manifest: %s",
                        paste(missing_in_manifest, collapse = ", "),
                        paste(names(totals), collapse = ", ")))
        }
        
        # Ensure count_df has all manifest groups (fill missing with 0)
        missing_in_counts <- setdiff(names(totals), count_groups)
        if (length(missing_in_counts) > 0) {
            for (group in missing_in_counts) {
                count_df[[group]] <- 0
            }
            warning(sprintf("Added missing groups with 0 counts: %s",
                            paste(missing_in_counts, collapse = ", ")))
        }
        
        groups <- detect_groups_from_counts(count_df)
        
        # Run Fisher tests
        if (length(groups) == 2) {
            results <- run_fisher_2group(count_df, totals, rows)
        } else {
            results <- run_fisher_pairwise(count_df, totals, rows)
        }
        
        # Apply appropriate FDR corrections for single dataset
        results <- apply_fdr_corrections_single(results, fdr_scope)
            # Clean names if requested
        if (clean_names && !is.null(results) && nrow(results) > 0) {
            results <- clean_row_names(results)
        }
        
        # Add metadata
        attr(results, "manifest_groups") <- groups
        attr(results, "group_totals") <- totals
        attr(results, "fdr_scope") <- fdr_scope
        attr(results, "is_single_dataset") <- TRUE
    }
    

    
    return(results)
}

#' Apply FDR corrections for single dataset
#'
#' @param results Results dataframe
#' @param fdr_scope "global", "within", or "both"
#'
#' @return Dataframe with appropriate FDR corrections
apply_fdr_corrections_single <- function(results, fdr_scope = "both") {
  if (is.null(results) || nrow(results) == 0) {
    return(results)
  }
  
  # For single dataset, p_adj (from run_fisher_2group/pairwise) is the within-dataset FDR
  if ("p_adj" %in% colnames(results)) {
    if (fdr_scope %in% c("within", "both")) {
      results$p_adj_within <- results$p_adj
    }
    
    if (fdr_scope %in% c("global", "both")) {
      # For single dataset, global and within are the same
      results$p_adj_global <- results$p_adj
    }
    
    # Remove the original p_adj to avoid confusion
    results <- results %>% select(-p_adj)
  }
  
  return(results)
}

#' Run Fisher tests on multiple count datasets
#'
#' @param count_list Named list of count dataframes
#' @param manifest Manifest dataframe
#' @param rows Which rows to test (default: "all")
#' @param fdr_scope "global", "within", or "both" for FDR correction
#' @param clean_names Clean row names (TRUE/FALSE)
#' 
#' @return List of results with metadata
#' @export
run_fisher_on_multiple_datasets <- function(count_list, manifest, rows = "all", 
                                           fdr_scope = "both", clean_names = TRUE) {
  
  if (!is.list(count_list)) {
    stop("count_list must be a named list of dataframes")
  }
  
  if (is.null(names(count_list))) {
    names(count_list) <- paste0("Dataset_", seq_along(count_list))
  }
  
  # Extract manifest info once
  manifest_info <- extract_groups_from_manifest(manifest)
  groups <- manifest_info$groups
  totals <- manifest_info$counts
  
  # Process each dataset
  all_results <- list()
  
  for (dataset_name in names(count_list)) {
    cat(sprintf("Processing dataset: %s\n", dataset_name))
    
    tryCatch({
      # Prepare count_df (check and fill missing groups)
      count_df <- count_list[[dataset_name]]
      
      # Check that all groups in count_df exist in manifest
      count_groups <- detect_groups_from_counts(count_df)
      missing_in_manifest <- setdiff(count_groups, names(totals))
      if (length(missing_in_manifest) > 0) {
        stop(sprintf("Groups in count_df not found in manifest: %s\nAvailable in manifest: %s",
                     paste(missing_in_manifest, collapse = ", "),
                     paste(names(totals), collapse = ", ")))
      }
      
      # Ensure count_df has all manifest groups (fill missing with 0)
      missing_in_counts <- setdiff(names(totals), count_groups)
      if (length(missing_in_counts) > 0) {
        for (group in missing_in_counts) {
          count_df[[group]] <- 0
        }
        warning(sprintf("Dataset %s: Added missing groups with 0 counts: %s",
                        dataset_name, paste(missing_in_counts, collapse = ", ")))
      }
      
      # Run Fisher tests
      if (length(groups) == 2) {
        results <- run_fisher_2group(count_df = count_df, totals = totals, rows = rows)
      } else {
        results <- run_fisher_pairwise(count_df = count_df, totals = totals, rows = rows)
      }
      
      # Add dataset name
      results$Dataset <- dataset_name
      
      # Store results
      all_results[[dataset_name]] <- results
      
    }, error = function(e) {
      warning(sprintf("Error processing dataset %s: %s", dataset_name, e$message))
    })
  }
  
  # Combine all results
  if (length(all_results) == 0) {
    stop("No datasets were successfully processed")
  }
  
  combined_results <- bind_rows(all_results)
  
  # Apply FDR corrections
  if (fdr_scope %in% c("within", "both")) {
    combined_results <- apply_fdr_corrections(combined_results, test_col = "Dataset")
  }
  
  if (fdr_scope %in% c("global", "both")) {
    # Ensure global FDR is always applied (it's already in apply_fdr_corrections)
    if (!"p_adj_global" %in% colnames(combined_results)) {
      combined_results$p_adj_global <- p.adjust(combined_results$p_value, method = "fdr")
    }
  }
  
  # Clean names if requested
  if (clean_names) {
    combined_results <- clean_row_names(combined_results)
  }
  
  # Add comprehensive metadata
  attr(combined_results, "datasets") <- names(count_list)
  attr(combined_results, "n_datasets") <- length(count_list)
  attr(combined_results, "manifest_groups") <- manifest_info$groups
  attr(combined_results, "group_totals") <- manifest_info$counts
  attr(combined_results, "fdr_scope") <- fdr_scope
  attr(combined_results, "is_single_dataset") <- FALSE
  
  return(combined_results)
}

#' Clean row names in results
#' 
#' @param results_df Results dataframe with Row column
#' @param keep_original Keep original row names (TRUE/FALSE)
#' 
#' @return Cleaned dataframe
#' @export
clean_row_names <- function(results_df, keep_original = TRUE) {
    if (is.null(results_df) || nrow(results_df) == 0) {
        return(results_df)
    }
    
    # If Row column doesn't exist, return as is
    if (!"Row" %in% colnames(results_df)) {
        return(results_df)
    }
    
    # Create clean version
    clean_names <- results_df$Row
    
    # Remove everything before and including "with"
    clean_names <- gsub(".*with ", "", clean_names)
    
    # Remove the word before "variants" but keep "variants"
    clean_names <- gsub("\\S+ variants", "variants", clean_names)
    
    # Final cleanup
    clean_names <- trimws(clean_names)
    
    # Update the dataframe
    if (keep_original) {
      results_df$Row_Original <- results_df$Row
      results_df$Row <- clean_names
    } else {
      results_df$Row <- clean_names
    }

    return(results_df)
}

detect_groups_from_counts <- function(count_df) {
  group_cols <- colnames(count_df)[-1]
  if (length(group_cols) < 2) {
    stop("At least 2 groups are required for statistical testing")
  }
  group_cols
}

select_rows <- function(df, rows = "all") {
  if (identical(rows, "all")) {
    return(df)
  }
  
  missing <- setdiff(rows, df[[1]])
  if (length(missing) > 0) {
    stop("Unknown row names: ", paste(missing, collapse = ", "))
  }
  
  df[df[[1]] %in% rows, , drop = FALSE]
}

perform_fisher_core <- function(a, c, total_a, total_c) {
  b <- total_a - a
  d <- total_c - c
  
  if (any(c(a, b, c, d) < 0)) {
    return(NULL)
  }
  
  tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  test <- fisher.test(tbl)
  
  tibble::tibble(
    OR = unname(test$estimate),
    CI_low = test$conf.int[1],
    CI_high = test$conf.int[2],
    p_value = test$p.value
  )
}

run_fisher_2group <- function(count_df, totals, rows = "all") {
  
  count_df <- select_rows(count_df, rows)
  groups <- detect_groups_from_counts(count_df)
  
  if (length(groups) != 2) {
    stop("Fisher test requires exactly 2 groups")
  }
  
  g1 <- groups[1]
  g2 <- groups[2]
  
  purrr::map_dfr(seq_len(nrow(count_df)), function(i) {
    row <- count_df[i, ]
    
    res <- perform_fisher_core(
      a = row[[g1]],
      c = row[[g2]],
      total_a = totals[[g1]],
      total_c = totals[[g2]]
    )
    
    if (is.null(res)) return(NULL)
    
    dplyr::bind_cols(
      tibble::tibble(
        Row = row[[1]],
        Group1 = g1,
        Group2 = g2
      ),
      res
    )
  }) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "fdr"))
}

run_fisher_pairwise <- function(count_df, totals, rows = "all") {
  
  count_df <- select_rows(count_df, rows)
  groups <- detect_groups_from_counts(count_df)
  
  pairs <- combn(groups, 2, simplify = FALSE)
  
  purrr::map_dfr(seq_len(nrow(count_df)), function(i) {
    row <- count_df[i, ]
    
    purrr::map_dfr(pairs, function(p) {
      g1 <- p[1]; g2 <- p[2]
      
      res <- perform_fisher_core(
        a = row[[g1]],
        c = row[[g2]],
        total_a = totals[[g1]],
        total_c = totals[[g2]]
      )
      
      if (is.null(res)) return(NULL)
      
      dplyr::bind_cols(
        tibble::tibble(
          Row = row[[1]],
          Group1 = g1,
          Group2 = g2
        ),
        res
      )
    })
  }) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "fdr"))
}

#' Apply multiple FDR corrections
#'
#' @param results_df Results dataframe
#' @param test_col Column indicating test/dataset (for within-test FDR)
#' 
#' @return Dataframe with multiple FDR corrections
#' @export
apply_fdr_corrections <- function(results_df, test_col = NULL) {
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(results_df)
  }
  
  # Global FDR across all tests
  results_df$p_adj_global <- p.adjust(results_df$p_value, method = "fdr")
  
  # Within-test FDR (if test column provided)
  if (!is.null(test_col) && test_col %in% colnames(results_df)) {
    results_df <- results_df %>%
      group_by(!!sym(test_col)) %>%
      mutate(p_adj_within = p.adjust(p_value, method = "fdr")) %>%
      ungroup()
  }
  
  # Format p-values for display
  results_df$p_value_formatted <- format.pval(results_df$p_value, digits = 4, eps = 1e-10)
  results_df$p_adj_global_formatted <- format.pval(results_df$p_adj_global, digits = 4, eps = 1e-10)
  
  if ("p_adj_within" %in% colnames(results_df)) {
    results_df$p_adj_within_formatted <- format.pval(results_df$p_adj_within, digits = 4, eps = 1e-10)
  }
  
  return(results_df)
}

#!WILCOXON KRUSKAL
#' Run Wilcoxon/Kruskal-Wallis analysis on data
#'
#' @param data_df Dataframe or list of dataframes with 'Status' column and Number_ columns
#' @param fdr_scope "global", "within", or "both" for FDR correction
#' @param columns_to_test Specific columns to test (default: all Number_ columns)
#' @param run_dunn_posthoc Run Dunn's post-hoc tests for Kruskal-Wallis (TRUE/FALSE)
#'
#' @return Results with metadata
#' @export
run_wilcoxon_kruskal_analysis <- function(data_df, fdr_scope = "both", columns_to_test = NULL, 
                                         run_dunn_posthoc = TRUE) {
  
  # Handle both single dataframe and list
  if (is.list(data_df) && !is.data.frame(data_df)) {
    # Multiple datasets
    results <- run_nonparametric_multiple_datasets(
      data_list = data_df,
      fdr_scope = fdr_scope,
      columns_to_test = columns_to_test,
      run_dunn_posthoc = run_dunn_posthoc
    )
  } else {
    # Single dataset - detect number of groups
    if (!"Status" %in% colnames(data_df)) {
      stop("Dataframe must contain a 'Status' column")
    }
    
    groups <- unique(data_df$Status)
    groups <- groups[!is.na(groups)]
    
    if (length(groups) < 2) {
      stop(sprintf("Need at least 2 groups for statistical testing. Found: %s", 
                   paste(groups, collapse = ", ")))
    }
    
    # Choose test based on number of groups
    if (length(groups) == 2) {
      # Wilcoxon/Mann-Whitney for 2 groups
      results <- run_wilcoxon_single_dataset(
        data_df = data_df,
        fdr_scope = fdr_scope,
        columns_to_test = columns_to_test
      )
    } else {
      # Kruskal-Wallis for >2 groups
      results <- run_kruskal_single_dataset(
        data_df = data_df,
        fdr_scope = fdr_scope,
        columns_to_test = columns_to_test,
        run_dunn_posthoc = run_dunn_posthoc
      )
    }
    
    # Add metadata
    attr(results, "is_single_dataset") <- TRUE
    attr(results, "groups") <- groups
    attr(results, "test_type") <- ifelse(length(groups) == 2, "Wilcoxon", "Kruskal-Wallis")
  }
  
  return(results)
}

#' Run Wilcoxon/Mann-Whitney analysis on single dataset
#'
#' @param data_df Dataframe with 'Status' column and Number_ columns
#' @param fdr_scope FDR correction scope
#' @param columns_to_test Specific columns to test
#'
#' @return Results dataframe
run_wilcoxon_single_dataset <- function(data_df, fdr_scope = "both", columns_to_test = NULL) {
  
  # Get columns to test
  if (is.null(columns_to_test)) {
    columns_to_test <- grep("^Number_", names(data_df), value = TRUE)
  }
  
  if (length(columns_to_test) == 0) {
    stop("No columns found to test")
  }
  
  # Get groups
  groups <- unique(data_df$Status)
  groups <- groups[!is.na(groups)]
  
  if (length(groups) != 2) {
    stop(sprintf("Wilcoxon test requires exactly 2 groups. Found: %s", 
                 paste(groups, collapse = ", ")))
  }
  
  # Run analysis
  results <- analyze_wilcoxon_metrics(data_df, columns_to_test, groups)
  
  # Apply FDR corrections
  results <- apply_wilcoxon_fdr_corrections(results, fdr_scope)
  
  # Add group info
  attr(results, "groups") <- groups
  attr(results, "group1") <- groups[1]
  attr(results, "group2") <- groups[2]
  
  return(results)
}

#' Run Kruskal-Wallis analysis on single dataset
#'
#' @param data_df Dataframe with 'Status' column
#' @param fdr_scope FDR correction scope
#' @param columns_to_test Columns to test
#' @param run_dunn_posthoc Run Dunn's post-hoc tests
#'
#' @return Results list with Kruskal and Dunn results
run_kruskal_single_dataset <- function(data_df, fdr_scope = "both", columns_to_test = NULL, 
                                       run_dunn_posthoc = TRUE) {
  
  # Get columns to test
  if (is.null(columns_to_test)) {
    columns_to_test <- grep("^Number_", names(data_df), value = TRUE)
  }
  
  if (length(columns_to_test) == 0) {
    stop("No columns found to test")
  }
  
  # Get groups
  groups <- unique(data_df$Status)
  groups <- groups[!is.na(groups)]
  
  # Run Kruskal-Wallis tests
  kruskal_results <- analyze_kruskal_metrics(data_df, columns_to_test, groups)
  
  # Apply FDR corrections
  kruskal_results <- apply_kruskal_fdr_corrections(kruskal_results, fdr_scope)
  
  # Run Dunn's post-hoc tests if requested
  dunn_results <- NULL
  if (run_dunn_posthoc) {
    dunn_results <- run_dunn_posthoc_tests(data_df, kruskal_results, "none", groups)

     if (!is.null(dunn_results)) {
      # Apply FDR corrections to Dunn results (consistent with multiple datasets)
      if (fdr_scope %in% c("global", "both")) {
        dunn_results <- dunn_results %>%
          dplyr::mutate(dunn_p_adj_global = p.adjust(dunn_p, method = "fdr"))
      }
      
      if (fdr_scope %in% c("within", "both")) {
        # Within each metric (all comparisons for that metric)
        dunn_results <- dunn_results %>%
          dplyr::group_by(metric) %>%
          dplyr::mutate(dunn_p_adj_within = p.adjust(dunn_p, method = "fdr")) %>%
          dplyr::ungroup()
      }
    }
  }
  
  # Combine results
  results <- list(
    kruskal_results = kruskal_results,
    dunn_results = dunn_results
  )
  
  # Add metadata
  attr(results, "groups") <- groups
  attr(results, "n_groups") <- length(groups)
  class(results) <- "kruskal_results"
  
  return(results)
}

#' Analyze metrics using Wilcoxon/Mann-Whitney test
#'
#' @param df Dataframe with 'Status' column
#' @param columns_to_test Columns to analyze
#' @param groups Vector of 2 group names from Status column
#'
#' @return Results dataframe
analyze_wilcoxon_metrics <- function(df, columns_to_test, groups) {
  
  if (length(groups) != 2) {
    stop("Wilcoxon test requires exactly 2 groups")
  }
  
  group1 <- groups[1]
  group2 <- groups[2]
  
  results <- purrr::map_dfr(columns_to_test, function(col) {
    # Initialize result row
    result_row <- tibble::tibble(
      metric = col,
      group1_zero = FALSE,
      group2_zero = FALSE,
      group1 = group1,
      group2 = group2,
      shapiro_group1_p = NA_real_,
      shapiro_group2_p = NA_real_,
      primary_test = NA_character_,
      primary_p = NA_real_,
      secondary_test = NA_character_,
      secondary_p = NA_real_,
      note = NA_character_
    )
    
    # Split data using actual group names
    group1_data <- df[[col]][df$Status == group1]
    group2_data <- df[[col]][df$Status == group2]
    
    # Check for all-zero cases
    result_row$group1_zero <- all(group1_data == 0)
    result_row$group2_zero <- all(group2_data == 0)
    
    if (result_row$group1_zero || result_row$group2_zero) {
      result_row$note <- dplyr::case_when(
        result_row$group1_zero && result_row$group2_zero ~ sprintf("All zeros in both groups (%s, %s)", group1, group2),
        result_row$group1_zero ~ sprintf("All zeros in %s group", group1),
        TRUE ~ sprintf("All zeros in %s group", group2)
      )
      return(result_row)
    }
    
    # Handle constant values for Shapiro
    safe_shapiro <- function(x) {
      if (length(unique(x)) < 3) return(list(p.value = NA))
      tryCatch(shapiro.test(x), error = function(e) list(p.value = NA))
    }
    
    # Normality checks
    shapiro_group1 <- safe_shapiro(group1_data)
    shapiro_group2 <- safe_shapiro(group2_data)
    result_row$shapiro_group1_p <- shapiro_group1$p.value
    result_row$shapiro_group2_p <- shapiro_group2$p.value
    
    # Determine test strategy
    if (!any(is.na(c(shapiro_group1$p.value, shapiro_group2$p.value))) &&
        shapiro_group1$p.value > 0.05 && 
        shapiro_group2$p.value > 0.05) {
      primary <- "t-test"
      secondary <- "wilcox"
    } else {
      primary <- "wilcox"
      secondary <- "t-test"
    }
    
    # Create formula for testing
    formula_test <- as.formula(paste(col, "~ Status"))
    
    # Perform tests
    tryCatch({
      t_res <- t.test(formula_test, data = df) %>% broom::tidy()
      w_res <- wilcox.test(formula_test, data = df) %>% broom::tidy()
      
      result_row$primary_test <- primary
      result_row$primary_p <- if(primary == "t-test") t_res$p.value else w_res$p.value
      result_row$secondary_test <- secondary
      result_row$secondary_p <- if(secondary == "t-test") t_res$p.value else w_res$p.value
    }, error = function(e) {
      result_row$note <- paste("Test error:", e$message)
    })
    
    result_row
  })
  
  return(results)
}

#' Apply FDR corrections to Wilcoxon results
#'
#' @param results Results dataframe
#' @param fdr_scope FDR correction scope
#'
#' @return Results with FDR corrections
apply_wilcoxon_fdr_corrections <- function(results, fdr_scope = "both") {
  
  if (fdr_scope %in% c("within", "both")) {
    # Within-dataset FDR
    results <- results %>%
      dplyr::mutate(
        within_primary_p_adj = p.adjust(primary_p, method = "fdr"),
        within_secondary_p_adj = p.adjust(secondary_p, method = "fdr")
      )
  }
  
  if (fdr_scope %in% c("global", "both")) {
    # For single dataset, global = within
    if (fdr_scope == "both") {
      results$global_primary_p_adj <- results$within_primary_p_adj
      results$global_secondary_p_adj <- results$within_secondary_p_adj
    } else if (fdr_scope == "global") {
      results$global_primary_p_adj <- p.adjust(results$primary_p, method = "fdr")
      results$global_secondary_p_adj <- p.adjust(results$secondary_p, method = "fdr")
    }
  }
  
  return(results)
}

#' Analyze metrics using Kruskal-Wallis test
#'
#' @param df Dataframe with 'Status' column
#' @param columns_to_test Columns to analyze
#' @param groups Vector of group names from Status column
#'
#' @return Kruskal-Wallis results
analyze_kruskal_metrics <- function(df, columns_to_test, groups) {
  
  # Filter data
  df_filtered <- df %>% 
    dplyr::filter(Status %in% groups)
  
  # Test each metric
  results <- purrr::map_dfr(columns_to_test, function(col) {
    # Initialize result row
    result_row <- tibble::tibble(
      metric = col,
      kruskal_statistic = NA_real_,
      kruskal_p = NA_real_,
      test_valid = TRUE,
      note = NA_character_,
      groups_tested = paste(groups, collapse = ", ")
    )
    
    # Extract data
    test_data <- df_filtered[c("Status", col)]
    colnames(test_data) <- c("group", "value")
    
    # Remove any NA values
    test_data <- test_data %>% dplyr::filter(!is.na(value))
    
    # Check for all-zero cases by group
    zero_groups <- test_data %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(all_zero = all(value == 0), .groups = "drop") %>%
      dplyr::filter(all_zero) %>%
      dplyr::pull(group)
    
    if(length(zero_groups) > 0) {
      result_row$note <- paste("All zeros in groups:", paste(zero_groups, collapse = ", "))
      result_row$test_valid <- FALSE
      return(result_row)
    }
    
    # Check if we have at least 2 groups with data
    groups_with_data <- test_data %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(n = n(), non_zero = any(value != 0), .groups = "drop") %>%
      dplyr::filter(non_zero) %>%
      dplyr::pull(group)
    
    if(length(groups_with_data) < 2) {
      result_row$note <- paste("Insufficient groups with non-zero data. Only found:", 
                               paste(groups_with_data, collapse = ", "))
      result_row$test_valid <- FALSE
      return(result_row)
    }
    
    # Perform Kruskal-Wallis test
    tryCatch({
      kruskal_res <- kruskal.test(value ~ group, data = test_data)
      result_row$kruskal_statistic <- kruskal_res$statistic
      result_row$kruskal_p <- kruskal_res$p.value
      result_row$groups_tested <- paste(groups_with_data, collapse = ", ")
    }, error = function(e) {
      result_row$note <- paste("Kruskal-Wallis test error:", e$message)
      result_row$test_valid <- FALSE
    })
    
    result_row
  })
  
  return(results)
}

#' Apply FDR corrections to Kruskal results
#'
#' @param results Kruskal results dataframe
#' @param fdr_scope FDR correction scope
#'
#' @return Results with FDR corrections
apply_kruskal_fdr_corrections <- function(results, fdr_scope = "both") {
  
  if (fdr_scope %in% c("within", "both")) {
    # Within-dataset FDR
    results <- results %>%
      dplyr::mutate(within_kruskal_p_adj = p.adjust(kruskal_p, method = "fdr"))
  }
  
  if (fdr_scope %in% c("global", "both")) {
    # For single dataset, global = within
    if (fdr_scope == "both") {
      results$global_kruskal_p_adj <- results$within_kruskal_p_adj
    } else if (fdr_scope == "global") {
      results$global_kruskal_p_adj <- p.adjust(results$kruskal_p, method = "fdr")
    }
  }
  
  return(results)
}

#' Run Dunn's post-hoc tests for significant Kruskal-Wallis results
#'
#' @param data_df Original dataframe
#' @param kruskal_results Kruskal-Wallis results
#' @param fdr_scope FDR correction scope
#' @param groups Vector of group names
#' @param significance_threshold Threshold for significance (default: 0.05)
#'
#' @return Dunn's test results
run_dunn_posthoc_tests <- function(data_df, kruskal_results, fdr_scope = "none", 
                                   groups = NULL, significance_threshold = 0.05) {
  
  # If groups not provided, get them from data
  if (is.null(groups)) {
    groups <- unique(data_df$Status)
    groups <- groups[!is.na(groups)]
  }
  
  # Filter significant results
  if (fdr_scope %in% c("within", "both") && "within_kruskal_p_adj" %in% names(kruskal_results)) {
    significant_tests <- kruskal_results %>%
      dplyr::filter(within_kruskal_p_adj < significance_threshold & test_valid)
  } else if (fdr_scope == "global" && "global_kruskal_p_adj" %in% names(kruskal_results)) {
    significant_tests <- kruskal_results %>%
      dplyr::filter(global_kruskal_p_adj < significance_threshold & test_valid)
  } else {
    # No FDR applied yet
    significant_tests <- kruskal_results %>%
      dplyr::filter(kruskal_p < significance_threshold & test_valid)
  }
  
  if (nrow(significant_tests) == 0) {
    return(NULL)
  }
  
  dunn_results_list <- list()
  
  for (i in 1:nrow(significant_tests)) {
    met <- significant_tests$metric[i]
    
    # Run Dunn's test
    dunn_res <- run_dunn_test(data_df, met, groups)
    
    if (!is.null(dunn_res) && nrow(dunn_res) > 0) {
      dunn_res <- dunn_res %>%
        dplyr::mutate(metric = met)
      dunn_results_list[[length(dunn_results_list) + 1]] <- dunn_res
    }
  }
  
  if (length(dunn_results_list) == 0) {
    return(NULL)
  }
  
  combined_dunn <- dplyr::bind_rows(dunn_results_list)
  
  return(combined_dunn)
}

#' Run Dunn's test
#'
#' @param df Dataframe with 'Status' column
#' @param metric Column to test
#' @param groups Vector of group names
#'
#' @return Dunn's test results
run_dunn_test <- function(df, metric, groups) {
  
  df_filtered <- df %>% 
    dplyr::filter(Status %in% groups)
  
  test_data <- df_filtered[c("Status", metric)]
  colnames(test_data) <- c("group", "value")
  
  # Remove NA values
  test_data <- test_data %>% 
    dplyr::filter(!is.na(value))
  
  # Check if we have enough groups
  unique_groups <- unique(test_data$group)
  if (length(unique_groups) < 2) {
    return(NULL)
  }
  
  tryCatch({
    # Perform Dunn's test
    dunn_res <- dunn.test::dunn.test(
      x = test_data$value,
      g = test_data$group,
      method = "none",  # Get raw p-values
      list = TRUE
    )
    
    # Return results as a tibble
    tibble::tibble(
      comparison = dunn_res$comparisons,
      dunn_z = dunn_res$Z,
      dunn_p = dunn_res$P,
      groups_in_test = paste(unique_groups, collapse = ", ")
    )
  }, error = function(e) {
    return(NULL)
  })
}

#' Run non-parametric analysis on multiple datasets
#'
#' @param data_list List of dataframes
#' @param fdr_scope FDR correction scope
#' @param columns_to_test Columns to test
#' @param run_dunn_posthoc Run Dunn's post-hoc tests
#'
#' @return Combined results
run_nonparametric_multiple_datasets <- function(data_list, fdr_scope = "both", 
                                                columns_to_test = NULL, 
                                                run_dunn_posthoc = TRUE) {
  
  if (is.null(names(data_list))) {
    names(data_list) <- paste0("Dataset_", seq_along(data_list))
  }
  
  all_results <- list()
  all_groups_info <- list()
  all_test_types <- list()
  
  for (dataset_name in names(data_list)) {
    cat(sprintf("Processing dataset: %s\n", dataset_name))
    
    tryCatch({
      # Get groups for this dataset
      groups <- unique(data_list[[dataset_name]]$Status)
      groups <- groups[!is.na(groups)]
      
      if (length(groups) < 2) {
        warning(sprintf("Dataset %s: Need at least 2 groups. Found: %s", 
                       dataset_name, paste(groups, collapse = ", ")))
        next
      }
      
      # Run appropriate test
      if (length(groups) == 2) {
        # Wilcoxon
        results <- run_wilcoxon_single_dataset(
          data_df = data_list[[dataset_name]],
          fdr_scope = "none",  # We'll apply FDR later
          columns_to_test = columns_to_test
        )
        
        # Add dataset name
        results$Dataset <- dataset_name
        
        all_results[[dataset_name]] <- results
        all_test_types[[dataset_name]] <- "Wilcoxon"
        
      } else {
        # Kruskal-Wallis
        results <- run_kruskal_single_dataset(
          data_df = data_list[[dataset_name]],
          fdr_scope = "none",  # We'll apply FDR later
          columns_to_test = columns_to_test,
          run_dunn_posthoc = FALSE  # We'll handle Dunn tests separately
        )
        
        # Add dataset name to kruskal results
        kruskal_df <- results$kruskal_results
        kruskal_df$Dataset <- dataset_name
        
        # Store results
        results_list <- list(
          kruskal_results = kruskal_df,
          dunn_results = NULL
        )
        
        # Run Dunn's tests if requested
        if (run_dunn_posthoc) {
          dunn_df <- run_dunn_posthoc_tests(
            data_df = data_list[[dataset_name]],
            kruskal_results = kruskal_df,
            fdr_scope = "none",
            groups = groups
          )
          
          if (!is.null(dunn_df)) {
            dunn_df$Dataset <- dataset_name
            results_list$dunn_results <- dunn_df
          }
        }
        
        all_results[[dataset_name]] <- results_list
        all_test_types[[dataset_name]] <- "Kruskal-Wallis"
      }
      
      all_groups_info[[dataset_name]] <- groups
      
    }, error = function(e) {
      warning(sprintf("Error processing dataset %s: %s", dataset_name, e$message))
    })
  }
  
  if (length(all_results) == 0) {
    stop("No datasets were successfully processed")
  }
  
  # Separate Wilcoxon and Kruskal results
  wilcoxon_dfs <- list()
  kruskal_dfs <- list()
  dunn_dfs <- list()
  
  for (dataset_name in names(all_results)) {
    if (all_test_types[[dataset_name]] == "Wilcoxon") {
      wilcoxon_dfs[[dataset_name]] <- all_results[[dataset_name]]
    } else {
      kruskal_dfs[[dataset_name]] <- all_results[[dataset_name]]$kruskal_results
      if (!is.null(all_results[[dataset_name]]$dunn_results)) {
        dunn_dfs[[dataset_name]] <- all_results[[dataset_name]]$dunn_results
      }
    }
  }
  
  # Combine results
  final_results <- list()
  
  # Combine Wilcoxon results
  if (length(wilcoxon_dfs) > 0) {
    combined_wilcoxon <- dplyr::bind_rows(wilcoxon_dfs)
    
    # Apply global FDR corrections
    if (fdr_scope %in% c("global", "both")) {
      combined_wilcoxon <- combined_wilcoxon %>%
        dplyr::mutate(
          global_primary_p_adj = p.adjust(primary_p, method = "fdr"),
          global_secondary_p_adj = p.adjust(secondary_p, method = "fdr")
        )
    }
    
    # Apply within-dataset FDR corrections
    if (fdr_scope %in% c("within", "both")) {
      combined_wilcoxon <- combined_wilcoxon %>%
        dplyr::group_by(Dataset) %>%
        dplyr::mutate(
          within_primary_p_adj = p.adjust(primary_p, method = "fdr"),
          within_secondary_p_adj = p.adjust(secondary_p, method = "fdr")
        ) %>%
        dplyr::ungroup()
    }
    
    final_results$wilcoxon_results <- combined_wilcoxon
  }
  
  # Combine Kruskal results
  if (length(kruskal_dfs) > 0) {
    combined_kruskal <- dplyr::bind_rows(kruskal_dfs)
    
    # Apply global FDR corrections
    if (fdr_scope %in% c("global", "both")) {
      combined_kruskal <- combined_kruskal %>%
        dplyr::mutate(global_kruskal_p_adj = p.adjust(kruskal_p, method = "fdr"))
    }
    
    # Apply within-dataset FDR corrections
    if (fdr_scope %in% c("within", "both")) {
      combined_kruskal <- combined_kruskal %>%
        dplyr::group_by(Dataset) %>%
        dplyr::mutate(within_kruskal_p_adj = p.adjust(kruskal_p, method = "fdr")) %>%
        dplyr::ungroup()
    }
    
    final_results$kruskal_results <- combined_kruskal
  }
  
  # Combine Dunn results if available
  if (length(dunn_dfs) > 0) {
    combined_dunn <- dplyr::bind_rows(dunn_dfs)
    
    # Apply FDR corrections to Dunn results
    if (fdr_scope %in% c("global", "both")) {
      combined_dunn <- combined_dunn %>%
        dplyr::mutate(dunn_p_adj_global = p.adjust(dunn_p, method = "fdr"))
    }
    
    if (fdr_scope %in% c("within", "both")) {
      # Within each dataset-metric combination
      combined_dunn <- combined_dunn %>%
        dplyr::group_by(Dataset, metric) %>%
        dplyr::mutate(dunn_p_adj_within = p.adjust(dunn_p, method = "fdr")) %>%
        dplyr::ungroup()
    }
    
    final_results$dunn_results <- combined_dunn
  }
  
  # Add metadata
  attr(final_results, "datasets") <- names(all_results)
  attr(final_results, "n_datasets") <- length(all_results)
  attr(final_results, "test_types") <- all_test_types
  attr(final_results, "groups_info") <- all_groups_info
  attr(final_results, "fdr_scope") <- fdr_scope
  attr(final_results, "is_single_dataset") <- FALSE
  
  class(final_results) <- "nonparametric_results"
  
  return(final_results)
}


#' Convert results to dataframe format
#' 
#' @param results Output from run_wilcoxon_kruskal_analysis
#' @return Dataframe with all results
#' @export
results_to_dataframe <- function(results) {
  
  # If it's already a dataframe (Wilcoxon), return as is
  if (is.data.frame(results)) {
    return(results)
  }
  
  # If it's a list with kruskal_results (single dataset)
  if (is.list(results) && "kruskal_results" %in% names(results)) {
    kruskal_df <- results$kruskal_results
    
    # Add Dunn results if available
    if (!is.null(results$dunn_results)) {
      dunn_df <- results$dunn_results
      
      # Add indicator column
      kruskal_df$test_type <- "Kruskal-Wallis"
      dunn_df$test_type <- "Dunn_posthoc"
      
      # Combine
      combined <- dplyr::bind_rows(
        list(kruskal = kruskal_df, dunn = dunn_df),
        .id = "analysis_level"
      )
      
      return(combined)
    } else {
      # Just Kruskal results
      kruskal_df$test_type <- "Kruskal-Wallis"
      kruskal_df$analysis_level <- "kruskal"
      return(kruskal_df)
    }
  }
  
  # If it's the multiple datasets list structure
  if (is.list(results) && ("wilcoxon_results" %in% names(results) || 
                           "kruskal_results" %in% names(results))) {
    
    result_parts <- list()
    
    # Add Wilcoxon results if present
    if (!is.null(results$wilcoxon_results)) {
      wilcoxon_df <- results$wilcoxon_results
      wilcoxon_df$test_type <- "Wilcoxon"
      wilcoxon_df$analysis_level <- "primary"
      result_parts$wilcoxon <- wilcoxon_df
    }
    
    # Add Kruskal results if present
    if (!is.null(results$kruskal_results)) {
      kruskal_df <- results$kruskal_results
      kruskal_df$test_type <- "Kruskal-Wallis"
      kruskal_df$analysis_level <- "primary"
      result_parts$kruskal <- kruskal_df
    }
    
    # Add Dunn results if present
    if (!is.null(results$dunn_results)) {
      dunn_df <- results$dunn_results
      dunn_df$test_type <- "Dunn_posthoc"
      dunn_df$analysis_level <- "posthoc"
      result_parts$dunn <- dunn_df
    }
    
    # Combine everything
    if (length(result_parts) > 0) {
      return(dplyr::bind_rows(result_parts, .id = "dataset_type"))
    }
  }
  
  # If we get here, return the original results
  warning("Could not convert results to dataframe. Returning original format.")
  return(results)
}

# #' Print method for kruskal_results
# #'
# #' @param x kruskal_results object
# #' @param ... Additional arguments
# #'
# #' @export
# print.kruskal_results <- function(x, ...) {
#   cat("Kruskal-Wallis Analysis Results\n")
#   cat("===============================\n\n")
  
#   if (is.list(x) && "kruskal_results" %in% names(x)) {
#     cat("Kruskal-Wallis tests:", nrow(x$kruskal_results), "\n")
#     if (!is.null(x$dunn_results)) {
#       cat("Dunn's post-hoc tests:", nrow(x$dunn_results), "\n")
#     }
#   } else {
#     NextMethod()
#   }
# }

#' Print method for nonparametric_results
#'
#' @param x nonparametric_results object
#' @param ... Additional arguments
#'
#' @export
print.nonparametric_results <- function(x, ...) {
  cat("Non-parametric Analysis Results\n")
  cat("===============================\n\n")
  
  if (!is.null(x$wilcoxon_results)) {
    cat("Wilcoxon/Mann-Whitney tests:", nrow(x$wilcoxon_results), "\n")
  }
  if (!is.null(x$kruskal_results)) {
    cat("Kruskal-Wallis tests:", nrow(x$kruskal_results), "\n")
  }
  if (!is.null(x$dunn_results)) {
    cat("Dunn's post-hoc tests:", nrow(x$dunn_results), "\n")
  }
  
  cat("\nDatasets processed:", attr(x, "n_datasets", exact = TRUE), "\n")
  cat("FDR scope:", attr(x, "fdr_scope", exact = TRUE), "\n")
}

#!' Calculate rates for count data
calculate_rates <- function(data_df, manifest, rows = "all", clean_names = TRUE, 
                           row_selection = "variants", private = FALSE, 
                           name = NULL, gene_lists = NULL, ...) {
  
  # Extract manifest info for group totals
  manifest_info <- extract_groups_from_manifest(manifest)
  group_totals <- manifest_info$counts
  
  # Check if it's a single dataframe or a list of dataframes
  if (is.data.frame(data_df)) {
    # Single dataframe
    cat(sprintf("Processing single dataframe with row_selection = '%s'...\n", row_selection))
    
    if (is.null(name)) {
      stop("For raw variant data, you must provide a 'name' parameter")
    }
    if (is.null(gene_lists)) {
      stop("For raw variant data, you must provide 'gene_lists' parameter")
    }
    
    # Create count dataframe first
    cat(sprintf("Creating count dataframe with name='%s'...\n", name))
    count_df <- create_count_dataframes(
      df = data_df,
      name = name,
      gene_lists = gene_lists,
      manifest_df = manifest,
      row_selection = row_selection,
      private = private,
      ...
    )
    
    # Calculate rates
    results <- calculate_rates_from_counts(
      count_df = count_df,
      group_totals = group_totals,
      rows = rows,
      clean_names = clean_names
    )
    
    # Add metadata
    attr(results, "is_single_dataset") <- TRUE
    attr(results, "row_selection") <- row_selection
    attr(results, "dataset_name") <- name
    
    return(results)
    
  } else if (is.list(data_df) && !is.data.frame(data_df)) {
    # List of dataframes
    cat(sprintf("Processing list of %d dataframes with row_selection = '%s'...\n", 
                length(data_df), row_selection))
    
    all_results <- list()
    
    # Process each dataframe in the list
    for (i in seq_along(data_df)) {
      dataset_name <- names(data_df)[i]
      if (is.null(dataset_name) || dataset_name == "") {
        dataset_name <- paste0("Dataset_", i)
      }
      
      cat(sprintf("  Processing dataset %d: %s\n", i, dataset_name))
      
      tryCatch({
        # Get the dataframe
        df_to_process <- data_df[[i]]
        
        # Check if it's a list with parameters
        if (is.list(df_to_process) && !is.data.frame(df_to_process)) {
          # This is a list containing parameters
          if ("data" %in% names(df_to_process)) {
            df_data <- df_to_process$data
          } else {
            # Assume first element is dataframe
            df_data <- df_to_process[[1]]
          }
          
          # Get name from parameters or use default
          dataset_name_param <- df_to_process$name
          if (!is.null(dataset_name_param)) {
            dataset_name <- dataset_name_param
          }
          
          # Get gene_lists from parameters or use global
          dataset_gene_lists <- df_to_process$gene_lists
          if (is.null(dataset_gene_lists)) {
            dataset_gene_lists <- gene_lists
          }
        } else {
          # It's just a dataframe
          df_data <- df_to_process
          dataset_gene_lists <- gene_lists
        }
        
        # Check if dataframe is already processed
        if ("row_names" %in% names(df_data) || 
            (ncol(df_data) > 0 && is.character(df_data[[1]]) && 
             grepl("Number of", df_data[1, 1]))) {
          # Already processed count dataframe
          cat("    Data appears to be already processed count dataframe\n")
          count_df <- df_data
        } else {
          # Raw variant data - need to process it
          if (is.null(dataset_gene_lists)) {
            stop(sprintf("For dataset '%s', you must provide 'gene_lists'", dataset_name))
          }
          
          # Create count dataframe
          count_df <- create_count_dataframes(
            df = df_data,
            name = dataset_name,  # Use dataset name
            gene_lists = dataset_gene_lists,
            manifest_df = manifest,
            row_selection = row_selection,
            private = private,
            ...
          )
        }
        
        # Calculate rates for this dataset
        single_result <- calculate_rates_from_counts(
          count_df = count_df,
          group_totals = group_totals,
          rows = rows,
          clean_names = clean_names
        )
        
        if (!is.null(single_result)) {
          single_result$Dataset <- dataset_name
          all_results[[dataset_name]] <- single_result
        }
        
      }, error = function(e) {
        warning(sprintf("Error processing dataset %s: %s", dataset_name, e$message))
      })
    }
    
    if (length(all_results) == 0) {
      stop("No datasets were successfully processed")
    }
    
    # Combine all results
    combined_results <- do.call(rbind, all_results)
    
    # Add metadata
    attr(combined_results, "datasets") <- names(all_results)
    attr(combined_results, "n_datasets") <- length(all_results)
    attr(combined_results, "row_selection") <- row_selection
    attr(combined_results, "is_single_dataset") <- FALSE
    
    return(combined_results)
    
  } else {
    stop("data_df must be either a dataframe or a list of dataframes")
  }
}


#' Calculate rates from already processed count dataframe
#'
#' @param count_df Processed count dataframe
#' @param group_totals List of group totals
#' @param rows Which rows to calculate rates for
#' @param clean_names Clean row names
#'
#' @return Rate calculation results
calculate_rates_from_counts <- function(count_df, group_totals, rows = "all", clean_names = TRUE) {
  
  # Validate group totals
  if (!is.list(group_totals) || length(group_totals) == 0) {
    stop("group_totals must be a non-empty list")
  }
  
  # Check that all groups in count_df exist in group_totals
  data_groups <- names(count_df)[-1]  # Assuming first column is row names
  missing_in_totals <- setdiff(data_groups, names(group_totals))
  if (length(missing_in_totals) > 0) {
    stop(sprintf("Groups in data_df not found in group_totals: %s\nAvailable in group_totals: %s",
                 paste(missing_in_totals, collapse = ", "),
                 paste(names(group_totals), collapse = ", ")))
  }
  
  # Select rows
  if (identical(rows, "all")) {
    selected_rows <- 1:nrow(count_df)
  } else {
    selected_rows <- rows
    # Validate row indices
    invalid_rows <- setdiff(selected_rows, 1:nrow(count_df))
    if (length(invalid_rows) > 0) {
      stop(sprintf("Invalid row indices: %s. Valid range: 1-%d",
                   paste(invalid_rows, collapse = ", "), nrow(count_df)))
    }
  }
  
  # Calculate rates for selected rows
  results_list <- list()
  
  for (row_num in selected_rows) {
    row_data <- count_df[row_num, ]
    row_name <- row_data[[1]]  # Get row name from first column
    
    # Calculate rates for all groups
    rate_values <- numeric()
    count_values <- numeric()
    
    for (group in names(group_totals)) {
      count <- row_data[[group]]
      total <- group_totals[[group]]
      
      if (is.na(count) || is.na(total)) {
        warning(sprintf("Missing count or total for group '%s' in row %d", group, row_num))
        rate <- NA
      } else if (total == 0) {
        warning(sprintf("Total is 0 for group '%s' in row %d", group, row_num))
        rate <- NA
      } else {
        rate <- count / total
      }
      
      rate_values[group] <- rate
      count_values[group] <- count
    }
    
    # Create result row
    result_row <- data.frame(
      Row_Name = row_name,
      stringsAsFactors = FALSE
    )
    
    # Add rate, count, and total columns for each group
    for (group in names(group_totals)) {
      result_row[[paste0(group, "_Rate")]] <- round(rate_values[group], 6)
      result_row[[paste0(group, "_Count")]] <- count_values[group]
      result_row[[paste0(group, "_Total")]] <- group_totals[[group]]
    }
    
    results_list[[length(results_list) + 1]] <- result_row
  }
  
  # Combine all results
  if (length(results_list) == 0) {
    return(NULL)
  }
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  
  # Clean names if requested
  if (clean_names && !is.null(results) && nrow(results) > 0) {
    results <- clean_rate_row_names(results)
  }
  
  return(results)
}

#' Clean row names in rate results
#'
#' @param results Rate results dataframe
#'
#' @return Cleaned dataframe
clean_rate_row_names <- function(results) {
  if (is.null(results) || nrow(results) == 0) {
    return(results)
  }
  
  # Create a clean version
  results$Row_Name_Clean <- results$Row_Name
  
  # Remove everything before and including "with"
  results$Row_Name_Clean <- gsub(".*with ", "", results$Row_Name_Clean)
  
  # Remove the word before "variants" but keep "variants"
  results$Row_Name_Clean <- gsub("\\S+ variants", "variants", results$Row_Name_Clean)
  
  # Final cleanup
  results$Row_Name_Clean <- trimws(results$Row_Name_Clean)

  results$Row_Name <- results$Row_Name_Clean 
  results$Row_Name_Clean <- NULL
  return(results)
}

