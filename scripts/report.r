
#for genelists 
#' Create Detailed Gene List Summary
#'
#' Creates a comprehensive summary table showing which genes are in which gene lists,
#' their association with TRs, outliers, and other metadata.
#'
#' @param df The input dataframe containing TR information
#' @param gene_lists A named list of gene vectors for analysis
#' @param group_column Column name for outlier groups (default: "outlier_label")
#' @param gene_column Column name for gene names (default: "Gene")
#' @param sample_column Column name for sample IDs/outliers (default: "outliers")
#' @param region_column Column name for regions (default: "region")
#' @param motif_column Column name for motifs (default: "motif")
#' @param tr_id_column Column name for TR identifiers (default: "repeatID")
#' 
#' @return A list with detailed and wide format summaries
#' 
#' @import dplyr
#' @import tidyr
#'
create_gene_list_summary <- function(df,
                                     gene_lists,
                                     output_dir = getwd(),
                                     group_column = "sample_label",
                                     gene_column = "SYMBOL",
                                     sample_column = "SAMPLES"
                      ) {
  
  # Validate input
  required_columns <- c(group_column, gene_column, sample_column)
  missing_cols <- setdiff(required_columns, names(df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure gene lists have unique names
  if (is.null(names(gene_lists))) {
    names(gene_lists) <- paste0("list_", seq_along(gene_lists))
  }
  
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # # Function to split and count samples correctly
  # split_and_count_samples <- function(sample_string) {
  #   if (is.na(sample_string) || sample_string == "") {
  #     return(list(samples = character(0), count = 0))
  #   }
    
  #   # Split by common delimiters: comma, semicolon, or space
  #   samples <- unlist(strsplit(sample_string, split = "[,; ]+"))
    
  #   # Remove empty strings
  #   samples <- samples[samples != ""]
    
  #   # Get unique samples
  #   unique_samples <- unique(samples)
    
  #   return(list(samples = unique_samples, count = length(unique_samples)))
  # }

  split_and_count_samples <- function(sample_string) {

  if (is.na(sample_string) || sample_string == "") {
    return(list(samples = character(0), count = 0))
  }

  samples <- unlist(strsplit(sample_string,  split = "[,; ]+"))
  samples <- trimws(samples)

  list(
    samples = samples,
    count = length(samples)
  )
}

  
  # --- Main Processing -----------------------------------------------------
  
  message("\n===Creating detailed gene list summary...===\n")
  message(sprintf("Total rows in dataset: %d", nrow(df)))
  message(sprintf("Total unique genes: %d", n_distinct(df[[gene_column]])))
  message(sprintf("Gene lists to analyze: %s", paste(names(gene_lists), collapse = ", ")))
  
  # Get all unique groups
  groups <- unique(df[[group_column]])
  groups <- groups[!is.na(groups)]
  message(sprintf("Outlier groups found: %s", paste(groups, collapse = ", ")))
  
  all_details <- list()
  detail_counter <- 1
  
  df$Variant_ID <- paste(df$CHROM, df$POS, df$REF, df$ALT, sep = "_")

  # Process each group
  for (group in groups) {
    message(sprintf("  Processing group: %s", group))
    
    # Subset dataframe for this group
    df_group <- df %>% filter(!!sym(group_column) == group)
    
    if (nrow(df_group) == 0) next
    
    # Process each gene list for this group
    for (list_name in names(gene_lists)) {
      message(sprintf("    Analyzing gene list: %s", list_name))
      
      gene_list <- gene_lists[[list_name]]
      
      # Skip empty gene lists
      if (length(gene_list) == 0) next
      
      # Filter dataframe for genes in this list
      df_filtered <- df_group %>% 
        filter(!!sym(gene_column) %in% gene_list)
      
      if (nrow(df_filtered) == 0) next
      
      # Add gene list and process samples for each row
      df_with_list <- df_filtered %>%
        mutate(
          gene_list = list_name,
          sample_info = purrr::map(.data[[sample_column]], split_and_count_samples),
          sample_count = purrr::map_int(sample_info, "count"),
          all_samples = purrr::map_chr(sample_info, ~ paste(.x$samples, collapse = ", "))
        ) %>%
        select(-sample_info)

      
      # Add to details list
      all_details[[detail_counter]] <- df_with_list
      detail_counter <- detail_counter + 1
    }
  }
  
  # Combine all details
  if (length(all_details) == 0) {
    message("No genes found in any gene list/group combination.")
    return(list(
      detailed_summary = data.frame(),
      wide_summary = data.frame()
    ))
  }
  
  detailed_summary <- bind_rows(all_details) %>%
    select(
      gene_list,
      outlier_label = !!sym(group_column),
      gene = !!sym(gene_column),
      # tr_id = !!sym(tr_id_column),
      # region = !!sym(region_column),
      # motif = !!sym(motif_column),
      Variant_ID,
      sample_count,
      all_samples
    ) %>%
    arrange(gene_list, outlier_label, gene, Variant_ID)
  
  # --- Create Wide Format Summary (aggregated) -----------------------------
  
  message("\nCreating wide format summary...")
  
  # Create aggregated wide format summary
  wide_summary <- detailed_summary %>%
    group_by(gene_list, outlier_label) %>%
    summarise(
      # All unique genes
      all_genes = paste(sort(unique(gene)), collapse = ", "),
      num_genes = n_distinct(gene),
      
      # All unique Variant IDs
      all_Variant_IDS = paste(sort(unique(Variant_ID)), collapse = ", "),
      num_variants = n_distinct(Variant_ID),
      
      # All unique samples
      all_sample_ids = paste(unique(unlist(strsplit(all_samples, ","))), collapse = ","),
      num_samples = sum(sample_count),
      
      # # All unique regions
      # all_region_types = paste(sort(unique(region)), collapse = "; "),
      
      # # All unique motifs
      # all_motif_types = paste(sort(unique(motif)), collapse = "; "),
      
      .groups = 'drop'
    ) %>%
    arrange(gene_list, outlier_label)
  
  # --- Create Gene-List-Group Details (one row per gene) -------------------
  
  message("\nCreating gene-level summary...")
  
  gene_level_summary <- detailed_summary %>%
    group_by(gene_list, outlier_label, gene) %>%
    summarise(
      tr_count = n_distinct(Variant_ID),
      sample_count = sum(sample_count),
      all_Variant_IDS = paste(sort(unique(Variant_ID)), collapse = ", "),
      all_samples = paste(unique(unlist(strsplit(all_samples, ", "))), collapse = ", "),
      # all_regions = paste(sort(unique(region)), collapse = "; "),
      # all_motifs = paste(sort(unique(motif)), collapse = "; "),
      .groups = 'drop'
    ) %>%
    arrange(gene_list, outlier_label, desc(tr_count), gene)
  
  # --- Output Results ------------------------------------------------------
  
  message("\n=== SUMMARY STATISTICS ===")
  message(sprintf("Total TR rows in detailed summary: %d", nrow(detailed_summary)))
  message(sprintf("Total unique genes: %d", n_distinct(detailed_summary$gene)))
  message(sprintf("Total unique TRs: %d", n_distinct(detailed_summary$Variant_ID)))
  
  # Print a quick overview
  for (list_name in names(gene_lists)) {
    list_data <- detailed_summary %>% filter(gene_list == list_name)
    if (nrow(list_data) > 0) {
      message(sprintf("\nGene list: %s", list_name))
      for (group in groups) {
        group_data <- list_data %>% filter(outlier_label == group)
        if (nrow(group_data) > 0) {
          message(sprintf("  Group %s: %d genes, %d TRs, %d sample associations", 
                         group, n_distinct(group_data$gene),
                         n_distinct(group_data$Variant_ID), sum(group_data$sample_count)))
        }
      }
    }
  }
  
  # Write to Excel with three sheets
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    
    wb <- createWorkbook()
    
    # 1. Detailed Summary sheet (one row per TR)
    addWorksheet(wb, "Detailed Summary (per TR)")
    writeData(wb, "Detailed Summary (per TR)", detailed_summary)
    setColWidths(wb, "Detailed Summary (per TR)", cols = 1:ncol(detailed_summary), 
                widths = "auto")
    
    # 2. Gene-Level Summary (one row per gene)
    addWorksheet(wb, "Gene-List-Group Details")
    writeData(wb, "Gene-List-Group Details", gene_level_summary)
    setColWidths(wb, "Gene-List-Group Details", cols = 1:ncol(gene_level_summary), 
                 widths = "auto")
    
    # 3. Wide Format Summary sheet (aggregated by gene list and group)
    addWorksheet(wb, "Wide Format Summary")
    writeData(wb, "Wide Format Summary", wide_summary)
    setColWidths(wb, "Wide Format Summary", cols = 1:ncol(wide_summary), 
                 widths = "auto")
    
    # Save the workbook
    output_file <- file.path(output_dir, "gene_list_summary.xlsx")
    saveWorkbook(wb, output_file, overwrite = TRUE)
    message(sprintf("\nSummary saved to: %s", output_file))
    message("Sheets included:")
    message("  1. Detailed Summary (per TR) - One row per TR")
    message("  2. Gene-List-Group Details - One row per gene (aggregated)")
    message("  3. Wide Format Summary - Aggregated by gene list and group")
    message("\n")
  }
  
  # Return a list with all summaries
  return(list(
    detailed_summary = detailed_summary,      # One row per TR
    gene_level_summary = gene_level_summary,  # One row per gene
    wide_summary = wide_summary               # Aggregated by gene list and group
  ))
}

create_gene_lists_excel <- function(df,
                                    gene_lists,
                                    output_dir = getwd(),
                                    dataset_name = "dataset",
                                    output_prefix = "group_datasets_",
                                    group_column = "sample_label",
                                    gene_column = "SYMBOL",
                                    dedup_columns = c("SYMBOL", "CHROM", "REF", "ALT", "POS", "SAMPLES"),
                                    selected_columns_gene_list = c("CHROM", "POS", "REF", "ALT", "Allele", "IMPACT", "SYMBOL", 
                                "Consequence", "am_pathogenicity", "LoF", "CANONICAL", "Gene",
                                "BIOTYPE", "Feature", "SAMPLES", "GENOTYPES",
                                "lof.oe_ci.upper", "cds_length", "num_coding_exons",
                                "sample_label_2", "sample_label", "count",
                                "type", "group"),
                                    selected_columns = c("CHROM", "POS", "REF", "ALT", "Allele", "IMPACT", "SYMBOL",
                      "Consequence", "am_pathogenicity", "LoF", "CANONICAL", "Gene",
                      "BIOTYPE", "Feature", "SAMPLES", "GENOTYPES", 
                      "lof.oe_ci.upper", "cds_length", "num_coding_exons",
                      "sample_label_2", "sample_label",  "count")) {
  
  # --- Handle both single dataframe and list of dataframes ---
  if (is.list(df) && !is.data.frame(df)) {
    # It's a list of dataframes - process each one
    message(sprintf("\n==== Processing list of %d datasets ===\n", length(df)))
    
    all_results <- list()
    
    for (i in seq_along(df)) {
      current_name <- names(df)[i]
      if (is.null(current_name) || current_name == "") {
        current_name <- paste0("Dataset_", i)
      }
      
      message(sprintf("\nProcessing dataset %d: %s", i, current_name))
      
      tryCatch({
        # Get the dataframe
        df_to_process <- df[[i]]
        
        # If it's a list with parameters, extract the dataframe
        if (is.list(df_to_process) && !is.data.frame(df_to_process)) {
          if ("data" %in% names(df_to_process)) {
            df_data <- df_to_process$data
            # Use name from parameters if available
            if (!is.null(df_to_process$name)) {
              current_name <- df_to_process$name
            }
          } else {
            # Assume first element is dataframe
            df_data <- df_to_process[[1]]
          }
        } else {
          # It's just a dataframe
          df_data <- df_to_process
        }
        
        # Process single dataframe
        result <- create_gene_lists_excel(
          df = df_data,
          dataset_name = current_name,
          gene_lists = gene_lists,
          output_dir = output_dir,
          output_prefix = output_prefix,
          group_column = group_column,
          gene_column = gene_column,
          dedup_columns = dedup_columns,
          selected_columns_gene_list = selected_columns_gene_list,
          selected_columns = selected_columns
        )
        
        all_results[[current_name]] <- result
        
      }, error = function(e) {
        warning(sprintf("Error processing dataset %s: %s", current_name, e$message))
      })
    }
    
    return(invisible(all_results))
    
  } else if (!is.data.frame(df)) {
    # Not a dataframe or list
    stop("Input 'df' must be a dataframe or a list of dataframes")
  }

  # --- Single dataframe processing (original logic) ---
  # Validate input
  if (!group_column %in% names(df)) {
    stop(sprintf("Group column '%s' not found in dataframe", group_column))
  }
  
  if (!gene_column %in% names(df)) {
    stop(sprintf("Gene column '%s' not found in dataframe", gene_column))
  }
  
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }


  # --- Helper Functions ----------------------------------------------------
  
  # Split dataframe by group column
  split_by_group <- function(df, group_col = group_column) {
    # Remove any NA values in the group column
    df <- df %>% filter(!is.na(!!sym(group_col)))
    
    # Split by group
    groups <- unique(df[[group_col]])
    result <- list()
    
    for (group in groups) {
      group_df <- df %>% filter(!!sym(group_col) == group)
      # Deduplicate
      group_df <- distinct(group_df, across(all_of(dedup_columns)), .keep_all = TRUE)
      result[[as.character(group)]] <- group_df
    }
    
    return(result)
  }
  
  # # Save unique genes for each group
  # save_unique_genes_all <- function(split_data, prefix = output_prefix) {
  #   # Set working directory to output_dir
  #   old_wd <- getwd()
  #   setwd(output_dir)
  #   on.exit(setwd(old_wd))  # Reset on exit
    
  #   for (group in names(split_data)) {
  #     df_group <- split_data[[group]]
  #     if (nrow(df_group) > 0) {
  #       unique_df <- df_group %>% 
  #         distinct(across(all_of(gene_column)), .keep_all = TRUE) %>% 
  #         select(any_of(intersect(selected_columns, names(df_group))))
        
  #       file_name <- paste0(prefix, tolower(group), "_", dataset_name, ".tsv")
  #       write.table(unique_df, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  #       message(sprintf("Saved: %s (%d rows, %d unique genes)", 
  #                      file_name, nrow(unique_df), 
  #                      length(unique(unique_df[[gene_column]]))))
  #     }
  #   }
  # }
  
# Save unique genes for each group
save_unique_genes_all <- function(split_data, prefix = output_prefix) {
  # Set working directory to output_dir
  old_wd <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_wd))  # Reset on exit
  
  for (group in names(split_data)) {
    df_group <- split_data[[group]]
    if (nrow(df_group) > 0) {
      unique_df <- df_group %>% 
        distinct(across(all_of(gene_column)), .keep_all = TRUE) %>% 
        select(any_of(intersect(selected_columns, names(df_group))))
      
      # Make sure file_name is a single string
      file_name <- paste0(prefix, tolower(group), "_", dataset_name, ".tsv")
      
      # Debug print
      message(sprintf("Debug: file_name = '%s' (length: %d)", 
                      file_name, length(file_name)))
      
      # Ensure it's a single string
      if (length(file_name) > 1) {
        file_name <- file_name[1]
        message("Warning: file_name was a vector, using first element")
      }
      
      # Write the file
      write.table(unique_df, file = file_name, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
      message(sprintf("Saved: %s (%d rows, %d unique genes)", 
                     file_name, nrow(unique_df), 
                     length(unique(unique_df[[gene_column]]))))
    }
  }
}

  # Filter dataset by gene lists
  filter_by_gene_lists <- function(split_data, gene_lists, gene_col = gene_column) {
    results <- list()
    
    for (group in names(split_data)) {
      df_group <- split_data[[group]]
      
      if (nrow(df_group) > 0) {
        for (list_name in names(gene_lists)) {
          genes <- gene_lists[[list_name]]
          if (length(genes) > 0) {
            filtered <- df_group %>%
              filter(!!sym(gene_col) %in% genes) %>%
              mutate(type = list_name, group = group) %>% 
              select(any_of(intersect(selected_columns_gene_list, names(df_group))))
            
            if (nrow(filtered) > 0) {
              key <- paste(dataset_name, group, list_name, sep = "_")
              results[[key]] <- filtered
            }
          }
        }
      }
    }
    
    return(results)
  }
  

  # Write list of dataframes to Excel with safe sheet names
  write_list_to_excel <- function(data_list, filename) {
    # Set working directory to output_dir
    old_wd <- getwd()
    setwd(output_dir)
    on.exit(setwd(old_wd))  # Reset on exit

    wb <- createWorkbook()
    
    # Filter out empty dataframes
    data_list <- data_list[sapply(data_list, nrow) > 0]
    
    if (length(data_list) == 0) {
      message(sprintf("No data to write for %s", filename))
      return(NULL)
    }
    
    original_names <- names(data_list)
    
    # Generate safe sheet names: truncate to 31 chars and make unique
    safe_names <- substr(original_names, 1, 31)
    safe_names <- make.unique(safe_names, sep = "_")
    
    # Create mapping: original name → safe sheet name
    name_mapping <- data.frame(
      Original_Name = original_names,
      Sheet_Name = safe_names,
      Total_Rows = sapply(data_list, nrow),
      Unique_Genes = sapply(data_list, function(df) {
        if (nrow(df) > 0 && gene_column %in% names(df)) {
          n_distinct(df[[gene_column]])
        } else {
          0
        }
      }),
      stringsAsFactors = FALSE
    )
    
    # Write summary sheet
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", name_mapping)
    
    # Write each dataframe with safe sheet name
    for (i in seq_along(data_list)) {
      if (nrow(data_list[[i]]) > 0) {
        addWorksheet(wb, safe_names[i])
        writeData(wb, safe_names[i], data_list[[i]])
      }
    }
    
    saveWorkbook(wb, filename, overwrite = TRUE)
    message(sprintf("Saved Excel file: %s", filename))
    message(sprintf("\n"))
  }
  
  # Filter to unique genes
  filter_unique_genes <- function(df) {
    if (nrow(df) > 0 && gene_column %in% names(df)) {
      df %>% distinct(across(all_of(gene_column)), .keep_all = TRUE)
    } else {
      df
    }
  }
  
  # --- Main Processing -----------------------------------------------------
  message(sprintf("\n==== Excel reports ===\n" ))
  # message(sprintf("Processing dataset: %s", dataset_name))
  message(sprintf("Total rows: %d", nrow(df)))
  message(sprintf("Unique groups in '%s': %s", group_column, 
                  paste(unique(df[[group_column]]), collapse = ", ")))
  
  # 1. Split dataframe by group
  split_data <- split_by_group(df)
  message(sprintf("Split into %d groups", length(split_data)))
  
  # 2. Save unique genes to TSV files
  save_unique_genes_all(split_data)

  # 3. Filter by gene lists
  filtered_results <- filter_by_gene_lists(split_data, gene_lists)
  message(sprintf("Created %d filtered datasets", length(filtered_results)))
  
  # 4. Create combined filtered dataset
  combined_filtered <- list()
  if (length(filtered_results) > 0) {
    combined_df <- bind_rows(filtered_results)
    if (nrow(combined_df) > 0) {
      combined_filtered[[paste0("combined_", dataset_name)]] <- combined_df
      message(sprintf("Combined filtered dataset: %d rows", nrow(combined_df)))
    }
  }
  
  # 5. Create unique versions
  unique_split_list <- lapply(split_data, filter_unique_genes)
  
  if (length(combined_filtered) > 0) {
    unique_combined_filtered <- lapply(combined_filtered, filter_unique_genes)
  } else {
    unique_combined_filtered <- list()
  }
  
  # 6. Write to Excel files
  if (length(split_data) > 0) {
    write_list_to_excel(split_data, sprintf("%s_grouped_datasets.xlsx", dataset_name))
  }
  
  if (length(unique_split_list) > 0) {
    write_list_to_excel(unique_split_list, sprintf("%s_unique_grouped_datasets.xlsx", dataset_name))
  }
  
  if (length(combined_filtered) > 0) {
    write_list_to_excel(combined_filtered, sprintf("%s_combined_filtered_by_gene_lists.xlsx", dataset_name))
  }
  
  if (length(unique_combined_filtered) > 0) {
    write_list_to_excel(unique_combined_filtered, sprintf("%s_unique_combined_filtered_by_gene_lists.xlsx", dataset_name))
  }
  
  # Return processed data invisibly
  invisible(list(
    split_data = split_data,
    filtered_results = filtered_results,
    combined_filtered = combined_filtered,
    unique_split_list = unique_split_list,
    unique_combined_filtered = unique_combined_filtered
  ))
}



#!excel for enrihcment 
create_enrichment_excel <- function(
  datasets = NULL,
  output_dir = getwd(),
  gene_lists = NULL,
  group_column = "sample_label",
  gene_column = "SYMBOL",
  dedup_columns = c("SYMBOL", "CHROM", "REF", "ALT", "POS", "SAMPLES"),
  selected_columns = c("CHROM", "POS", "REF", "ALT", "Allele", "IMPACT", "SYMBOL", 
                      "Consequence", "am_pathogenicity", "LoF", "CANONICAL", "Gene",
                      "BIOTYPE", "Feature", "SAMPLES", "GENOTYPES", "mis.oe_ci.upper",
                      "lof.oe_ci.upper", "cds_length", "num_coding_exons",
                      "sample_label_2", "sample_label", "sample_label_3", "count")
) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Set working directory
  setwd(output_dir)
  
  # Load required packages
  if (!require(dplyr)) stop("dplyr package is required")
  if (!require(openxlsx)) stop("openxlsx package is required")
  
  # --- DEFAULT DATASETS (if not provided) ---
  if (is.null(datasets)) {
    datasets <- list(
      Missense = Missense,
      Missense_canonical = Missense_canonical,
      MPC_only = MPC_only,
      PTV_canonic = PTV_canonic
    )
  }
  
  # --- DEFAULT GENE LISTS (if not provided) ---
  if (is.null(gene_lists)) {
    gene_lists <- list(
      bipolar = genes_bipolar,
      schema_pval = genes_schema_pval,
      schema_qval = genes_schema_qval,
      gwas = genes_gwas,
      brain_ntm = unique(brain_gene_consensus_ntm_consensus_no_pitular$Gene.name),
      brain_2sd = unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene.name)
    )
  }
  
  # --- VALIDATE INPUTS ---
  if (!is.list(datasets) || length(datasets) == 0) {
    stop("datasets must be a non-empty named list of data frames")
  }
  
  if (!is.list(gene_lists) || length(gene_lists) == 0) {
    stop("gene_lists must be a non-empty named list of gene vectors")
  }
  
  # Check if required columns exist in datasets
  required_cols <- c(group_column, gene_column, dedup_columns, "sample_label_2")
  for (ds_name in names(datasets)) {
    missing_cols <- setdiff(required_cols, names(datasets[[ds_name]]))
    if (length(missing_cols) > 0) {
      warning(sprintf("Dataset '%s' missing columns: %s", ds_name, paste(missing_cols, collapse = ", ")))
    }
  }
  
  # --- CORE FUNCTIONS ---
  
  # Split dataset by group
  split_by_group <- function(df, group_col = group_column) {
    df %>%
      split(.[[group_col]]) %>%
      lapply(function(x) {
        x %>% 
          select(any_of(selected_columns)) %>%
          distinct(across(all_of(dedup_columns)), .keep_all = TRUE)
      })
  }
  
  # Create safe sheet names
  create_safe_sheet_name <- function(full_name) {
    # Use abbreviations
    short_name <- gsub("pure", "p", full_name)
    short_name <- gsub("mixed", "m", short_name)
    short_name <- gsub("bipolar", "bip", short_name)
    short_name <- gsub("schema_pval", "spv", short_name)
    short_name <- gsub("schema_qval", "sqv", short_name)
    short_name <- gsub("brain_ntm", "ntm", short_name)
    short_name <- gsub("brain_2sd", "2sd", short_name)
    short_name <- gsub("all genes", "all", short_name)
    short_name <- gsub("Control", "Ctrl", short_name)
    short_name <- gsub("Case", "Cs", short_name)
    short_name <- gsub("Converter", "Conv", short_name)
    short_name <- gsub("Non_Converter", "NonC", short_name)
    short_name <- gsub("SCZ", "SCZ", short_name)
    short_name <- gsub("BD", "BD", short_name)
    short_name <- gsub("\\+", "&", short_name)  # Replace + with & for better readability
    
    # Ensure length is <= 31
    substr(short_name, 1, 31)
  }
  
  # Main function to create Excel sheets for a single dataset
  create_dataset_excel <- function(dataset_name, dataset, gene_lists, output_dir) {
    message("Processing dataset: ", dataset_name)
    
    # Create workbook for this specific dataset
    wb <- createWorkbook()
    
    # Split the dataset
    split_dataset <- split_by_group(dataset)
    
    # Get all available groups from group_column
    available_groups <- names(split_dataset)
    message("Available groups in ", group_column, ": ", paste(available_groups, collapse = ", "))
    
    # Initialize sheets list and summary data
    sheets_list <- list()
    summary_data <- data.frame()
    
    # Step 1: Create pure sheets for each group
    for (group in available_groups) {
      group_data <- split_dataset[[group]]
      
      # For PURE sheets (distinct genes) for each gene list
      for (gene_list_name in names(gene_lists)) {
        gene_list <- gene_lists[[gene_list_name]]
        
        pure_sheet_name <- paste(group, "pure", gene_list_name)
        sheet_data <- group_data %>%
          filter(!!sym(gene_column) %in% gene_list) %>%
          distinct(!!sym(gene_column), .keep_all = TRUE) %>%
          select(any_of(selected_columns))
        
        sheets_list[[pure_sheet_name]] <- sheet_data
        
        # Add to summary
        safe_name <- create_safe_sheet_name(pure_sheet_name)
        summary_data <- bind_rows(summary_data, data.frame(
          Original_Name = pure_sheet_name,
          Sheet_Name = safe_name,
          Rows = nrow(sheet_data),
          Unique_Genes = n_distinct(sheet_data[[gene_column]]),
          Group = group,
          Type = "Pure",
          Gene_List = gene_list_name,
          stringsAsFactors = FALSE
        ))
      }
      
      # Pure all genes sheet (no filtering)
      pure_all_sheet_name <- paste(group, "pure all genes")
      sheet_data <- group_data %>%
        distinct(!!sym(gene_column), .keep_all = TRUE) %>%
        select(any_of(selected_columns))
      
      sheets_list[[pure_all_sheet_name]] <- sheet_data
      
      # Add to summary
      safe_name <- create_safe_sheet_name(pure_all_sheet_name)
      summary_data <- bind_rows(summary_data, data.frame(
        Original_Name = pure_all_sheet_name,
        Sheet_Name = safe_name,
        Rows = nrow(sheet_data),
        Unique_Genes = n_distinct(sheet_data[[gene_column]]),
        Group = group,
        Type = "Pure",
        Gene_List = "All genes",
        stringsAsFactors = FALSE
      ))
    }
    
    # Step 2: Handle mixed group - collapse sample_label_2 and match to groups
    if ("mixed" %in% available_groups) {
      # Get all mixed data
      mixed_data <- dataset %>% filter(!!sym(group_column) == "mixed")
      
      # For each sample in mixed group, check which group it belongs to based on sample_label_2
      mixed_samples_by_group <- list()
      
      for (group in available_groups) {
        if (group != "mixed") {
          # Find samples in mixed group that have this group in their sample_label_2
          group_mixed_samples <- mixed_data %>%
            filter(grepl(group, sample_label_2, ignore.case = TRUE))
          
          if (nrow(group_mixed_samples) > 0) {
            mixed_samples_by_group[[group]] <- group_mixed_samples
            message("Found ", n_distinct(group_mixed_samples$SAMPLES), 
                   " samples in mixed group for ", group)
          }
        }
      }
      
      # Step 3: Create pure + mixed sheets for each group
      for (group in names(mixed_samples_by_group)) {
        group_mixed_data <- mixed_samples_by_group[[group]]
        group_pure_data <- split_dataset[[group]]
        
        # Combine pure group data with matched mixed data
        group_combined_data <- bind_rows(group_pure_data, group_mixed_data)
        
        # For each gene list, create pure + mixed sheets
        for (gene_list_name in names(gene_lists)) {
          gene_list <- gene_lists[[gene_list_name]]
          
          # Pure + mixed sheet (distinct genes)
          pure_mixed_sheet_name <- paste(group, "pure + mixed", gene_list_name)
          sheet_data <- group_combined_data %>%
            filter(!!sym(gene_column) %in% gene_list) %>%
            distinct(!!sym(gene_column), .keep_all = TRUE) %>%
            select(any_of(selected_columns))
          
          sheets_list[[pure_mixed_sheet_name]] <- sheet_data
          
          # Add to summary
          safe_name <- create_safe_sheet_name(pure_mixed_sheet_name)
          summary_data <- bind_rows(summary_data, data.frame(
            Original_Name = pure_mixed_sheet_name,
            Sheet_Name = safe_name,
            Rows = nrow(sheet_data),
            Unique_Genes = n_distinct(sheet_data[[gene_column]]),
            Group = group,
            Type = "Pure + Mixed",
            Gene_List = gene_list_name,
            stringsAsFactors = FALSE
          ))
        }
        
        # Pure + mixed all genes sheet (no filtering)
        pure_mixed_all_sheet_name <- paste(group, "pure + mixed all genes")
        sheet_data <- group_combined_data %>%
          distinct(!!sym(gene_column), .keep_all = TRUE) %>%
          select(any_of(selected_columns))
        
        sheets_list[[pure_mixed_all_sheet_name]] <- sheet_data
        
        # Add to summary
        safe_name <- create_safe_sheet_name(pure_mixed_all_sheet_name)
        summary_data <- bind_rows(summary_data, data.frame(
          Original_Name = pure_mixed_all_sheet_name,
          Sheet_Name = safe_name,
          Rows = nrow(sheet_data),
          Unique_Genes = n_distinct(sheet_data[[gene_column]]),
          Group = group,
          Type = "Pure + Mixed",
          Gene_List = "All genes",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Create safe sheet names and ensure uniqueness
    safe_sheet_names <- sapply(names(sheets_list), create_safe_sheet_name)
    safe_sheet_names <- make.unique(safe_sheet_names, sep = "_")
    
    # Update summary with final safe names
    for (i in seq_along(safe_sheet_names)) {
      orig_name <- names(sheets_list)[i]
      safe_name <- safe_sheet_names[i]
      summary_data$Sheet_Name[summary_data$Original_Name == orig_name] <- safe_name
    }
    
    # --- CREATE EXCEL FILE ---
    
    # Add summary sheet as FIRST sheet
    addWorksheet(wb, "Summary", tabColour = "lightblue")
    writeData(wb, "Summary", summary_data)
    
    # Add each data sheet
    for (i in seq_along(sheets_list)) {
      addWorksheet(wb, safe_sheet_names[i])
      writeData(wb, safe_sheet_names[i], sheets_list[[i]])
    }
    
    # Save workbook with dataset name as prefix
    output_file <- paste0(dataset_name, "_enrichment_results.xlsx")
    saveWorkbook(wb, output_file, overwrite = TRUE)
    
    message("Created: ", output_file)
    message("Total sheets: ", length(sheets_list))
    
    return(list(
      file = output_file,
      sheets_created = length(sheets_list),
      summary = summary_data
    ))
  }
  
  # --- MAIN PROCESSING ---
  
  # Process each dataset and create separate Excel files
  results <- list()
  for (dataset_name in names(datasets)) {
    result <- create_dataset_excel(
      dataset_name = dataset_name,
      dataset = datasets[[dataset_name]],
      gene_lists = gene_lists,
      output_dir = output_dir
    )
    results[[dataset_name]] <- result
  }
  
  # Final summary
  message("\n=== PROCESSING COMPLETE ===")
  message("Output directory: ", output_dir)
  message("Files created:")
  for (dataset_name in names(results)) {
    message("  - ", results[[dataset_name]]$file, " (", results[[dataset_name]]$sheets_created, " sheets)")
  }
  
  return(invisible(results))
}
