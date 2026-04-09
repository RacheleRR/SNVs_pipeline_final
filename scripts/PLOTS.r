#plots 

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(ggsignif)
library(purrr)
library(ggrepel)
library(cowplot)  # For better themes and color palettes



# Main plotting function with manifest-based detection
create_statistical_plots <- function(results_df, 
                                    manifest = NULL,
                                    plot_type = c("volcano", "parallel", "dot", "manhattan", "forest"),
                                    significant_only = FALSE,
                                    max_labels = 20,
                                    p_value_threshold = 0.05,
                                    p_value_col = NULL,
                                    color_palette = "Set2",
                                    base_size = 14) {
    
    # Check for required columns
    required_cols <- c("p_value", "Row")
    if (!all(required_cols %in% names(results_df))) {
        stop("Results dataframe must contain at least 'p_value' and 'Row' columns")
    }

    # SIMPLE: Determine if this is pairwise by checking Group1 and Group2 columns
    is_pairwise <- FALSE
    
    if ("Group1" %in% names(results_df) && "Group2" %in% names(results_df)) {
        # Get ALL unique group names from BOTH columns
        all_groups <- unique(c(results_df$Group1, results_df$Group2))
        
        # If we have MORE THAN 2 unique groups total, it's pairwise comparisons
        if (length(all_groups) > 2) {
            is_pairwise <- TRUE
            
            # Create Comparison column with alphabetical ordering
            results_df <- results_df %>%
                rowwise() %>%
                mutate(
                    # Sort group names alphabetically for consistent "vs" format
                    sorted = list(sort(c(Group1, Group2))),
                    Comparison = paste(sorted[1], "vs", sorted[2])
                ) %>%
                ungroup() %>%
                select(-sorted) %>%
                mutate(
                    # Make it a factor with alphabetical ordering
                    Comparison = factor(Comparison, levels = sort(unique(Comparison)))
                )
        }
        # If exactly 2 unique groups total, it's a simple 2-group comparison
        else if (length(all_groups) == 2) {
            results_df <- results_df %>%
                mutate(Comparison = paste(all_groups[1], "vs", all_groups[2]))
        }
    }
    
    # Determine p-value column to use
    if (is.null(p_value_col)) {
        # Determine p-value column to use (prioritize adjusted p-values)
        p_value_col <- ifelse("p_adj_global" %in% names(results_df), "p_adj_global", 
                            ifelse("p_adj_within" %in% names(results_df), "p_adj_within", "p_value"))
    }
    
    # Filter significant results if requested
    if (significant_only) {
        n_before <- nrow(results_df)
        results_df <- results_df %>% 
            filter(.data[[p_value_col]] < p_value_threshold)
        
        if (nrow(results_df) == 0) {
            message("No significant results to plot (", p_value_col, " < ", p_value_threshold, ")")
            return(NULL)
        }
        message("Filtered to ", nrow(results_df), " significant results (from ", n_before, ")")
    }
    
    # Create the requested plot type
    plot_type <- match.arg(plot_type)
    
    plot_func <- switch(plot_type,
           volcano = create_volcano_plot_cowplot,
           parallel = create_parallel_plot_cowplot,
           dot = create_dot_plot_cowplot,
           manhattan = create_manhattan_plot_cowplot,
           forest = create_forest_plot_cowplot
    )
    
    # Create the plot
    plot_func(results_df, is_pairwise, p_value_col, p_value_threshold, 
             max_labels, color_palette, base_size)
}

# Helper function to get colors from cowplot
get_cowplot_colors <- function(n, palette = "Set2") {
  if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    n_request <- max(3, n)  # request at least 3 to avoid warning
    if (n_request <= max_colors) {
      return(RColorBrewer::brewer.pal(n_request, palette)[1:n])  # trim to actual n needed
    }
  }
  
  default_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                      "#0072B2", "#D55E00", "#CC79A7", "#999999")
  if (n <= length(default_colors)) {
    return(default_colors[1:n])
  }
  return(viridis::viridis(n))
}

# Volcano plot with cowplot styling
create_volcano_plot_cowplot <- function(results_df, is_pairwise, p_value_col, p_value_threshold, 
                                       max_labels, color_palette, base_size) {
    
    # Prepare plot data
    plot_data <- results_df %>%
        mutate(
            log_p = -log10(.data[[p_value_col]]),
            is_significant = .data[[p_value_col]] < p_value_threshold
        )
    
    # Add effect size - ONLY OR
    if ("OR" %in% names(plot_data)) {
        plot_data$effect_size <- log2(plot_data$OR)
        x_lab <- "log₂(Odds Ratio)"
    } else {
        stop("OR column required for volcano plot")
    }
    
    # Determine coloring and labeling
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        # For pairwise comparisons, use both color AND shape
        color_var <- "Comparison"
        shape_var <- "Comparison"  # Use shape for different comparisons
        n_comparisons <- length(unique(plot_data$Comparison))
        
        title_text <- "Pairwise Associations Volcano Plot"
        
        # Limit shapes to available ones (max 6 distinct shapes)
        if (n_comparisons > 6) {
            warning("More than 6 comparisons: some shapes will be repeated")
            # Use color as primary, shape only for first 6
            shape_var <- ifelse(as.numeric(factor(plot_data$Comparison)) <= 6, 
                               "Comparison", "1")
        }
        
    } else if ("Dataset" %in% names(plot_data) && length(unique(plot_data$Dataset)) > 1) {
        # Multiple datasets
        color_var <- "Dataset"
        shape_var <- "is_significant"
        n_colors <- length(unique(plot_data$Dataset))
        title_text <- "Volcano Plot by Dataset"
    } else {
        # Single dataset or no grouping
        color_var <- "is_significant"
        shape_var <- "is_significant"
        n_colors <- 2
        title_text <- "Volcano Plot"
    }
    
    # Create base plot
    p <- ggplot(plot_data, aes(x = effect_size, y = log_p)) +
        theme_cowplot(font_size = base_size) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "grey90", size = 0.2)
        )
    
    # Add points with appropriate aesthetics
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        # For pairwise: color AND shape by comparison
        p <- p + geom_point(aes(color = Comparison, shape = Comparison, 
                               alpha = is_significant), 
                           size = 3, stroke = 0.5)
        
        # Use custom palette for colors
        n_comparisons <- length(unique(plot_data$Comparison))
        colors <- get_cowplot_colors(n_comparisons, color_palette)
        
        # Use different shapes (max 6 available)
        available_shapes <- c(16, 17, 15, 18, 8, 3)
        shapes_to_use <- available_shapes[1:min(n_comparisons, 6)]
        if (n_comparisons > 6) {
            shapes_to_use <- c(shapes_to_use, rep(16, n_comparisons - 6))
        }
        
        p <- p +
            scale_color_manual(values = colors) +
            scale_shape_manual(values = shapes_to_use)
            
    } else if (color_var == "is_significant") {
        # Simple significance coloring
        p <- p + geom_point(aes(color = is_significant, shape = is_significant), 
                           size = 3, stroke = 0.5) +
            scale_color_manual(
                values = c("FALSE" = "grey70", "TRUE" = "#E41A1C"),
                labels = c("FALSE" = "Not significant", "TRUE" = paste0("p < ", p_value_threshold))
            ) +
            scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17))
    } else {
        # Color by dataset
        p <- p + geom_point(aes_string(color = color_var, shape = shape_var), 
                           size = 3, stroke = 0.5)
        
        n_colors <- length(unique(plot_data[[color_var]]))
        colors <- get_cowplot_colors(n_colors, color_palette)
        p <- p + scale_color_manual(values = colors) +
            scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17))
    }
    
    p <- p +
        scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 0.9)) +
        guides(alpha = "none")
    
    # Add labels for significant points
    if (any(plot_data$is_significant)) {
        sig_data <- plot_data %>% filter(is_significant)
        
        # Limit labels if too many
        if (nrow(sig_data) > max_labels) {
            sig_data <- sig_data %>% 
                arrange(desc(log_p)) %>%
                slice_head(n = max_labels)
            message("Labeling top ", max_labels, " most significant points")
        }
        
        # Use Row column (no cleaning)
        label_col <- "Row"
        
        p <- p + ggrepel::geom_text_repel(
            data = sig_data,
            aes_string(label = label_col),
            max.overlaps = max_labels,
            size = base_size / 3.5,
            box.padding = 0.35,
            point.padding = 0.3,
            segment.color = "grey50",
            segment.alpha = 0.5,
            min.segment.length = 0.2
        )
    }
    
    # Add significance threshold line
    p <- p + 
        geom_hline(yintercept = -log10(p_value_threshold), 
                  linetype = "dashed", 
                  color = "red3", 
                  alpha = 0.7,
                  size = 0.8) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", size = 0.5) +
        labs(
            title = title_text,
            x = x_lab,
            y = paste0("-log₁₀(", p_value_col, ")")
        )
    
    # Adjust legend
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        p <- p + guides(color = guide_legend(title = "Comparison", nrow = 2, byrow = TRUE),
                       shape = guide_legend(title = "Comparison", nrow = 2, byrow = TRUE))
    } else if (color_var == "is_significant") {
        p <- p + labs(color = "Significance", shape = "Significance")
    } else {
        p <- p + labs(color = color_var)
        if (shape_var == "is_significant") {
            p <- p + guides(shape = "none")  # Hide shape legend if it's just significance
        }
    }
    
    return(p)
}

# Dot plot with cowplot styling - ADDED LABELS FOR SIGNIFICANT POINTS
create_dot_plot_cowplot <- function(results_df, is_pairwise, p_value_col, p_value_threshold, 
                                   max_labels, color_palette, base_size) {
    
    # Prepare plot data
    plot_data <- results_df %>%
        mutate(
            log_p = -log10(.data[[p_value_col]]),
            is_significant = .data[[p_value_col]] < p_value_threshold
        )
    
    # Order rows by p-value (most significant first) - using original Row
    plot_data <- plot_data %>%
        mutate(Row = factor(Row, levels = unique(Row[order(log_p, decreasing = TRUE)])))
    
    # Determine coloring and shaping
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        # For pairwise: color and shape by comparison
        color_var <- "Comparison"
        shape_var <- "Comparison"
        n_comparisons <- length(unique(plot_data$Comparison))
        title_text <- "Pairwise Comparison Significance"
        
        # Limit shapes to available ones
        if (n_comparisons > 6) {
            shape_var <- ifelse(as.numeric(factor(plot_data$Comparison)) <= 6, 
                               "Comparison", "1")
        }
        
    } else if ("Dataset" %in% names(plot_data) && length(unique(plot_data$Dataset)) > 1) {
        color_var <- "Dataset"
        shape_var <- "is_significant"
        n_colors <- length(unique(plot_data$Dataset))
        title_text <- "Statistical Test Results by Dataset"
    } else {
        color_var <- "is_significant"
        shape_var <- "is_significant"
        n_colors <- 2
        title_text <- "Statistical Test Results"
    }
    
    # Create base plot
    p <- ggplot(plot_data, aes(x = log_p, y = Row)) +
        theme_cowplot(font_size = base_size) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = base_size * 0.8),
            panel.grid.major.x = element_line(color = "grey90", size = 0.2)
        )
    
    # Add points with appropriate aesthetics
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        p <- p + geom_point(aes(color = Comparison, shape = Comparison, 
                               size = is_significant), 
                           alpha = 0.8)
        
        # Use custom palette for colors
        n_comparisons <- length(unique(plot_data$Comparison))
        colors <- get_cowplot_colors(n_comparisons, color_palette)
        
        # Use different shapes
        available_shapes <- c(16, 17, 15, 18, 8, 3)
        shapes_to_use <- available_shapes[1:min(n_comparisons, 6)]
        if (n_comparisons > 6) {
            shapes_to_use <- c(shapes_to_use, rep(16, n_comparisons - 6))
        }
        
        p <- p +
            scale_color_manual(values = colors) +
            scale_shape_manual(values = shapes_to_use) +
            scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3)) +
            guides(size = "none")
            
    } else if (color_var == "is_significant") {
        p <- p + geom_point(aes(color = is_significant, shape = is_significant, 
                               size = is_significant), 
                           alpha = 0.8) +
            scale_color_manual(
                values = c("FALSE" = "grey70", "TRUE" = "#E41A1C"),
                labels = c("FALSE" = "Not significant", "TRUE" = paste0("p < ", p_value_threshold))
            ) +
            scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
            scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3)) +
            guides(size = "none")
    } else {
        p <- p + geom_point(aes_string(color = color_var, shape = shape_var,
                                      size = "is_significant"), 
                           alpha = 0.8)
        
        n_colors <- length(unique(plot_data[[color_var]]))
        colors <- get_cowplot_colors(n_colors, color_palette)
        p <- p + scale_color_manual(values = colors) +
            scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
            scale_size_manual(values = c("FALSE" = 2, "TRUE" = 3)) +
            guides(size = "none")
    }
    
    # ADD LABELS FOR SIGNIFICANT POINTS (like in volcano plot)
    if (any(plot_data$is_significant)) {
        sig_data <- plot_data %>% filter(is_significant)
        
        # Limit labels if too many
        if (nrow(sig_data) > max_labels) {
            sig_data <- sig_data %>% 
                arrange(desc(log_p)) %>%
                slice_head(n = max_labels)
            message("Labeling top ", max_labels, " most significant points in dot plot")
        }
        
        # Add text labels for significant points
        p <- p + geom_text(
            data = sig_data,
            aes(label = Row, x = log_p),
            hjust = -0.1,
            vjust = 0.5,
            size = base_size / 3.5,
            check_overlap = TRUE
        )
    }
    
    # Add significance threshold line
    p <- p +
        geom_vline(xintercept = -log10(p_value_threshold), 
                  linetype = "dashed", 
                  color = "red3", 
                  alpha = 0.7,
                  size = 0.8) +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +  # More space for labels
        labs(
            title = title_text,
            x = paste0("-log₁₀(", p_value_col, ")"),
            y = "Feature"
        )
    
    # Adjust legend
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        p <- p + guides(color = guide_legend(title = "Comparison", nrow = 2, byrow = TRUE),
                       shape = guide_legend(title = "Comparison", nrow = 2, byrow = TRUE))
    } else if (color_var == "is_significant") {
        p <- p + labs(color = "Significance", shape = "Significance")
    } else {
        p <- p + labs(color = color_var)
        if (shape_var == "is_significant") {
            p <- p + guides(shape = "none")
        }
    }
    
    return(p)
}

# Manhattan plot - UPDATED TITLE AND ADDED LABELS
create_manhattan_plot_cowplot <- function(results_df, is_pairwise, p_value_col, p_value_threshold, 
                                         max_labels, color_palette, base_size) {
    
    # Prepare plot data
    plot_data <- results_df %>%
        mutate(
            log_p = -log10(.data[[p_value_col]]),
            is_significant = .data[[p_value_col]] < p_value_threshold
        )
    
    # Create x-axis variable
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        # Alphabetical ordering of comparisons
        plot_data <- plot_data %>%
            mutate(x_var = factor(Comparison, levels = sort(unique(Comparison))))
        x_lab <- "Group Comparisons"
        # UPDATED TITLE: Different for pairwise
        title_text <- "Pairwise Manhattan Plot"
    } else {
        # Use original Row values
        plot_data <- plot_data %>%
            mutate(x_var = factor(Row, levels = unique(Row)))
        x_lab <- "Features"
        title_text <- "Manhattan Plot"
    }
    
    # Determine coloring
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        color_var <- "Comparison"
        n_colors <- length(unique(plot_data$Comparison))
    } else if ("Dataset" %in% names(plot_data) && length(unique(plot_data$Dataset)) > 1) {
        color_var <- "Dataset"
        n_colors <- length(unique(plot_data$Dataset))
        title_text <- "Manhattan Plot by Dataset"
    } else {
        color_var <- "is_significant"
        n_colors <- 2
    }
    
    # Create plot
    p <- ggplot(plot_data, aes(x = x_var, y = log_p)) +
        geom_point(aes_string(color = color_var), 
                  size = 2.5, 
                  alpha = 0.7) +
        geom_hline(yintercept = -log10(p_value_threshold), 
                  linetype = "dashed", 
                  color = "red3", 
                  alpha = 0.7,
                  size = 0.8) +
        theme_cowplot(font_size = base_size) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
            panel.grid.major.y = element_line(color = "grey90", size = 0.2)
        )
    
    # Set colors
    if (color_var == "is_significant") {
        p <- p + scale_color_manual(
            values = c("FALSE" = "grey60", "TRUE" = "#E41A1C"),
            labels = c("FALSE" = "Not significant", "TRUE" = paste0("p < ", p_value_threshold))
        )
    } else {
        colors <- get_cowplot_colors(n_colors, color_palette)
        p <- p + scale_color_manual(values = colors)
    }
    
    # ADD LABELS FOR SIGNIFICANT POINTS IN MANHATTAN PLOT
    if (any(plot_data$is_significant)) {
        sig_data <- plot_data %>% filter(is_significant)
        
        # Limit labels if too many
        if (nrow(sig_data) > max_labels) {
            sig_data <- sig_data %>% 
                arrange(desc(log_p)) %>%
                slice_head(n = max_labels)
            message("Labeling top ", max_labels, " most significant points in Manhattan plot")
        }
        
        # Add text labels for significant points
        p <- p + geom_text(
            data = sig_data,
            aes(label = Row, y = log_p),
            hjust = -0.1,
            size = base_size / 4,
            check_overlap = TRUE
        )
    }

    p <- p +
        labs(
            title = title_text,  # UPDATED: Different title for pairwise
            x = x_lab,
            y = paste0("-log₁₀(", p_value_col, ")"),
            color = ifelse(color_var == "is_significant", "Significance", color_var)
        )
    
    return(p)
}

# Parallel coordinates plot (for pairwise comparisons with effect sizes)
create_parallel_plot_cowplot <- function(results_df, is_pairwise, p_value_col, p_value_threshold, 
                                        max_labels, color_palette, base_size) {
    
    # Check for required columns
    if (!"OR" %in% names(results_df)) {
        stop("Parallel plot requires OR column")
    }
    
    # Prepare plot data
    plot_data <- results_df %>%
        mutate(
            log_OR = log2(OR),
            is_significant = .data[[p_value_col]] < p_value_threshold
        )
    
    # Ensure Comparison is factor with alphabetical ordering
    if (is_pairwise && "Comparison" %in% names(plot_data)) {
        plot_data <- plot_data %>%
            mutate(Comparison = factor(Comparison, levels = sort(unique(Comparison))))
    } else {
        stop("Parallel plot requires pairwise comparisons")
    }
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Comparison, y = log_OR)) +
        geom_line(aes(group = Row, color = Row, alpha = is_significant), 
                 size = 0.8) +
        geom_point(aes(color = Row, size = -log10(.data[[p_value_col]]), alpha = is_significant)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", size = 0.5) +
        scale_size_continuous(range = c(1.5, 4),
                             name = paste0("-log₁₀(", p_value_col, ")")) +
        scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8)) +
        guides(alpha = "none") +
        theme_cowplot(font_size = base_size) +
        theme(
            legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.y = element_line(color = "grey90", size = 0.2)
        ) +
        labs(
            title = "Effect Size Patterns Across Comparisons",
            x = "Group Comparison",
            y = "log₂(Odds Ratio)",
            color = "Feature"
        )
    
    # Color features using custom palette
    n_features <- length(unique(plot_data$Row))
    colors <- get_cowplot_colors(n_features, color_palette)
    p <- p + scale_color_manual(values = colors)
    
    # Facet by dataset if available
    if ("Dataset" %in% names(plot_data) && length(unique(plot_data$Dataset)) > 1) {
        p <- p + facet_wrap(~ Dataset, scales = "free_x", ncol = 2)
    }
    
    return(p)
}
create_forest_plot_cowplot <- function(results_df, is_pairwise, p_value_col, p_value_threshold,
                                       max_labels, color_palette, base_size) {
  
  if (!all(c("OR", "CI_low", "CI_high") %in% names(results_df))) {
    stop("Forest plot requires OR, CI_low, and CI_high columns")
  }
  
  not_sig_label <- paste0("Not significant (", p_value_col, " \u2265 ", p_value_threshold, ")")
  sig_label     <- paste0("Significant (", p_value_col, " < ",  p_value_threshold, ")")
  
  plot_data <- results_df %>%
    mutate(
      log_OR      = log2(OR),
      log_CI_low  = log2(CI_low),
      log_CI_high = log2(CI_high),
      is_significant      = .data[[p_value_col]] < p_value_threshold,
      significance_status = ifelse(is_significant, sig_label, not_sig_label),
      # make factor so scale only contains levels present in data
      significance_status = factor(significance_status, 
                                   levels = c(not_sig_label, sig_label))
    ) %>%
    mutate(Row = factor(Row, levels = unique(Row[order(log_OR, decreasing = TRUE)])))
  
  # Drop unused factor levels — critical for shape/fill scales
  plot_data$significance_status <- droplevels(plot_data$significance_status)
  levels_present <- levels(plot_data$significance_status)
  
  # Determine color variable
  if (is_pairwise && "Comparison" %in% names(plot_data)) {
    color_var    <- "Comparison"
    n_colors     <- length(unique(plot_data$Comparison))
    title_text   <- "Pairwise Forest Plot of Effect Sizes"
  } else if ("Dataset" %in% names(plot_data) && length(unique(plot_data$Dataset)) > 1) {
    color_var    <- "Dataset"
    n_colors     <- length(unique(plot_data$Dataset))
    title_text   <- "Forest Plot of Effect Sizes by Dataset"
  } else {
    color_var    <- "is_significant"
    n_colors     <- 2
    title_text   <- "Forest Plot of Effect Sizes"
  }
  
  # Build color palette
  if (color_var == "is_significant") {
    colors <- c("FALSE" = "grey70", "TRUE" = "#E41A1C")
  } else {
    colors <- get_cowplot_colors(n_colors, color_palette)
    names(colors) <- unique(plot_data[[color_var]])
  }
  
  # Shape scale — only levels that exist in the data
  all_shapes        <- c(21, 16)
  names(all_shapes) <- c(not_sig_label, sig_label)
  shape_values      <- all_shapes[levels_present]   # subset to present levels only
  
  # Fill scale
  fill_colors <- c("white" = "white", colors)
  
  # Base plot
  p <- ggplot(plot_data, aes(x = log_OR, y = Row)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.8) +
    theme_cowplot(font_size = base_size) +
    theme(
      legend.position      = "bottom",
      plot.title           = element_text(hjust = 0.5),
      axis.text.y          = element_text(size = base_size * 0.8),
      panel.grid.major.x   = element_line(color = "grey90", linewidth = 0.2)
    )
  
  # Points — stroke as mapped aesthetic to avoid vector-as-fixed issue
  p <- p +
    geom_point(
      aes(color = .data[[color_var]],
          shape = significance_status,
          fill  = ifelse(is_significant, .data[[color_var]], "white"),
          stroke = ifelse(is_significant, 0.5, 0.8)),   # now inside aes
      size     = 3,
      position = position_dodge(width = 0.3)
    ) +
    geom_errorbar(   # replaces deprecated geom_errorbarh
      aes(xmin  = log_CI_low,
          xmax  = log_CI_high,
          color = .data[[color_var]],
          alpha = is_significant),
      height   = 0.2,
      linewidth = 0.8,
      orientation = "y",
      position = position_dodge(width = 0.3)
    ) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8), guide = "none") +
    scale_shape_manual(
      values = shape_values,   # only levels present
      name   = "Point Style"
    ) +
    scale_fill_manual(values = fill_colors, guide = "none")
  
  # Color scale
  if (color_var == "is_significant") {
    p <- p + scale_color_manual(
      values = colors,
      labels = c("FALSE" = "Not significant", "TRUE" = paste0("p < ", p_value_threshold)),
      name   = "Significance"
    )
  } else {
    p <- p + scale_color_manual(values = colors, name = color_var)
  }
  
  p <- p + labs(
    title = title_text,
    x     = "log\u2082(Odds Ratio) with 95% CI",
    y     = "Feature"
  )
  
  # Guide overrides — built dynamically from levels_present
  shape_fill_override <- ifelse(levels_present == sig_label, "black", "white")
  
  p <- p + guides(
    shape = guide_legend(
      title         = "Point Style",
      override.aes  = list(
        color = "black",
        fill  = shape_fill_override   # same length as levels_present
      )
    )
  )
  
  if (is_pairwise && "Comparison" %in% names(plot_data)) {
    p <- p + guides(
      color = guide_legend(
        title        = "Comparison",
        override.aes = list(shape = 16, fill = colors, stroke = 0.5),
        nrow = 2, byrow = TRUE
      )
    )
  }
  
  return(p)
}

# Batch plotting function (removed comparison_order parameter)
create_comprehensive_visualization <- function(results_df, 
                                              manifest = NULL,
                                              plot_types = c("volcano", "forest","parallel", "dot", "manhattan"),
                                              significant_only = FALSE,
                                              save_plots = FALSE,
                                              output_dir = ".",
                                              prefix = "statistical_plots",
                                              p_value_threshold = 0.05,
                                              p_value_col = NULL,
                                              max_labels = 20,
                                              color_palette = "Set2",
                                              base_size = 14,
                                              width = 12,
                                              height = 8) {
    
    plots <- list()
    
    for (plot_type in plot_types) {
        cat("Creating", plot_type, "plot...\n")
        
        p <- create_statistical_plots(
            results_df = results_df,
            manifest = manifest,
            plot_type = plot_type,
            significant_only = significant_only,
            max_labels = max_labels,
            p_value_threshold = p_value_threshold,
            p_value_col = p_value_col,
            color_palette = color_palette,
            base_size = base_size
        )
        
        if (!is.null(p)) {
            plots[[plot_type]] <- p
            
            if (save_plots) {
                if (!dir.exists(output_dir)) {
                    dir.create(output_dir, recursive = TRUE)
                }
                
                filename <- file.path(output_dir, 
                                     paste0(prefix, "_", plot_type, ".png"))
                ggsave(filename, p, width = width, height = height, dpi = 300)
                message("Saved: ", filename)
            }
        }
    }
    
    return(plots)
}


#! Kruskal-Wallis/Wilcoxon visualizations - Unified with Fisher logic

# Main function for variant plots
create_variant_plots <- function(plot_data, 
                                results_table,
                                significant_only = TRUE,
                                p_value_threshold = 0.05,
                                p_value_col = NULL,
                                title_suffix = "",
                                metrics_per_page = 4,
                                color_palette = "Set2",
                                base_size = 12,
                                output_width = 14,
                                output_height = 10) {
    
    # Handle both single dataframe and list of dataframes
    if (is.list(plot_data) && !is.data.frame(plot_data)) {
        # It's a list of datasets - process each one
        plot_list <- plot_data
        all_plots <- list()
        
        for (dataset_name in names(plot_list)) {
            cat("Processing:", dataset_name, "\n")
            
            # Get the corresponding results for this dataset
            id_col <- intersect(c("Dataset", "variant_type"), names(results_table))[1]

            if (is.na(id_col)) {
                warning("No recognized ID column found in results_table")
                next
            }

            dataset_results <- results_table %>%
                filter(.data[[id_col]] == dataset_name)
            
            if (nrow(dataset_results) == 0) {
                warning("No results found for dataset: ", dataset_name)
                next
            }
            
            # Recursive call for single dataset
            plot_result <- create_variant_plots(
                plot_data = plot_list[[dataset_name]],
                results_table = dataset_results,
                significant_only = significant_only,
                p_value_threshold = p_value_threshold,
                p_value_col = p_value_col,
                title_suffix = paste(title_suffix, "-", dataset_name),
                metrics_per_page = metrics_per_page,
                color_palette = color_palette,
                base_size = base_size,
                output_width = output_width,
                output_height = output_height
            )
            
            if (!is.null(plot_result)) {
                all_plots[[dataset_name]] <- plot_result
            }
        }
        
        return(if (length(all_plots) > 0) all_plots else NULL)
        
    } else {
        # Single dataframe - proceed with plotting
        return(create_single_variant_plots(
            plot_data = plot_data,
            results_table = results_table,
            significant_only = significant_only,
            p_value_threshold = p_value_threshold,
            p_value_col = p_value_col,
            title_suffix = title_suffix,
            metrics_per_page = metrics_per_page,
            color_palette = color_palette,
            base_size = base_size,
            output_width = output_width,
            output_height = output_height
        ))
    }
}

# Core plotting function for single datasets
# Core plotting function for single datasets - FIXED FILTERING
create_single_variant_plots <- function(plot_data, 
                                       results_table, 
                                       significant_only = TRUE,
                                       p_value_threshold = 0.05,
                                       p_value_col = NULL,
                                       title_suffix = "",
                                       metrics_per_page = 4,
                                       color_palette = "Set2",
                                       base_size = 12,
                                       output_width = 14,
                                       output_height = 10) {
    
    # Determine group column (use Status if available, otherwise group)
    group_col <- if ("Status" %in% names(plot_data)) "Status" else "group"
    
    if (!group_col %in% names(plot_data)) {
        stop("Plot data must contain either 'Status' or 'group' column")
    }
    
    # Get number of groups from plot_data
    groups <- unique(na.omit(plot_data[[group_col]]))
    n_groups <- length(groups)
    
    if (n_groups < 2) {
        stop("Need at least 2 groups for statistical testing. Found: ", n_groups)
    }
    
    # AUTO-DETECT test type based on number of groups
    test_type <- ifelse(n_groups == 2, "wilcoxon", "kruskal")
    cat("Detected", n_groups, "groups. Using", test_type, "test visualization.\n")
    
    # Determine p-value column to use
    if (is.null(p_value_col)) {
        if (test_type == "wilcoxon") {
            p_value_col <- ifelse("global_primary_p_adj" %in% names(results_table), 
                                 "global_primary_p_adj", "p_value")
        } else {
            p_value_col <- ifelse("global_kruskal_p_adj" %in% names(results_table), 
                                 "global_kruskal_p_adj", "kruskal_p")
        }
    } else if (!p_value_col %in% names(results_table)) {
        stop("Specified p_value_col '", p_value_col, "' not found in results_table")
    }
    
    # FIXED: Apply significance filter correctly
    if (significant_only) {
        # For Kruskal-Wallis, we need to:
        # 1. Keep metrics where Kruskal p < threshold
        # 2. Keep ALL Dunn results for those metrics (even if they have NA in p_value_col)
        
        if (test_type == "kruskal") {
            # Get significant metrics from Kruskal results
            sig_kruskal_metrics <- results_table %>%
                filter(analysis_level == "kruskal") %>%
                filter(.data[[p_value_col]] < p_value_threshold) %>%
                pull(metric) %>%
                unique()
            
            if (length(sig_kruskal_metrics) == 0) {
                message("No significant Kruskal-Wallis results found in ", title_suffix, 
                       " (", p_value_col, " < ", p_value_threshold, ")")
                return(NULL)
            }
            
            # Keep all rows (Kruskal AND Dunn) for significant metrics
            filtered_results <- results_table %>%
                filter(metric %in% sig_kruskal_metrics)
            
            message("Found ", length(sig_kruskal_metrics), 
                   " significant metrics (", nrow(filtered_results), " total rows)")
            
        } else {
            # For Wilcoxon, simple filtering works
            filtered_results <- results_table %>%
                filter(.data[[p_value_col]] < p_value_threshold)
            
            if (nrow(filtered_results) == 0) {
                message("No significant Wilcoxon results found in ", title_suffix, 
                       " (", p_value_col, " < ", p_value_threshold, ")")
                return(NULL)
            }
            message("Found ", nrow(filtered_results), " significant Wilcoxon results")
        }
    } else {
        filtered_results <- results_table
    }
    
    # Get metrics to plot
    if (test_type == "kruskal") {
        # For Kruskal, get metrics from Kruskal rows
        sig_metrics <- filtered_results %>% 
            filter(analysis_level == "kruskal") %>%
            filter(str_detect(metric, "^Number_of")) %>%
            pull(metric) %>%
            unique()
    } else {
        # For Wilcoxon
        sig_metrics <- filtered_results %>% 
            filter(str_detect(metric, "^Number_of")) %>%
            pull(metric) %>%
            unique()
    }
    
    if (length(sig_metrics) == 0) {
        message("No metrics found for ", test_type, " plotting in ", title_suffix)
        return(NULL)
    }
    
    # Create appropriate plots
    if (test_type == "wilcoxon") {
        return(create_wilcoxon_plots_cowplot(plot_data, filtered_results, group_col,
                                            sig_metrics, title_suffix, p_value_threshold,
                                            p_value_col, color_palette, base_size))
    } else {
        return(create_kruskal_plots_cowplot(plot_data, filtered_results, group_col,
                                           sig_metrics, title_suffix, p_value_threshold,
                                           p_value_col, color_palette, base_size,
                                           metrics_per_page))
    }
}

# Wilcoxon plotting function with cowplot styling
create_wilcoxon_plots_cowplot <- function(plot_data, results_table, group_col,
                                         sig_metrics, title_suffix, p_value_threshold,
                                         p_value_col, color_palette, base_size) {
    
    # Get group names and colors
    groups <- unique(na.omit(plot_data[[group_col]]))
    group_colors <- get_cowplot_colors(length(groups), color_palette)
    
    # Prepare data for plotting
    plot_data_long <- plot_data %>%
        select(all_of(c("sample_id", group_col)), any_of(sig_metrics)) %>%
        pivot_longer(
            cols = -c(sample_id, all_of(group_col)),
            names_to = "metric",
            values_to = "value"
        ) %>%
        mutate(
            metric_clean = gsub("^Number_of_", "", metric) %>%
                gsub("_Variants$", "", .) %>%
                gsub("_", " ", .)
        )
    
    # Get annotation data
    annotation_data <- results_table %>%
        filter(str_detect(metric, "^Number_of")) %>%
        mutate(
            metric_clean = gsub("^Number_of_", "", metric) %>%
                gsub("_Variants$", "", .) %>%
                gsub("_", " ", .),
            p_value_display = format.pval(.data[[p_value_col]], digits = 2)
        ) %>%
        select(metric, metric_clean, p_value_display) %>%
        distinct()
    
    # Create plots
    plots <- list()
    
    for (metric in sig_metrics) {
        metric_clean <- annotation_data$metric_clean[annotation_data$metric == metric][1]
        p_value_display <- annotation_data$p_value_display[annotation_data$metric == metric][1]
        
        # Filter data for this metric
        metric_data <- plot_data_long %>% filter(metric == !!metric)
        
        # Get groups present in data
        groups_present <- unique(metric_data[[group_col]])
        groups_present <- groups_present[!is.na(groups_present)]
        
        if (length(groups_present) < 2) {
            next
        }
        
        # Create Wilcoxon plot
        p <- ggplot(metric_data, aes(x = .data[[group_col]], y = value)) +
            geom_boxplot(aes(fill = .data[[group_col]]), 
                        outlier.shape = NA, width = 0.6, alpha = 0.7) +
            geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "grey30") +
            scale_fill_manual(values = group_colors) +
            labs(
                title = metric_clean,
                subtitle = paste("Wilcoxon p =", p_value_display),
                x = "Group",
                y = "Count"
            ) +
            theme_cowplot(font_size = base_size) +
            theme(
                plot.title = element_text(face = "bold"),
                plot.subtitle = element_text(color = "grey40"),
                legend.position = "none",
                panel.grid.major.y = element_line(color = "grey90", size = 0.2)
            )
        
        plots[[metric]] <- p
    }
    
    # Combine plots
    if (length(plots) == 0) {
        return(NULL)
    } else if (length(plots) == 1) {
        return(plots[[1]])
    } else {
        # Arrange in grid with cowplot
        plot_grid <- plot_grid(plotlist = plots, ncol = 2, 
                              labels = "AUTO", label_size = base_size)
        
        # Add title
        title_gg <- ggdraw() + 
            draw_label(paste("Wilcoxon Test Results:", title_suffix),
                      fontface = 'bold', size = base_size + 2)
        
        final_plot <- plot_grid(title_gg, plot_grid, 
                               ncol = 1, rel_heights = c(0.1, 1))
        
        return(final_plot)
    }
}

# Kruskal-Wallis plotting function with cowplot styling
# Kruskal-Wallis plotting function with cowplot styling - FIXED DATA HANDLING
create_kruskal_plots_cowplot <- function(plot_data, results_table, group_col,
                                        sig_metrics, title_suffix, p_value_threshold,
                                        p_value_col, color_palette, base_size,
                                        metrics_per_page = 4) {
    
    # Check if ggsignif is available
    if (!requireNamespace("ggsignif", quietly = TRUE)) {
        stop("Package 'ggsignif' is required for significance brackets. Please install it with: install.packages('ggsignif')")
    }
    
    # Get group names and colors
    groups <- unique(na.omit(plot_data[[group_col]]))
    group_colors <- get_cowplot_colors(length(groups), color_palette)
    
    # Create individual plots
    individual_plots <- list()
    
    for (metric in sig_metrics) {
        # Extract metric data
        metric_data <- plot_data %>%
            select(all_of(c("sample_id", group_col)), any_of(metric)) %>%
            rename(value = all_of(metric))
        
        # Get Kruskal stats for this metric
        kruskal_stats <- results_table %>% 
            filter(metric == !!metric & analysis_level == "kruskal")
        
        if (nrow(kruskal_stats) == 0) {
            cat("No Kruskal stats found for metric:", metric, "\n")
            next
        }
        
        # Get Dunn's test results for this metric
        dunn_stats <- results_table %>% 
            filter(metric == !!metric & analysis_level == "dunn")
        
        # Clean metric name for display
        metric_display <- gsub("^Number_of_", "", metric) %>%
            gsub("_Variants$", "", .) %>%
            gsub("_", " ", .)
        
        # Create base plot
        p <- ggplot(metric_data, aes(x = .data[[group_col]], y = value)) +
            geom_boxplot(aes(fill = .data[[group_col]]), 
                        outlier.shape = NA, alpha = 0.7) +
            geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.5, color = "grey30") +
            scale_fill_manual(values = group_colors) +
            labs(
                title = metric_display,
                x = "Group",
                y = "Count"
            ) +
            theme_cowplot(font_size = base_size) +
            theme(
                plot.title = element_text(face = "bold"),
                legend.position = "none",
                panel.grid.major.y = element_line(color = "grey90", size = 0.2)
            )
        
        # Add Kruskal-Wallis annotation
        if (!is.na(kruskal_stats[[p_value_col]][1])) {
            p_value_display <- format.pval(kruskal_stats[[p_value_col]][1], digits = 2)
            
            # Position the annotation
            max_value <- max(metric_data$value, na.rm = TRUE)
            annotation_y <- max_value * 1.1
            
            # Color based on significance
            annotation_color <- ifelse(kruskal_stats[[p_value_col]][1] < p_value_threshold, 
                                     "red3", "grey40")
            
            p <- p + annotate("text", x = 1, y = annotation_y,
                             label = paste("K-W p =", p_value_display),
                             hjust = 0, size = base_size / 3, 
                             color = annotation_color, fontface = "italic")
            
            # Add Dunn's test results if Kruskal is significant
            if (kruskal_stats[[p_value_col]][1] < p_value_threshold && nrow(dunn_stats) > 0) {
                # Prepare comparisons
                comparisons <- lapply(dunn_stats$comparison, function(x) {
                    strsplit(x, " - ")[[1]]
                })
                
                # Determine significance stars
                stars <- case_when(
                    dunn_stats$dunn_p < 0.001 ~ "***",
                    dunn_stats$dunn_p < 0.01 ~ "**",
                    dunn_stats$dunn_p < 0.05 ~ "*",
                    TRUE ~ "ns"
                )
                
                # Calculate y positions
                min_value <- min(metric_data$value, na.rm = TRUE)
                range_value <- max_value - min_value
                start_y <- max_value + (range_value * 0.12)
                step_y <- range_value * 0.15
                
                # Add significance brackets
                for (i in 1:length(comparisons)) {
                    if (length(comparisons[[i]]) == 2) {
                        group1 <- comparisons[[i]][1]
                        group2 <- comparisons[[i]][2]
                        
                        if (group1 %in% groups && group2 %in% groups) {
                            y_pos <- start_y + (i-1) * step_y
                            
                            cat("  Adding Dunn bracket for", metric, ":", 
                                group1, "vs", group2, "p =", 
                                format.pval(dunn_stats$dunn_p[i], digits = 2), "\n")
                            
                            p <- p + ggsignif::geom_signif(
                                comparisons = list(c(group1, group2)),
                                annotations = stars[i],
                                y_position = y_pos,
                                tip_length = 0.02,
                                textsize = base_size / 3.5,
                                vjust = -0.2,
                                size = 0.8,
                                color = "black",
                                map_signif_level = FALSE
                            )
                        }
                    }
                }
            }
        }
        
        individual_plots[[metric]] <- p
    }
    
    if (length(individual_plots) == 0) {
        return(NULL)
    }
    
    # Create paginated grids using cowplot
    plot_chunks <- split(
        individual_plots,
        ceiling(seq_along(individual_plots) / metrics_per_page)
    )
    
    grid_list <- list()
    
    for (page_num in seq_along(plot_chunks)) {
        # Create plot grid for this page
        plot_grid <- plot_grid(plotlist = plot_chunks[[page_num]], ncol = 2,
                              labels = "AUTO", label_size = base_size)
        
        # Add title
        title_gg <- ggdraw() + 
            draw_label(paste("Kruskal-Wallis Results:", title_suffix, "- Page", page_num),
                      fontface = 'bold', size = base_size + 2)
        
        final_plot <- plot_grid(title_gg, plot_grid, 
                               ncol = 1, rel_heights = c(0.1, 1))
        
        grid_list[[paste("Page", page_num)]] <- final_plot
    }
    
    return(grid_list)
}

# Save function for variant plots
save_variant_plots <- function(plots, prefix = "variant_plots", 
                              output_dir = ".",
                              output_width = 14, output_height = 10,
                              dpi = 300) {
    
    if (is.null(plots)) {
        message("No plots to save")
        return()
    }
    
    # Create output directory
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Handle different plot structures
    if (inherits(plots, "ggplot") || inherits(plots, "gg")) {
        # Single ggplot
        filename <- file.path(output_dir, paste0(prefix, ".png"))
        ggsave(filename, plots, width = output_width, height = output_height, dpi = dpi)
        message("Saved: ", filename)
        
    } else if (inherits(plots, "plot_grid") || inherits(plots, "gtable") || inherits(plots, "gTree")) {
        # Cowplot grid or other grid object
        filename <- file.path(output_dir, paste0(prefix, ".png"))
        png(filename, width = output_width, height = output_height, 
            units = "in", res = dpi)
        print(plots)
        dev.off()
        message("Saved: ", filename)
        
    } else if (is.list(plots)) {
        # List of plots
        for (plot_name in names(plots)) {
            plot_obj <- plots[[plot_name]]
            
            if (inherits(plot_obj, "ggplot") || inherits(plot_obj, "gg") || 
                inherits(plot_obj, "plot_grid") || inherits(plot_obj, "gtable")) {
                # Plot object
                filename <- file.path(output_dir, 
                                     paste0(prefix, "_", gsub(" ", "_", tolower(plot_name)), ".png"))
                
                if (inherits(plot_obj, "ggplot") || inherits(plot_obj, "gg")) {
                    ggsave(filename, plot_obj, width = output_width, height = output_height, dpi = dpi)
                } else {
                    png(filename, width = output_width, height = output_height, 
                        units = "in", res = dpi)
                    print(plot_obj)
                    dev.off()
                }
                message("Saved: ", filename)
                
            } else if (is.list(plot_obj)) {
                # Nested list (multiple datasets or pages)
                save_variant_plots(plot_obj, 
                                 paste0(prefix, "_", gsub(" ", "_", tolower(plot_name))),
                                 output_dir,
                                 output_width, output_height, dpi)
            }
        }
    }
}