create_deg_summary <- function(deg_dir) {
    # Get all DEG files
    deg_files <- list.files(deg_dir, pattern = ".csv$", full.names = TRUE)

    # Initialize list to store summaries
    deg_summaries <- list()

    for (file in deg_files) {

        # Extract metadata from filename
        file_name <- basename(file)
        parts <- strsplit(gsub(".csv", "", file_name), "_")[[1]]
        cell_type <- parts[1]
        comparison_name <- paste(parts[-1], collapse = "_")

        # Read DEG results
        deg_data <- read.csv(file)

        # Calculate summary statistics
        total_genes <- nrow(deg_data)
        up_regulated <- sum(deg_data$avg_log2FC > 0.2 & deg_data$p_val_adj < 0.05)
        down_regulated <- sum(deg_data$avg_log2FC < -0.2 & deg_data$p_val_adj < 0.05)
        significant_genes <- sum(deg_data$p_val_adj < 0.05)

        # Compute top genes
        top_up_genes <- deg_data$gene[deg_data$avg_log2FC > 0.2 & deg_data$p_val_adj < 0.05]
        top_down_genes <- deg_data$gene[deg_data$avg_log2FC < -0.2 & deg_data$p_val_adj < 0.05]
        top_sig_genes <- deg_data$gene[deg_data$p_val_adj < 0.05]

        # Get significant up and down genes separately
        sig_up_genes <- deg_data$gene[deg_data$avg_log2FC > 0 & deg_data$p_val_adj < 0.05]
        sig_down_genes <- deg_data$gene[deg_data$avg_log2FC < 0 & deg_data$p_val_adj < 0.05]

        # Create summary row
        summary_row <- data.frame(
            Cell_Type = cell_type,
            Comparison = comparison_name,
            Total_DEGs = total_genes,
            Upregulated = up_regulated,
            Downregulated = down_regulated,
            Significant_Genes = significant_genes,
            Sig_Upregulated = length(sig_up_genes),
            Sig_Downregulated = length(sig_down_genes),
            Top_Up_Genes = paste(head(top_up_genes, 10), collapse = ", "),
            Top_Down_Genes = paste(head(top_down_genes, 10), collapse = ", "),
            Top_Significant_Genes = paste(head(top_sig_genes, 10), collapse = ", "),
            Sig_Up_Genes_List = paste(head(sig_up_genes, 10), collapse = ", "),
            Sig_Down_Genes_List = paste(head(sig_down_genes, 10), collapse = ", ")
        )

        deg_summaries[[length(deg_summaries) + 1]] <- summary_row
    }

    # Combine all summaries into a single dataframe
    deg_summary_df <- do.call(rbind, deg_summaries)

    return(deg_summary_df)
}

# Generate summary
complete_summary <- create_deg_summary(deg_dir = bg_dir)

# Save summary
write.csv(complete_summary, file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/DEG_Complete_Analysis_Summary.csv"), row.names = FALSE)

complete_summary <- complete_summary %>%
    dplyr::rename(Cell_Type = Cell_Type,
                  Comparison = Comparison,
                  Total_DEGs = Total_DEGs,
                  Upregulated = Upregulated,
                  Downregulated = Downregulated,
                  Significant_Genes = Significant_Genes)

# Pivot `complete_summary` for stacked bar chart
deg_long <- complete_summary %>%
    pivot_longer(cols = c(Upregulated, Downregulated), names_to = "Regulation", values_to = "Count")

# 1. Stacked Bar Chart of Up/Down Regulated Genes by Cell Type and Comparison
p1_stacked <- ggplot(deg_long, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Comparison, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(10, 20, 40, 20)) +
    labs(
        title = "Up/Down Regulated Genes by Cell Type and Comparison",
        subtitle = "Significant genes (p < 0.05, |log2FC| > 0.)",
        y = "Number of Genes"
    ) +
    scale_fill_manual(values = c("Upregulated" = "#66c2a5", "Downregulated" = "#fc8d62"), name = "Regulation") +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3)

# 2. Side-by-Side Bar Chart for Up/Down Regulated Genes
p1 <- ggplot(deg_long, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_wrap(~ Comparison, scales = "free_y") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 40, 20),
        panel.grid.major.x = element_blank()
    ) +
    labs(
        title = "Up/Down Regulated Genes by Cell Type and Comparison",
        subtitle = "Significant genes (p < 0.05, |log2FC| > 0.)",
        y = "Number of Genes"
    ) +
    scale_fill_manual(values = c("Upregulated" = "#66c2a5", "Downregulated" = "#fc8d62"), name = "Regulation") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3)

ggsave(file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/Sig_UpDown_Distribution_stacked.png"), p1_stacked, width = 18, height = 12)
ggsave(file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/Sig_UpDown_Distribution.png"), p1, width = 16, height = 12)

# 3. Comparison of Total DEGs vs Significant Genes
long_summary <- complete_summary %>%
    pivot_longer(cols = c(Total_DEGs, Significant_Genes), names_to = "Gene_Category", values_to = "Count") %>%
    mutate(Gene_Category = recode(Gene_Category,
                                  "Total_DEGs" = "Total Genes",
                                  "Significant_Genes" = "DEGs"))

p2 <- ggplot(long_summary, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Gene_Category)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ Comparison, scales = "free_y") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 40, 20),
        panel.grid.major.x = element_blank()
    ) +
    labs(
        title = "Comparison of Total Genes vs DEGs",
        subtitle = "Significant genes (p < 0.05)",
        y = "Number of Genes"
    ) +
    scale_fill_manual(values = c("Total Genes" = "#66c2a5", "DEGs" = "#fc8d62"), name = "Gene_Category") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5)

p2_stacked <- ggplot(long_summary, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Gene_Category)) +
    geom_bar(stat = "identity", position = "stack") +  # Stacked bars
    facet_wrap(~ Comparison, scales = "free_y") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 40, 20)
    ) +
    labs(
        title = "Comparison of Total Genes vs DEGs",
        subtitle = "Significant genes (p < 0.05)",
        y = "Number of Genes"
    ) +
    scale_fill_manual(values = c("Total Genes" = "#66c2a5", "DEGs" = "#fc8d62"), name = "Gene Category") +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3)


ggsave(file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/DEG_vs_SigGenes.png"), p2, width = 16, height = 12)
ggsave(file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/DEG_vs_SigGenes_stacked.png"), p2_stacked, width = 16, height = 12)


# 4. Create Individual Plots for Each Comparison
for (comp in unique(deg_long$Comparison)) {
    p4 <- deg_long %>% filter(Comparison == comp) %>%
        ggplot(aes(x = reorder(Cell_Type, -Count), y = Count, fill = Regulation)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(10, 20, 40, 20),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        ) +
        labs(
            title = paste("Sigs in", comp),
            subtitle = "",
            y = "Number of Genes",
            x = "Cell Type"
        ) +
        scale_fill_manual(values = c("Downregulated" = "#2166AC", "Upregulated" = "#B2182B")) +
        geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5)

    ggsave(file.path(paste0("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/Sig_UpDown_", comp, "_from_summary.png")), p4, width = 10, height = 6, dpi = 300)
}

for (comp in unique(long_summary$Comparison)) {
    p5 <- long_summary %>% filter(Comparison == comp) %>%
        ggplot(aes(x = reorder(Cell_Type, -Count), y = Count, fill = Gene_Category)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(10, 20, 40, 20),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        ) +
        labs(
            title = paste("DEG vs Sig in", comp),
            subtitle = "",
            y = "Number of Genes",
            x = "Cell Type"
        ) +
        scale_fill_manual(values = c("Total Genes" = "#2166AC", "DEGs" = "#B2182B")) +
        geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5)

    ggsave(file.path(paste0("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/DEG_vs_Sig_", comp, "_from_summary.png")), p5, width = 10, height = 6, dpi = 300)
}

# stacked individual
for (comp in unique(deg_long$Comparison)) {
    p4_stacked <- deg_long %>% filter(Comparison == comp) %>%
        ggplot(aes(x = reorder(Cell_Type, -Count), y = Count, fill = Regulation)) +
        geom_bar(stat = "identity", position = "stack") +  # Stacked bars
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(10, 20, 40, 20),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        ) +
        labs(
            title = paste("Sigs in", comp),
            subtitle = "",
            y = "Number of Genes",
            x = "Cell Type"
        ) +
        scale_fill_manual(values = c("Downregulated" = "#2166AC", "Upregulated" = "#B2182B")) +
        geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3.5)

    ggsave(file.path(paste0("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/Sig_UpDown_", comp, "_stacked.png")),
           p4_stacked, width = 10, height = 6, dpi = 300)
}

for (comp in unique(long_summary$Comparison)) {
    p5_stacked <- long_summary %>% filter(Comparison == comp) %>%
        ggplot(aes(x = reorder(Cell_Type, -Count), y = Count, fill = Gene_Category)) +
        geom_bar(stat = "identity", position = "stack") +  # Stacked bars
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major.x = element_blank(),
            plot.margin = margin(10, 20, 40, 20),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        ) +
        labs(
            title = paste("DEG vs Sig in", comp),
            subtitle = "",
            y = "Number of Genes",
            x = "Cell Type"
        ) +
        scale_fill_manual(values = c("Total Genes" = "#2166AC", "DEGs" = "#B2182B")) +
        geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3.5)

    ggsave(file.path(paste0("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/Summary/DEG_vs_Sig_", comp, "_stacked.png")),
           p5_stacked, width = 10, height = 6, dpi = 300)
}

# Loop through each comparison and create a scatter plot
comparisons <- unique(complete_summary$Comparison)

# Define a function to calculate Pearson correlation and linear regression for each comparison
calculate_correlation <- function(data) {
    correlation <- cor.test(data$Total_DEGs, data$Significant_Genes, method = "pearson")

    # Perform linear regression
    linear_model <- lm(Significant_Genes ~ Total_DEGs, data = data)
    summary_lm <- summary(linear_model)

    # Extract correlation coefficient and p-value
    correlation_coefficient <- correlation$estimate
    p_value <- correlation$p.value
    lm_p_value <-
        summary_lm$coefficients[2, 4]  # p-value of the slope

    return(list(correlation_coefficient = correlation_coefficient,
                pearson_p_value = p_value,
                lm_p_value = lm_p_value))
}

# Initialize a list to store results
correlation_results <- list()

# Loop through each comparison and calculate correlation
comparisons <- unique(complete_summary$Comparison)

for (comp in comparisons) {
    # Filter data for the current comparison
    comp_data <- complete_summary %>% filter(Comparison == comp)
    # Calculate correlation and linear regression for this comparison
    result <- calculate_correlation(comp_data)
    # Store results in a data frame for easy viewing
    correlation_results[[comp]] <- data.frame(Comparison = comp,
                                              Pearson_Correlation = result$correlation_coefficient,
                                              Pearson_p_value = result$pearson_p_value,
                                              LM_p_value = result$lm_p_value)

    # Create the scatter plot for the current comparison with updated style
    scatter_plot <- ggplot(comp_data, aes(x = Total_DEGs, y = Significant_Genes, color = Cell_Type)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_text_repel(aes(label = Cell_Type), size = 3, max.overlaps = 10) +
        geom_smooth(method = "lm", color = "red", se = FALSE) +
        labs(title = paste("Total Genes vs DEGs -", comp),
             subtitle = paste0("Pearson Correlation: ", round(result$correlation_coefficient, 2), ", p-value: ", signif(result$pearson_p_value, 3)),
             x = "Total Number of Genes",
             y = "Number of DEGs") +
        theme_minimal() +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 12),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))

    # Save the plot
    plot_file_name <- paste0("Scatter_Total_vs_Significant_DEGs_", gsub("[:]", "_", comp), ".png")
    ggsave(filename = file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Plots/ScatterPlots", plot_file_name), plot = scatter_plot, width = 10, height = 8, dpi = 300)
}

# Combine results into a single data frame
correlation_results_df <- do.call(rbind, correlation_results)

# Save the results as a CSV file
write.csv(correlation_results_df, file = file.path("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/Correlation_Results.csv"), row.names = FALSE)

output_dir_base <-"~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
sig_gene_files <- list.files(path = file.path(output_dir_base, "DEG", "sig_csv"), pattern = ".csv$", full.names = TRUE)

# combine list
deg_list <- lapply(
    sig_gene_files,
    function(file) {
        read_csv(file) %>%
            mutate(condition = paste(cell_type, comparison, sep = "_"))
    }
)

deg_all <- bind_rows(deg_list)

deg_top10 <-
    deg_all %>%
    group_by(condition) %>%
    arrange(p_val_adj, .by_group = TRUE) %>%
    slice_head(n=10) %>%
    dplyr::select(avg_log2FC, p_val_adj, gene, condition, cell_type, comparison) %>%
    mutate(comparison = factor(comparison, levels = c(
        "PTZvsSAL_1hr",
        "PTZvsSAL_24hr",
        "24hrvs1hr_PTZ",
        "24hrvs1hr_SAL"
    )))

# plot
plot <- ggplot(deg_top10, aes(x = comparison, y = avg_log2FC, group = gene, color = gene), label = T) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "Gene-wise avg_log2FC Across Timepoints",
         x = "Comparison Group",
         y = "avg_log2FC (log2)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ cell_type, scales = "free_y") +
    guides(color = "none")

ggsave(filename = file.path(output_dir_base, "DEG", "Plots", "Gene_Behavior", "LinePlot_all_by_pattern.png"), plot, width = 12, height =12, dpi = 300)

# more
deg_wide <- deg_all %>%
    # filter(comparison %in% c("PTZvsSAL_1hr", "PTZvsSAL_24hr", "24hrvs1hr_PTZ")) %>%
    dplyr::select(gene, cell_type, comparison, avg_log2FC) %>%
    pivot_wider(
        names_from = comparison,
        values_from = avg_log2FC
    )

deg_wide <- deg_wide %>%
    mutate(
        pattern = case_when(
            PTZvsSAL_1hr > 0 & PTZvsSAL_24hr > 0 ~ "Increase maintain",
            PTZvsSAL_1hr > 0 & PTZvsSAL_24hr < 0 ~ "Increase then Reversal",
            PTZvsSAL_1hr > 0 & is.na(PTZvsSAL_24hr) ~ "Increase then Baseline",
            PTZvsSAL_1hr < 0 & PTZvsSAL_24hr < 0 ~ "Sustained Decrease",
            PTZvsSAL_1hr < 0 & PTZvsSAL_24hr > 0 ~ "Decrease then Reversal",
            PTZvsSAL_1hr < 0 & is.na(PTZvsSAL_24hr) ~ "Decrease then Baseline",
            is.na(PTZvsSAL_1hr) & PTZvsSAL_24hr > 0 ~ "Late Upregulation",
            is.na(PTZvsSAL_1hr) & PTZvsSAL_24hr < 0 ~ "Late Downregulation",
            TRUE ~ "Other"
        ))

bar <- deg_wide %>%
    dplyr::count(pattern, sort = TRUE) %>%
    ggplot(aes(x = reorder(pattern, n), y = n, fill = pattern)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(title = "Gene Count per Pattern", x = "Pattern", y = "Gene Count")

ggsave(filename = file.path(output_dir_base, "DEG", "Plots", "Gene_Behavior", "BarPlot_all_by_pattern.png"), bar, width = 12, height =12, dpi = 300)

bar_cell <- deg_wide %>%
    dplyr::count(cell_type, pattern, sort = TRUE) %>%
    ggplot(aes(x = reorder(pattern, n), y = n, fill = pattern)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    theme_minimal() +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(title = "Gene Count per Pattern by Cell Type", x = "Pattern", y = "Gene Count")

ggsave(filename = file.path(output_dir_base, "DEG", "Plots", "Gene_Behavior", "BarPlot_all_by_pattern_cell.png"), bar_cell, width = 16, height =12, dpi = 300)

#  Injection point (baseline)
baseline_point <- deg_wide %>%
    distinct(gene, cell_type, pattern) %>%
    mutate(
        timepoint = "Injection",
        avg_log2FC = 0
    )

# 2. Long-format main timepoints
deg_long <- deg_wide %>%
    pivot_longer(
        cols = c("PTZvsSAL_1hr", "PTZvsSAL_24hr"),
        names_to = "comparison",
        values_to = "avg_log2FC"
    ) %>%
    mutate(
        timepoint = case_when(
            comparison == "PTZvsSAL_1hr" ~ "PTZ_1hr",
            comparison == "PTZvsSAL_24hr" ~ "PTZ_24hr"
        )
    ) %>%
    dplyr::select(gene, cell_type, pattern, timepoint, avg_log2FC)

deg_long_full <- bind_rows(baseline_point, deg_long) %>%
    arrange(gene, cell_type, timepoint) %>%
    mutate(group_id = paste(gene, cell_type, sep = "_"))

deg_long_full_filtered <- deg_long_full %>%
    filter(pattern != "Other")

cell <- ggplot(deg_long_full_filtered, aes(x = timepoint, y = avg_log2FC, group = group_id, color = pattern)) +
    geom_line(alpha = 0.7, linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ cell_type, scales = "free_y") +
    theme_minimal() +
    labs(
        title = "Gene Expression Over Time by Cell Type",
        x = "Timepoint",
        y = "avg_log2FC (PTZ vs SAL)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray")

ggsave(filename = file.path(output_dir_base, "DEG", "Plots", "Gene_Behavior", "LinePlot_cell.png"), cell, width = 16, height =12, dpi = 300)

pattern <- ggplot(deg_long_full_filtered, aes(x = timepoint, y = avg_log2FC, group = pattern)) +
    geom_line(aes(color = pattern), alpha = 0.5) +
    facet_grid(cell_type ~ pattern, scales = "free_y") +
    theme_minimal() +
    labs(title = "Gene Expression by Pattern and Cell Type",
         x = "Timepoint", y = "avg_log2FC")

ggsave(filename = file.path(output_dir_base, "DEG", "Plots", "Gene_Behavior", "LinePlot_grid.png"), pattern, width = 16, height =12, dpi = 300)

# simplified
pattern_summary_by_cell <- deg_long_full_filtered %>%
    group_by(cell_type, pattern, timepoint) %>%
    summarise(
        mean_log2FC = mean(avg_log2FC, na.rm = TRUE),
        gene_count = n(),
        .groups = "drop"
    )

pattern_summary_by_cell <- pattern_summary_by_cell %>%
    mutate(timepoint = factor(timepoint, levels = c("Injection", "PTZ_1hr", "PTZ_24hr")))

pattern_plot_by_cell <- ggplot(pattern_summary_by_cell, aes(x = timepoint, y = mean_log2FC, group = pattern, color = pattern)) +
    geom_line(linewidth = 1) +
    geom_point(aes(size = gene_count)) +
    scale_size_continuous(range = c(2, 8)) +
    facet_wrap(~ cell_type, scales = "free_y") +
    theme_minimal() +
    labs(
        title = "Average Gene Expression per Pattern Over Time by Cell Type",
        x = "Timepoint",
        y = "Mean avg_log2FC",
        size = "Gene Count",
        color = "Pattern"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )

ggsave(filename = file.path(output_dir_base, "DEG", "Plots", "Gene_Behavior", "LinePlot_agg.png"), pattern_plot_by_cell, width = 16, height =12, dpi = 300)
