library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots

deg_dir <-  "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/sig_csv/"
da_dir <-  "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DAG/Sig_ATAC/"
output_dir <-  "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DAG/Plots/CoveragePlots_test/"
dir.create(output_dir)
combined_object <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/filtered_50_cells_atac_cleaned.qs")

gene_list <-  function(deg_dir, da_dir, top_n = 20, max_distance = 50000) {

    # Read all DEG and DA results
    deg_files <- list.files(deg_dir, pattern = ".csv$", full.names = TRUE)
    da_files <- list.files(da_dir, pattern = ".csv$", full.names = TRUE)[-43]

    # Priority gene lists
    gene_list <- c("Fos", "Jun", "Junb", "Fosb", "Egr1", "Egr2", "Arc", "Npas4",
                   "Nr4a1", "Bdnf", "Homer1", "Creb1", "Camk2a", "Grin1", "Dlg4")

    gene_list2 <- c("Dcx", "Neurod1", "Neurod6", "Prox1", "Tbr1", "Tbr2", "Satb2",
                    "Bcl11b", "Sox2", "Pax6")

    # Combine all significant genes from DEG analysis
    all_sig_genes <- c()
    deg_summary <- data.frame()

    for (file in deg_files) {
        deg_data <- read.csv(file)
        if ("gene" %in% colnames(deg_data) && nrow(deg_data) > 0) {
            all_sig_genes <- c(all_sig_genes, deg_data$gene)

            # Add metadata for tracking
            file_name <- basename(file)
            parts <- strsplit(gsub(".csv", "", file_name), "_")[[1]]
            cell_type <- parts[1]
            comparison <- paste(parts[-1], collapse = "_")

            deg_data$source_file <- file_name
            deg_data$cell_type <- cell_type
            deg_data$comparison <- comparison

            deg_summary <- rbind(deg_summary, deg_data)
        }
    }

    # Process DA results and extract genes near DA peaks
    all_da_genes <- c()
    da_summary <- data.frame()

    for (file in da_files) {
        da_data <- read.csv(file)

        # Check if the file has the expected columns
        if ("closestgene" %in% colnames(da_data) && nrow(da_data) > 0) {

            if (nrow(da_data) > 0) {
                # Extract genes near DA peaks
                all_da_genes <- c(all_da_genes, da_data$closestgene)

                # Add metadata for tracking
                file_name <- basename(file)
                parts <- strsplit(gsub("DA_results_|.csv", "", file_name), "_")[[1]]
                cell_type <- parts[1]
                comparison <- paste(parts[-1], collapse = "_")

                da_data$source_file <- file_name
                da_data$cell_type_da <- cell_type
                da_data$comparison_da <- comparison

                da_summary <- rbind(da_summary, da_data)
            }
        }
    }

    # Remove NA and empty gene names
    all_da_genes <- all_da_genes[!is.na(all_da_genes) & all_da_genes != ""]

    # Find overlaps (genes that have both expression and accessibility changes)
    overlapping_genes <- intersect(all_sig_genes, all_da_genes)
    priority_overlap <- intersect(overlapping_genes, c(gene_list, gene_list2))

    # Score genes by frequency and priority
    all_genes <- unique(c(all_sig_genes, all_da_genes))

    gene_scores <- data.frame(
        gene = all_genes,
        deg_count = sapply(all_genes, function(x) sum(all_sig_genes == x)),
        da_count = sapply(all_genes, function(x) sum(all_da_genes == x)),
        priority = sapply(all_genes, function(x) x %in% c(gene_list, gene_list2)),
        overlap = sapply(all_genes, function(x) x %in% overlapping_genes)
    )

    # Create detailed mapping for overlapping genes
    overlap_details <- data.frame()

    for (gene in overlapping_genes) {
        # Get DEG info
        deg_info <- deg_summary %>% filter(gene == !!gene)

        # Get DA info
        da_info <- da_summary %>% filter(closestgene == !!gene)

        # Combine information
        if (nrow(deg_info) > 0 && nrow(da_info) > 0) {
            for (i in 1:nrow(deg_info)) {
                for (j in 1:nrow(da_info)) {
                    combined_info <- data.frame(
                        gene = gene,
                        deg_cell_type = deg_info$cell_type[i],
                        deg_comparison = deg_info$comparison[i],
                        deg_log2FC = deg_info$avg_log2FC[i],
                        deg_padj = deg_info$p_val_adj[i],
                        da_cell_type = da_info$cell_type_da[j],
                        da_comparison = da_info$comparison_da[j],
                        da_log2FC = da_info$avg_log2FC[j],
                        da_padj = da_info$p_val_adj[j],
                        da_peak = da_info$ranges[j],
                        distance_to_gene = da_info$distance[j],
                        regulation_deg = ifelse("Regulation" %in% colnames(deg_info), deg_info$Regulation[i],
                                                ifelse(deg_info$avg_log2FC[i] > 0, "Up", "Down")),
                        regulation_da = da_info$Regulation[j],
                        matched_analysis = (deg_info$cell_type[i] == da_info$cell_type_da[j] &
                                                deg_info$comparison[i] == da_info$comparison_da[j])
                    )
                    overlap_details <- rbind(overlap_details, combined_info)
                }
            }
        }
    }
    # save
    saveRDS(gene_scores, file = file.path(output_dir, "gene_scores.rds"))
    saveRDS(priority_overlap, file = file.path(output_dir, "priority_overlap.rds"))
    saveRDS(overlap_details, file = file.path(output_dir, "overlap_details.rds"))
    saveRDS(overlapping_genes, file = file.path(output_dir, "overlapping_genes.rds"))
    saveRDS(deg_summary, file = file.path(output_dir, "deg_summary.rds"))
    saveRDS(da_summary, file = file.path(output_dir, "da_dummary.rds"))

    return(list(
        top_genes = head(gene_scores, top_n),
        priority_overlap = priority_overlap,
        overlapping_genes = overlapping_genes,
        overlap_details = overlap_details,
        deg_summary = deg_summary,
        da_summary = da_summary,
        recommendations = list(
            immediate_early = intersect(gene_list, gene_scores$gene[1:50]),
            developmental = intersect(gene_list2, gene_scores$gene[1:50]),
            high_confidence = gene_scores$gene[gene_scores$overlap == TRUE][1:10],
            top_scored = head(gene_scores$gene, 10),
            perfect_matches = overlap_details %>%
                filter(matched_analysis == TRUE) %>%
                arrange(deg_padj, da_padj) %>%
                pull(gene) %>%
                unique() %>%
                head(10)
        )
    ))
}

# Function to create comprehensive coverage plots
coverage_plot <- function(analysis_results, combined_merged, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # 1. PERFECT MATCHES (same cell type and comparison for both DEG and DA)
    cat("Creating coverage plots for perfect matches (DEG + DA in same analysis)...\n")

    perfect_matches <- analysis_results$overlap_details %>%
        filter(matched_analysis == TRUE) %>%
        arrange(deg_padj, da_padj) %>%
        head(10)

    for (i in 1:nrow(perfect_matches)) {
        match_info <- perfect_matches[i, ]

        if (match_info$deg_cell_type %in% combined_merged$cell_type) {

            cell_subset <- subset(combined_merged, subset = cell_type == match_info$deg_cell_type)
            DefaultAssay(cell_subset) <- "ATAC"

            # Define groups
            groups <- switch(match_info$deg_comparison,
                             "PTZvsSAL_1hr" = c("SAL_1hr", "PTZ_1hr"),
                             "PTZvsSAL_24hr" = c("SAL_24hr", "PTZ_24hr"),
                             "24hrvs1hr_PTZ" = c("PTZ_1hr", "PTZ_24hr"),
                             "24hrvs1hr_SAL" = c("SAL_1hr", "SAL_24hr")
            )

            if (match_info$gene %in% rownames(cell_subset)) {
                tryCatch({

                    # Gene-centered view (broader regulatory landscape)
                    gene_plot <- CoveragePlot(
                        object = cell_subset,
                        region = match_info$gene,
                        features = match_info$gene,
                        idents = groups,
                        peaks.group.by = "Sample",
                        expression.assay = "RNA",
                        extend.upstream = 100000,
                        extend.downstream = 50000
                    ) + ggtitle(paste("Gene View:", match_info$gene,
                                      "| DEG log2FC:", round(match_info$deg_log2FC, 2),
                                      "| DA log2FC:", round(match_info$da_log2FC, 2)))

                    # Peak-centered view (focused on specific DA region)
                    peak_plot <- CoveragePlot(
                        object = cell_subset,
                        region = match_info$da_peak,
                        idents = groups,
                        peaks.group.by = "Sample",
                        extend.upstream = 25000,
                        extend.downstream = 25000
                    ) + ggtitle(paste("Peak View:", match_info$da_peak,
                                      "| Distance to", match_info$gene, ":", match_info$distance_to_gene, "bp"))

                    # Save both plots
                    gene_filename <- paste0("Coverage_GENE_", match_info$gene, "_",
                                            match_info$deg_cell_type, "_", match_info$deg_comparison, ".png")
                    peak_filename <- paste0("Coverage_PEAK_", gsub("[-:]", "_", match_info$da_peak), "_",
                                            match_info$deg_cell_type, "_", match_info$deg_comparison, ".png")

                    ggsave(file.path(output_dir, gene_filename), gene_plot, width = 16, height = 10, dpi = 300)
                    ggsave(file.path(output_dir, peak_filename), peak_plot, width = 14, height = 8, dpi = 300)

                    cat(sprintf("✓ Created plots for %s (DEG: %s, DA: %s)\n",
                                match_info$gene, match_info$regulation_deg, match_info$regulation_da))

                }, error = function(e) {
                    cat(sprintf("✗ Error plotting %s: %s\n", match_info$gene, e$message))
                })
            }
        }
    }

    # 2. PRIORITY GENES (regardless of perfect matching)
    cat("\nCreating coverage plots for priority genes...\n")

    priority_genes <- c("Fos", "Jun", "Egr1", "Arc", "Npas4", "Bdnf", "Dcx", "Prox1")
    available_priority <- intersect(priority_genes, analysis_results$overlapping_genes)

    for (gene in available_priority) {
        # Get all DEG instances of this gene
        gene_deg_info <- analysis_results$deg_summary %>% filter(gene == !!gene)

        for (i in 1:nrow(gene_deg_info)) {
            deg_info <- gene_deg_info[i, ]

            if (deg_info$cell_type %in% combined_merged$cell_type) {

                cell_subset <- subset(combined_merged, subset = cell_type == deg_info$cell_type)
                DefaultAssay(cell_subset) <- "ATAC"

                groups <- switch(deg_info$comparison,
                                 "PTZvsSAL_1hr" = c("SAL_1hr", "PTZ_1hr"),
                                 "PTZvsSAL_24hr" = c("SAL_24hr", "PTZ_24hr"),
                                 "24hrvs1hr_PTZ" = c("PTZ_1hr", "PTZ_24hr"),
                                 "24hrvs1hr_SAL" = c("SAL_1hr", "SAL_24hr")
                )

                if (gene %in% rownames(cell_subset)) {
                    tryCatch({
                        plot <- CoveragePlot(
                            object = cell_subset,
                            region = gene,
                            features = gene,
                            idents = groups,
                            peaks.group.by = "Sample",
                            expression.assay = "RNA",
                            extend.upstream = 75000,
                            extend.downstream = 50000
                        ) + ggtitle(paste("Priority Gene:", gene, "| Cell Type:", deg_info$cell_type,
                                          "| Comparison:", deg_info$comparison))

                        filename <- paste0("Coverage_PRIORITY_", gene, "_", deg_info$cell_type, "_", deg_info$comparison, ".png")
                        ggsave(file.path(output_dir, filename), plot, width = 16, height = 10, dpi = 300)

                        cat(sprintf("✓ Created priority plot for %s\n", gene))

                    }, error = function(e) {
                        cat(sprintf("✗ Error plotting priority gene %s: %s\n", gene, e$message))
                    })
                }
            }
        }
    }

    # 3. TOP DA PEAKS (highest significance)
    cat("\nCreating coverage plots for top DA peaks...\n")

    top_da_peaks <- analysis_results$da_summary %>%
        arrange(p_val_adj) %>%
        head(10)

    for (i in 1:nrow(top_da_peaks)) {
        peak_info <- top_da_peaks[i, ]

        if (peak_info$cell_type_da %in% combined_merged$cell_type) {

            cell_subset <- subset(combined_merged, subset = cell_type == peak_info$cell_type_da)
            DefaultAssay(cell_subset) <- "ATAC"

            groups <- switch(peak_info$comparison_da,
                             "PTZvsSAL_1hr" = c("SAL_1hr", "PTZ_1hr"),
                             "PTZvsSAL_24hr" = c("SAL_24hr", "PTZ_24hr"),
                             "24hrvs1hr_PTZ" = c("PTZ_1hr", "PTZ_24hr"),
                             "24hrvs1hr_SAL" = c("SAL_1hr", "SAL_24hr")
            )

            tryCatch({
                plot <- CoveragePlot(
                    object = cell_subset,
                    region = peak_info$ranges,
                    idents = groups,
                    peaks.group.by = "Sample",
                    extend.upstream = 20000,
                    extend.downstream = 20000
                ) + ggtitle(paste("Top DA Peak near", peak_info$closestgene,
                                  "| Distance:", peak_info$distance, "bp",
                                  "| log2FC:", round(peak_info$avg_log2FC, 2)))

                safe_peak_name <- gsub("[-:]", "_", peak_info$ranges)
                filename <- paste0("Coverage_TOP_DA_", safe_peak_name, "_", peak_info$cell_type_da, "_", peak_info$comparison_da, ".png")
                ggsave(file.path(output_dir, filename), plot, width = 14, height = 8, dpi = 300)

                cat(sprintf("✓ Created top DA peak plot near %s\n", peak_info$closestgene))

            }, error = function(e) {
                cat(sprintf("✗ Error plotting DA peak near %s: %s\n", peak_info$closestgene, e$message))
            })
        }
    }
}

# Function to create summary tables and plots
analysis_summary <- function(analysis_results, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # 1. Summary statistics
    summary_stats <- data.frame(
        metric = c("Total DEG genes", "Total genes near DA peaks", "Overlapping genes",
                   "Priority overlapping genes", "Perfect matches"),
        count = c(
            length(unique(analysis_results$deg_summary$gene)),
            length(unique(analysis_results$da_summary$closestgene)),
            length(analysis_results$overlapping_genes),
            length(analysis_results$priority_overlap),
            sum(analysis_results$overlap_details$matched_analysis)
        )
    )

    write.csv(summary_stats, file.path(output_dir, "analysis_summary_stats.csv"), row.names = FALSE)

    # 2. Detailed overlap information
    write.csv(analysis_results$overlap_details, file.path(output_dir, "overlap_details.csv"), row.names = FALSE)

    # 3. Top gene scores
    write.csv(analysis_results$top_genes, file.path(output_dir, "top_gene_scores.csv"), row.names = FALSE)

    # 4. Recommendations
    writeLines(
        paste0("Coverage Plot Recommendations:\n",
               "Immediate Early Genes: ", paste(analysis_results$recommendations$immediate_early, collapse = ", "), "\n",
               "Developmental Genes: ", paste(analysis_results$recommendations$developmental, collapse = ", "), "\n",
               "High Confidence: ", paste(analysis_results$recommendations$high_confidence, collapse = ", "), "\n",
               "Perfect Matches: ", paste(analysis_results$recommendations$perfect_matches, collapse = ", ")),
        file.path(output_dir, "recommendations.txt")
    )

    cat("Summary files saved to:", output_dir, "\n")
}


updated_results <- gene_list(
    deg_dir = deg_dir,
    da_dir = da_dir,
    top_n = 20,
    max_distance = 50000
)
saveRDS(updated_results, file = file.path(output_dir, "updated_results.rds"))

# Create coverage plots and summaries
coverage_plot(
    analysis_results = updated_results,
    combined_merged = combined_object,
    output_dir = output_dir
)

analysis_summary(
    analysis_results = updated_results,
    output_dir = output_dir
)
