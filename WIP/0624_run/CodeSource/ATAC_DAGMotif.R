library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots

combined_object <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/filtered_50_cells_atac_cleaned.qs")

# check
summary(combined_object@assays$ATAC@meta.features$GC.percent)
head(combined_object[["ATAC"]]@motifs)
dim(combined_object[['ATAC']][[]])
head(combined_object[['ATAC']][[]])
table(seqnames(Annotation(combined_object)))
seqlevels(Annotation(combined_object))
seqinfo(Annotation(combined_object))



# Initialize summary statistics
summary_stats <- data.frame()

# Directory setup
output_dir_base <-"~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
dir.create(output_dir_base)
dir.create(file.path(output_dir_base,  "DAG", "DA_ATAC"))
dir.create(file.path(output_dir_base, "DAG", "Sig_ATAC"))
dir.create(file.path(output_dir_base, "DAG", "Motif_Enrichment"), showWarnings = FALSE)


# cell type and comparison group list
cell_types <- unique(combined_object$cell_type)

comparisons <- list(
    PTZvsSAL_24hr = list(name = "PTZvsSAL_24hr", group1 = "PTZ_24hr", group2 = "SAL_24hr"),
    PTZvsSAL_1hr = list(name = "PTZvsSAL_1hr", group1 = "PTZ_1hr", group2 = "SAL_1hr"),
    `24hrvs1hr_PTZ` = list(name = "24hrvs1hr_PTZ", group1 = "PTZ_24hr", group2 = "PTZ_1hr"),
    `24hrvs1hr_SAL` = list(name = "24hrvs1hr_SAL", group1 = "SAL_24hr", group2 = "SAL_1hr")
)

# Motif anlaysis loop
for (ct in cell_types) {
    for (comp in comparisons) {
        tryCatch({
            message(sprintf("\nAnalyzing %s - %s", ct, comp$name))

            # subsetting cell type in each condition
            cell_type_data <- subset(combined_object,
                                     subset = cell_type == ct & Sample %in% c(comp$group1, comp$group2)) # cell_type is metadata column name for cell type, Sample is for group

            # cell_type_data <- PrepSCTFindMarkers(cell_type_data)

            # Set identities
            Idents(cell_type_data) <- "Sample"
            DefaultAssay(combined_object) <- "ATAC"

            # Perform DA peak analysis
            da_results <- FindMarkers(
                object = cell_type_data,
                # assay = "ATAC",
                # subset.ident = cell_type,
                group.by = "Sample",
                ident.1 = comp$group1,
                ident.2 = comp$group2,
                min.pct = 0.05,
                test.use = "LR",
                latent.vars = "nCount_ATAC"
            )

            # Add metadata
            da_results$gene <- rownames(da_results)
            da_results$cell_type <- ct
            da_results$comparison <- comp$name
            da_results$closestgene <- ClosestFeature(cell_type_data, regions =rownames(da_results))$gene_name
            da_results$distance <- ClosestFeature(cell_type_data, regions = rownames(da_results))$distance
            da_results <- da_results %>%
                mutate(Regulation = case_when(
                    avg_log2FC > 0 ~ "Up",
                    avg_log2FC < 0 ~ "Down",
                    TRUE ~ "NoChange"))

            # Filter significant DA peaks
            sig_atac_peaks <- da_results %>%
                mutate(Regulation = case_when(
                    avg_log2FC > 0 ~ "Up",
                    avg_log2FC < 0 ~ "Down",
                    TRUE ~ "NoChange")) %>%
                dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.2) %>%
                arrange(p_val_adj)


            # Save DA results
            safe_comp_name <- gsub("[:]", "_", comp$name)
            csv_filename <- sprintf("DA_results_%s_%s.csv", ct, safe_comp_name)
            write.csv(
                da_results,
                file = file.path(output_dir_base, "DAG", "DA_ATAC", csv_filename),
                row.names = FALSE
            )

            csv_filename2 <- sprintf("Sig_results_%s_%s.csv", ct, safe_comp_name)
            write.csv(
                sig_atac_peaks,
                file = file.path(output_dir_base, "DAG", "Sig_ATAC", csv_filename2),
                row.names = FALSE
            )

            message(sprintf("Analysis complete. Found %d significant DA peaks", length(sig_atac_peaks)))

            # **Perform Motif Analysis on Significant DA Peaks**

            temp_subset <- subset(
                combined_object,
                subset = cell_type == !!ct & Sample %in% c(comp$group1, comp$group2)
            )
            if (length(sig_atac_peaks) > 0) {
                # Split DA peaks into Upregulated & Downregulated
                up_peaks <- da_results$gene[da_results$p_val_adj < 0.05 & da_results$avg_log2FC > 0.2]
                down_peaks <- da_results$gene[da_results$p_val_adj < 0.05 & da_results$avg_log2FC < -0.2]

                # **Motif analysis for upregulated DA peaks**
                if (length(up_peaks) > 10) {
                    motif_results_up <- FindMotifs(
                        object = temp_subset,
                        assay = "ATAC",
                        features = up_peaks,
                        background = 500
                    )
                    motif_results_up$direction <- "up"
                    write.csv(
                        motif_results_up,
                        file = file.path(output_dir_base,"DAG", "Motif_Enrichment",
                                         sprintf("Motif_results_%s_%s_UP.csv", ct, safe_comp_name)),
                        row.names = FALSE
                    )
                } else {
                    message(sprintf("Skipping upregulated motif analysis: too few DA peaks (%d found)", length(up_peaks)))
                }

                # **Motif analysis for downregulated DA peaks**
                if (length(down_peaks) > 10) {
                    motif_results_down <- FindMotifs(
                        object = temp_subset,
                        assay = "ATAC",
                        features = down_peaks,
                        background = 500
                    )
                    motif_results_down$direction <- "down"
                    write.csv(
                        motif_results_down,
                        file = file.path(output_dir_base,"DAG", "Motif_Enrichment",
                                         sprintf("Motif_results_%s_%s_DOWN.csv", ct, safe_comp_name)),
                        row.names = FALSE
                    )
                } else {
                    message(sprintf("Skipping downregulated motif analysis: too few DA peaks (%d found)", length(down_peaks)))
                }

            } else {
                message("Skipping motif analysis: No significant DA peaks found.")
            }


        }, error = function(e) {
            message(sprintf("Error in analysis: %s", e$message))
        })
    }
}

# file read
csv_files <- list.files(file.path(output_dir_base,"DAG", "DA_ATAC"), pattern = "DA_results_.+\\.csv$", full.names = TRUE)

# loop
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DA_results_|.csv", "", file_name), "_")[[1]]
    cell_type <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Creating volcano plot for", cell_type, "-", comparison_name))

    # Read CSV file
    da_results <- read.csv(csv_file)

    # Create volcano plot
    volcano_plot <- EnhancedVolcano(da_results,
                                    lab = da_results$gene,
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    title = paste('Volcano plot -', cell_type),
                                    subtitle = safe_comparison_name,
                                    pCutoff = 0.05,
                                    FCcutoff = 0.2,
                                    pointSize = 1.0,
                                    labSize = 3.0
    )

    # Save volcano plot with safe filename
    plot_file_name <- paste0("Volcano_plot_", cell_type, "_", safe_comparison_name, ".png")
    ggsave(file.path(output_dir_base, "DAG", "Plots", "VolcanoPlots_ATAC", plot_file_name), volcano_plot, width = 12, height = 10)
    print(paste("Saved volcano plot:", plot_file_name))
}

