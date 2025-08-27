library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots

combined_object <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/filtered_50_cells_atac_cleaned.qs")

DefaultAssay(combined_object) <- "ATAC"
 add_closest_gene <- function(file_list, output_dir) {
     for (file in file_list) {
         da_data <- read_csv(csv_files[10])
         file_name <- basename(csv_files[10])

         da_data$closestgene <- ClosestFeature(da_data, regions = da_data$gene)$gene_name
         da_data$distance <- ClosestFeature(da_data, regions = rownames(da_data))$distance

         write_csv(da_data, file = paste0(output_dir, file_name))
     }
 }

 # cell type and comparison group list
 cell_types <- unique(combined_object$cell_type)

 comparisons <- list(
     PTZvsSAL_24hr = list(name = "PTZvsSAL_24hr", group1 = "PTZ_24hr", group2 = "SAL_24hr"),
     PTZvsSAL_1hr = list(name = "PTZvsSAL_1hr", group1 = "PTZ_1hr", group2 = "SAL_1hr"),
     `24hrvs1hr_PTZ` = list(name = "24hrvs1hr_PTZ", group1 = "PTZ_24hr", group2 = "PTZ_1hr"),
     `24hrvs1hr_SAL` = list(name = "24hrvs1hr_SAL", group1 = "SAL_24hr", group2 = "SAL_1hr")
 )
 output_dir_base <-"~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"

 for (ct in cell_types) {
     for (comp in comparisons) {
         tryCatch({
             message(sprintf("\nAnalyzing %s - %s", ct, comp$name))

             # subsetting cell type in each condition
             cell_type_data <- subset(combined_object,
                                      subset = cell_type == ct & Sample %in% c(comp$group1, comp$group2)) # cell_type is metadata column name for cell type, Sample is for group

             # test
             # cell_type_data <- subset(combined_object,
             #                          subset = cell_types[1] == ct & Sample %in% c(comparisons$PTZvsSAL_24hr$group1, comparisons$PTZvsSAL_24hr$group2)) # cell_type is metadata column name for cell type, Sample is for group
             #
             # cell_type_data <- PrepSCTFindMarkers(cell_type_data)
             # da_results <- FindMarkers(
             #     object = cell_type_data,
             #     group.by = "Sample",
             #     ident.1 = comparisons$PTZvsSAL_24hr$group1,
             #     ident.2 = comparisons$PTZvsSAL_24hr$group2,
             #     min.pct = 0.05,
             #     test.use = "LR",
             #     latent.vars = "nCount_ATAC"
             # )


             # Set identities
             Idents(cell_type_data) <- "Sample"
             DefaultAssay(combined_object) <- "ATAC"

             # Perform DA peak analysis
             da_results <- FindMarkers(
                 object = cell_type_data,
                 group.by = "Sample",
                 ident.1 = comp$group1,
                 ident.2 = comp$group2,
                 min.pct = 0.05,
                 test.use = "LR",
                 latent.vars = "nCount_ATAC"
             )

             # Add metadata
             da_results$ranges <- rownames(da_results)
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


         }, error = function(e) {
             message(sprintf("Error in analysis: %s", e$message))
         })
     }
 }
