# Parallelization for process
library(future)
options(future.globals.maxSize = 128000 * 1024 ^ 2)  # Increase memory limit
library(parallel)  # Base R parallelization package

# Single-cell analysis (scRNA-seq & scATAC-seq)
library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(sctransform)    # Normalization for scRNA-seq
library(DoubletFinder)  # Detects doublets in single-cell datasets
library(scCustomize)    # Enhances Seurat visualization and analysis
library(SingleCellExperiment)

# Data manipulation & utilities
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(data.table)  # High-speed alternative to data frames
library(Matrix)      # Sparse matrix operations (important for single-cell data)
library(jsonlite)

# Visualization & Plotting
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots
library(cowplot)          # Additional ggplot utilities
library(ggrepel)          # Prevents text overlap in plots
library(ggvenn)           # Venn diagrams
library(ggVennDiagram)    # Another Venn diagram package
library(VennDiagram)      # Classic Venn diagrams
library(ggpubr)           # Publication-ready plots
library(ComplexHeatmap)   # Advanced heatmaps
library(circlize)         # Circular visualizations
library(pheatmap)         # Heatmaps
library(RColorBrewer)     # Color palettes
library(fields)           # Spatial plots and interpolation
library(enrichplot)       # Visualization for enrichment analysis
library(eulerr)           # Euler diagrams
library(gridExtra)        # Arrange multiple plots
library(corrplot)         # Correlation heatmaps
library(viridis)          # Colorblind-friendly palettes
library(scales)           # scale plot
library(EnhancedVolcano)
library(RColorBrewer)

# Clustering & Dimensionality Reduction
library(clustree)  # Clustering trees for resolution selection
library(igraph)    # Graph-based clustering

# Differential Expression & Statistical Analysis
library(DESeq2)   # Differential expression analysis
library(MAST)
library(rstatix)  # Statistical functions for ggplot2
library(ROCR)     # ROC curves for model evaluation

# Gene Ontology & Pathway Analysis
library(clusterProfiler)  # Functional enrichment analysis
library(org.Mm.eg.db)     # Mouse genome annotations
library(DOSE)             # Disease ontology enrichment
library(rrvgo)            # Reduce GO terms redundancy
library(treemap)          # Tree visualization for GO terms
library(AnnotationDbi)    # Database interface
library(GO.db)            # Gene Ontology database
library(ontologyIndex)    # General ontology functions
library(fgsea)            # fGSEA analysis
library(msigdbr)          # ref for fGSEA

# Pseudotime trajectory analysis
library(slingshot)
library(tradeSeq)

# Genomics & Motif Analysis
library(BSgenome.Mmusculus.UCSC.mm10)  # Mouse genome
library(rtracklayer)                   # Import/export genomic data
library(JASPAR2020)                    # Transcription factor motifs
library(TFBSTools)                      # Analyze transcription factor binding sites

# File Handling & Misc
library(openxlsx)  # Read/write Excel files
dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")  # Load HDF5 library
library(hdf5r)     # HDF5 file format for large datasets

# Interactive Applications
library(shiny)  # Web-based interactive applications
library(qs)

set.seed(42)

module <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/ModuleTest/qsave/filtered_50_sct.qs")

# Set the default assay to RNA
DefaultAssay(module) <- "RNA"

# Define lists of cell types and comparisons
cell_types <- unique(module$cell_type)
comparisons_to_analyze <- list(
    "PTZvsSAL_24hr" = list(name = "PTZvsSAL_24hr", group1 = "PTZ_24hr", group2 = "SAL_24hr"),
    "PTZvsSAL_1hr" = list(name = "PTZvsSAL_1hr", group1 = "PTZ_1hr", group2 = "SAL_1hr"),
    "24hrvs1hr_PTZ" = list(name = "24hrvs1hr_PTZ", group1 = "PTZ_24hr", group2 = "PTZ_1hr"),
    "24hrvs1hr_SAL" = list(name = "24hrvs1hr_SAL", group1 = "SAL_24hr", group2 = "SAL_1hr")
)

output_dir_base <- "~/PTZ_ATAC_scRNA_072024/WIP/ModuleTest/"

# Main analysis loop with error handling
for (ct in cell_types) {
    for (comp in comparisons_to_analyze) {
        # Subset data
        # subsetting cell type in each condition
        cell_type_data <- subset(module,
                                 subset = cell_type == ct & Sample %in% c(comp$group1, comp$group2))

        # Set identities
        Idents(cell_type_data) <- "Sample"

        # Get cell counts
        cell_count_group1 <- sum(Idents(cell_type_data) == comp$group1)
        cell_count_group2 <- sum(Idents(cell_type_data) == comp$group2)

        message(sprintf("Cell counts - %s: %d, %s: %d", # %s for string, %d for integer
                        comp$group1, cell_count_group1,
                        comp$group2, cell_count_group2))

        # Check minimum cell threshold (e.g., 3 cells per group)
        if (cell_count_group1 >= 3 && cell_count_group2 >= 3) {
            # Perform DE analysis
            de_results <- FindMarkers(
                object = cell_type_data,
                assay = "SCT", # RNA for MAST, SCT for wilcox
                slot = "data",
                ident.1 = comp$group1,
                ident.2 = comp$group2,
                min.pct = 0.25,
                logfc.threshold = 0.0,
                # return.thresh = 1.01, # depreciated
                test.use = "wilcox"
            )

            # Add metadata
            de_results$gene <- rownames(de_results)
            de_results$cell_type <- ct
            de_results$comparison <- comp$name

            # Save results
            csv_filename <- sprintf("DE_results_%s_%s.csv", ct, comp$name)
            write.csv(
                de_results,
                file = file.path(output_dir_base, "DE_gene", csv_filename),
                row.names = FALSE
            )

        }
    }
}


# file read
csv_files <- list.files(file.path(output_dir_base,  "DE_gene"), pattern = "DE_results_.+\\.csv$", full.names = TRUE)

# loop
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Creating volcano plot for", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    if (nrow(de_results) == 0) {
        next
    } else {
        # Create volcano plot
        volcano_plot <- EnhancedVolcano(de_results,
                                        lab = de_results$gene,
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        title = paste('Volcano plot -', ct),
                                        subtitle = safe_comparison_name,
                                        pCutoff = 0.05,
                                        FCcutoff = 0.2,
                                        pointSize = 1.0,
                                        labSize = 3.0
        )

        # Save volcano plot with safe filename
        plot_file_name <- paste0("Volcano_plot_", ct, "_", safe_comparison_name, ".png")
        ggsave(file.path(output_dir_base, "VolcanoPlots", plot_file_name), volcano_plot, width = 12, height = 10)
        print(paste("Saved volcano plot:", plot_file_name))
    }}

# sig calculation function
get_significant_genes <- function(de_results, p_cutoff = 0.05, fc_cutoff = 0.2) {
    significant_genes <- de_results %>%
        mutate(Regulation = case_when(
            avg_log2FC > 0 ~ "Up",
            avg_log2FC < 0 ~ "Down",
            TRUE ~ "NoChange")) %>%
        filter(p_val_adj < p_cutoff & abs(avg_log2FC) >= fc_cutoff) %>%
        arrange(p_val_adj)
    return(significant_genes)
}

# Initialize a list to store summaries
summary_list <- list()

# loop for sig gene and save
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Processing", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    # Get significant genes
    significant_genes <- get_significant_genes(de_results)

    # Save significant genes with safe filename
    sig_genes_filename <- paste0("Significant_genes_", ct, "_", safe_comparison_name, ".csv")

    write.csv(significant_genes, file = file.path(output_dir_base, "Sig_gene", sig_genes_filename), row.names = FALSE)
    print(paste("Saved significant genes:", sig_genes_filename))

}

create_deg_summary <- function(deg_dir) {
    # Get all DEG files
    deg_files <- list.files(deg_dir, pattern = "DE_results_.+\\.csv$", full.names = TRUE)

    # Initialize list to store summaries
    deg_summaries <- list()

    for (file in deg_files) {

        # Extract metadata from filename
        file_name <- basename(file)
        parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
        ct <- parts[1]
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
            Cell_Type = ct,
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
complete_summary <- create_deg_summary(deg_dir = file.path(output_dir_base, "DE_gene"))

# Save summary
write.csv(complete_summary, file.path(output_dir_base, "DEG_Complete_Analysis_Summary.csv"), row.names = FALSE)
rm("module")

script <- readRDS("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/RDS_mid/combined_merged_processed_multimodal_processed_linked_annotated_50_cells_sct.rds")

# Set the default assay to RNA
DefaultAssay(script) <- "RNA"

# Define lists of cell types and comparisons
cell_types <- unique(script$cell_type)
comparisons_to_analyze <- list(
    "PTZvsSAL_24hr" = list(name = "PTZvsSAL_24hr", group1 = "PTZ_24hr", group2 = "SAL_24hr"),
    "PTZvsSAL_1hr" = list(name = "PTZvsSAL_1hr", group1 = "PTZ_1hr", group2 = "SAL_1hr"),
    "24hrvs1hr_PTZ" = list(name = "24hrvs1hr_PTZ", group1 = "PTZ_24hr", group2 = "PTZ_1hr"),
    "24hrvs1hr_SAL" = list(name = "24hrvs1hr_SAL", group1 = "SAL_24hr", group2 = "SAL_1hr")
)

output_dir_base <- "~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/"

# Main analysis loop with error handling
for (ct in cell_types) {
    for (comp in comparisons_to_analyze) {
        # Subset data
        # subsetting cell type in each condition
        cell_type_data <- subset(script,
                                 subset = cell_type == ct & sample %in% c(comp$group1, comp$group2))

        # Set identities
        Idents(cell_type_data) <- "sample"

        # Get cell counts
        cell_count_group1 <- sum(Idents(cell_type_data) == comp$group1)
        cell_count_group2 <- sum(Idents(cell_type_data) == comp$group2)

        message(sprintf("Cell counts - %s: %d, %s: %d", # %s for string, %d for integer
                        comp$group1, cell_count_group1,
                        comp$group2, cell_count_group2))

        # Check minimum cell threshold (e.g., 3 cells per group)
        if (cell_count_group1 >= 3 && cell_count_group2 >= 3) {
            # Perform DE analysis
            de_results <- FindMarkers(
                object = cell_type_data,
                assay = "SCT", # RNA for MAST, SCT for wilcox
                slot = "data",
                ident.1 = comp$group1,
                ident.2 = comp$group2,
                min.pct = 0.25,
                logfc.threshold = 0.0,
                # return.thresh = 1.01, # depreciated
                test.use = "wilcox"
            )

            # Add metadata
            de_results$gene <- rownames(de_results)
            de_results$cell_type <- ct
            de_results$comparison <- comp$name

            # Save results
            csv_filename <- sprintf("DE_results_%s_%s.csv", ct, comp$name)
            write.csv(
                de_results,
                file = file.path(output_dir_base, "DE_gene", csv_filename),
                row.names = FALSE
            )

        }
    }
}

# file read
csv_files <- list.files(file.path(output_dir_base,  "DE_gene"), pattern = "DE_results_.+\\.csv$", full.names = TRUE)

# loop
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Creating volcano plot for", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    if (nrow(de_results) == 0) {
        next
    } else {
        # Create volcano plot
        volcano_plot <- EnhancedVolcano(de_results,
                                        lab = de_results$gene,
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        title = paste('Volcano plot -', ct),
                                        subtitle = safe_comparison_name,
                                        pCutoff = 0.05,
                                        FCcutoff = 0.2,
                                        pointSize = 1.0,
                                        labSize = 3.0
        )

        # Save volcano plot with safe filename
        plot_file_name <- paste0("Volcano_plot_", ct, "_", safe_comparison_name, ".png")
        ggsave(file.path(output_dir_base, "VolcanoPlots", plot_file_name), volcano_plot, width = 12, height = 10)
        print(paste("Saved volcano plot:", plot_file_name))
    }}


# Initialize a list to store summaries
summary_list <- list()

# loop for sig gene and save
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Processing", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    # Get significant genes
    significant_genes <- get_significant_genes(de_results)

    # Save significant genes with safe filename
    sig_genes_filename <- paste0("Significant_genes_", ct, "_", safe_comparison_name, ".csv")

    write.csv(significant_genes, file = file.path(output_dir_base, "Sig_gene", sig_genes_filename), row.names = FALSE)
    print(paste("Saved significant genes:", sig_genes_filename))

}

# Generate summary
complete_summary <- create_deg_summary(deg_dir = file.path(output_dir_base, "DE_gene"))

# Save summary
write.csv(complete_summary, file.path(output_dir_base, "DEG_Complete_Analysis_Summary.csv"), row.names = FALSE)
rm(script)


module <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/ModuleTest/qsave/filtered_50_sctall.qs")

output_dir_base <- "~/PTZ_ATAC_scRNA_072024/WIP/ModuleTest/"

# Main analysis loop with error handling
for (ct in cell_types) {
    for (comp in comparisons_to_analyze) {
        # Subset data
        # subsetting cell type in each condition
        cell_type_data <- subset(module,
                                 subset = cell_type == ct & Sample %in% c(comp$group1, comp$group2))

        # Set identities
        Idents(cell_type_data) <- "Sample"

        # Get cell counts
        cell_count_group1 <- sum(Idents(cell_type_data) == comp$group1)
        cell_count_group2 <- sum(Idents(cell_type_data) == comp$group2)

        message(sprintf("Cell counts - %s: %d, %s: %d", # %s for string, %d for integer
                        comp$group1, cell_count_group1,
                        comp$group2, cell_count_group2))

        # Check minimum cell threshold (e.g., 3 cells per group)
        if (cell_count_group1 >= 3 && cell_count_group2 >= 3) {
            # Perform DE analysis
            de_results <- FindMarkers(
                object = cell_type_data,
                assay = "SCT", # RNA for MAST, SCT for wilcox
                slot = "data",
                ident.1 = comp$group1,
                ident.2 = comp$group2,
                min.pct = 0.25,
                logfc.threshold = 0.0,
                # return.thresh = 1.01, # depreciated
                test.use = "wilcox"
            )

            # Add metadata
            de_results$gene <- rownames(de_results)
            de_results$cell_type <- ct
            de_results$comparison <- comp$name

            # Save results
            csv_filename <- sprintf("DE_results_%s_%s.csv", ct, comp$name)
            write.csv(
                de_results,
                file = file.path(output_dir_base, "DE_gene_all", csv_filename),
                row.names = FALSE
            )

        }
    }
}


# file read
csv_files <- list.files(file.path(output_dir_base,  "DE_gene_all"), pattern = "DE_results_.+\\.csv$", full.names = TRUE)

# loop
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Creating volcano plot for", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    if (nrow(de_results) == 0) {
        next
    } else {
        # Create volcano plot
        volcano_plot <- EnhancedVolcano(de_results,
                                        lab = de_results$gene,
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        title = paste('Volcano plot -', ct),
                                        subtitle = safe_comparison_name,
                                        pCutoff = 0.05,
                                        FCcutoff = 0.2,
                                        pointSize = 1.0,
                                        labSize = 3.0
        )

        # Save volcano plot with safe filename
        plot_file_name <- paste0("Volcano_plot_", ct, "_", safe_comparison_name, ".png")
        ggsave(file.path(output_dir_base, "VolcanoPlots_all", plot_file_name), volcano_plot, width = 12, height = 10)
        print(paste("Saved volcano plot:", plot_file_name))
    }}

# Initialize a list to store summaries
summary_list <- list()

# loop for sig gene and save
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Processing", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    # Get significant genes
    significant_genes <- get_significant_genes(de_results)

    # Save significant genes with safe filename
    sig_genes_filename <- paste0("Significant_genes_", ct, "_", safe_comparison_name, ".csv")

    write.csv(significant_genes, file = file.path(output_dir_base, "Sig_gene_all", sig_genes_filename), row.names = FALSE)
    print(paste("Saved significant genes:", sig_genes_filename))

}

create_deg_summary <- function(deg_dir) {
    # Get all DEG files
    deg_files <- list.files(deg_dir, pattern = "DE_results_.+\\.csv$", full.names = TRUE)

    # Initialize list to store summaries
    deg_summaries <- list()

    for (file in deg_files) {

        # Extract metadata from filename
        file_name <- basename(file)
        parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
        ct <- parts[1]
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
            Cell_Type = ct,
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
complete_summary <- create_deg_summary(deg_dir = file.path(output_dir_base, "DE_gene_all"))

# Save summary
write.csv(complete_summary, file.path(output_dir_base, "DEG_Complete_Analysis_Summary_all.csv"), row.names = FALSE)
rm(module)

script <- readRDS("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/RDS_mid/combined_merged_processed_multimodal_processed_linked_annotated_50_cells_sctall.rds")

# Set the default assay to RNA
DefaultAssay(script) <- "RNA"

# Define lists of cell types and comparisons
cell_types <- unique(script$cell_type)
comparisons_to_analyze <- list(
    "PTZvsSAL_24hr" = list(name = "PTZvsSAL_24hr", group1 = "PTZ_24hr", group2 = "SAL_24hr"),
    "PTZvsSAL_1hr" = list(name = "PTZvsSAL_1hr", group1 = "PTZ_1hr", group2 = "SAL_1hr"),
    "24hrvs1hr_PTZ" = list(name = "24hrvs1hr_PTZ", group1 = "PTZ_24hr", group2 = "PTZ_1hr"),
    "24hrvs1hr_SAL" = list(name = "24hrvs1hr_SAL", group1 = "SAL_24hr", group2 = "SAL_1hr")
)

output_dir_base <- "~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/"

# Main analysis loop with error handling
for (ct in cell_types) {
    for (comp in comparisons_to_analyze) {
        # Subset data
        # subsetting cell type in each condition
        cell_type_data <- subset(script,
                                 subset = cell_type == ct & sample %in% c(comp$group1, comp$group2))

        # Set identities
        Idents(cell_type_data) <- "sample"

        # Get cell counts
        cell_count_group1 <- sum(Idents(cell_type_data) == comp$group1)
        cell_count_group2 <- sum(Idents(cell_type_data) == comp$group2)

        message(sprintf("Cell counts - %s: %d, %s: %d", # %s for string, %d for integer
                        comp$group1, cell_count_group1,
                        comp$group2, cell_count_group2))

        # Check minimum cell threshold (e.g., 3 cells per group)
        if (cell_count_group1 >= 3 && cell_count_group2 >= 3) {
            # Perform DE analysis
            de_results <- FindMarkers(
                object = cell_type_data,
                assay = "SCT", # RNA for MAST, SCT for wilcox
                slot = "data",
                ident.1 = comp$group1,
                ident.2 = comp$group2,
                min.pct = 0.25,
                logfc.threshold = 0.0,
                # return.thresh = 1.01, # depreciated
                test.use = "wilcox"
            )

            # Add metadata
            de_results$gene <- rownames(de_results)
            de_results$cell_type <- ct
            de_results$comparison <- comp$name

            # Save results
            csv_filename <- sprintf("DE_results_%s_%s.csv", ct, comp$name)
            write.csv(
                de_results,
                file = file.path(output_dir_base, "DE_gene_all", csv_filename),
                row.names = FALSE
            )

        }
    }
}

# file read
csv_files <- list.files(file.path(output_dir_base,  "DE_gene_all"), pattern = "DE_results_.+\\.csv$", full.names = TRUE)

# loop
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Creating volcano plot for", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    if (nrow(de_results) == 0) {
        next
    } else {
        # Create volcano plot
        volcano_plot <- EnhancedVolcano(de_results,
                                        lab = de_results$gene,
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        title = paste('Volcano plot -', ct),
                                        subtitle = safe_comparison_name,
                                        pCutoff = 0.05,
                                        FCcutoff = 0.2,
                                        pointSize = 1.0,
                                        labSize = 3.0
        )

        # Save volcano plot with safe filename
        plot_file_name <- paste0("Volcano_plot_", ct, "_", safe_comparison_name, ".png")
        ggsave(file.path(output_dir_base, "VolcanoPlots_all", plot_file_name), volcano_plot, width = 12, height = 10)
        print(paste("Saved volcano plot:", plot_file_name))
    }}


# Initialize a list to store summaries
summary_list <- list()

# loop for sig gene and save
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    ct <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Processing", ct, "-", comparison_name))

    # Read CSV file
    de_results <- read.csv(csv_file)

    # Get significant genes
    significant_genes <- get_significant_genes(de_results)

    # Save significant genes with safe filename
    sig_genes_filename <- paste0("Significant_genes_", ct, "_", safe_comparison_name, ".csv")

    write.csv(significant_genes, file = file.path(output_dir_base, "Sig_gene_all", sig_genes_filename), row.names = FALSE)
    print(paste("Saved significant genes:", sig_genes_filename))

}

# Generate summary
complete_summary <- create_deg_summary(deg_dir = file.path(output_dir_base, "DE_gene_all"))

# Save summary
write.csv(complete_summary, file.path(output_dir_base, "DEG_Complete_Analysis_Summary_all.csv"), row.names = FALSE)

