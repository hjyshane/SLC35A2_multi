library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate

deg_dir <-  "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/sig_csv/"
da_dir <-  "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DAG/Sig_ATAC/"
output_dir <-  "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DAG/Plots/CoveragePlots_test/"
dir.create(output_dir)

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