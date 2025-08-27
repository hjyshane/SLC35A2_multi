library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots

output_dir_base <-"~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
combined_object <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/filtered_50_cells_atac_cleaned.qs")

# Read list of significant gene files and DA peak files
sig_gene_files <- list.files(file.path(output_dir_base, "DEG", "sig_csv"),
                             pattern = ".csv$",
                             full.names = TRUE)

sig_peak_files <- list.files(file.path(output_dir_base, "DAG", "Sig_ATAC"),
                             pattern = "Sig_results_.+\\.csv$",
                             full.names = TRUE)

# Create output directory for ATAC-focused plots
output_dir_atac <- file.path(output_dir_base,"DAG", "Plots", "CoveragePlots")
dir.create(output_dir_atac, recursive = TRUE, showWarnings = FALSE)

# Ensure normalization
DefaultAssay(combined_object) <- "RNA"
combined_object <- NormalizeData(combined_object, verbose = FALSE)

# cell type and comparison group list
cell_types <- unique(combined_object$cell_type)

comparisons <- list(
    PTZvsSAL_24hr = list(name = "PTZvsSAL_24hr", group1 = "PTZ_24hr", group2 = "SAL_24hr"),
    PTZvsSAL_1hr = list(name = "PTZvsSAL_1hr", group1 = "PTZ_1hr", group2 = "SAL_1hr"),
    `24hrvs1hr_PTZ` = list(name = "24hrvs1hr_PTZ", group1 = "PTZ_24hr", group2 = "PTZ_1hr"),
    `24hrvs1hr_SAL` = list(name = "24hrvs1hr_SAL", group1 = "SAL_24hr", group2 = "SAL_1hr")
)

# Ensure data is properly normalized
DefaultAssay(combined_object) <- "RNA"
combined_object <- NormalizeData(combined_object, verbose = FALSE)

# Function to create coverage plots focused on differential ATAC peaks
create_atac_peak_coverage_plots <- function(seurat_obj, gene_file, peak_file, comparison_info, output_dir) {

    # Extract file info
    gene_file_name <- basename(gene_file)
    peak_file_name <- basename(peak_file)

    parts <- strsplit(gsub(".csv", "", gene_file_name), "_")[[1]]
    cell_type <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    print(paste("\n Processing:", cell_type, "-", comparison_name))

    # Read differential peaks and genes
    tryCatch({
        sig_peaks <- read.csv(peak_file)
        de_genes <- read.csv(gene_file)

        # Handle different peak file formats
        if("x" %in% colnames(sig_peaks)) {
            # If peaks are stored as single column
            peak_regions <- sig_peaks$x
        } else if(nrow(sig_peaks) > 0) {
            # If peaks are stored as row names or first column
            peak_regions <- rownames(sig_peaks)
            if(is.null(peak_regions) || all(peak_regions == as.character(1:nrow(sig_peaks)))) {
                peak_regions <- sig_peaks[,1]
            }
        } else {
            print("     No peaks found in file")
            return(NULL)
        }

        print(paste("Found", length(peak_regions), "differential peaks"))
        print(paste("Found", nrow(de_genes), "differential genes"))

        # Take top peaks (by significance if available)
        if(length(peak_regions) > 20) {
            peak_regions <- head(peak_regions, 20)  # Top 20 peaks
        }

        # Get top DEGs for context
        top_genes <- de_genes %>%
            dplyr::filter(!is.na(gene) & gene != "") %>%
            dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.3) %>%
            arrange(desc(abs(avg_log2FC))) %>%
            head(10)

        # Create cell subset
        cell_subset <- subset(
            seurat_obj,
            subset = cell_type == cell_type &
                Sample %in% c(comparison_info$group1, comparison_info$group2)
        )

        print(paste("Cells in subset:", ncol(cell_subset)))

        if(ncol(cell_subset) < 10) {
            print("Too few cells, skipping...")
            return(NULL)
        }

        # Normalize if needed
        DefaultAssay(cell_subset) <- "RNA"
        if(sum(cell_subset[["RNA"]]@layers$data) == 0) {
            cell_subset <- NormalizeData(cell_subset, verbose = FALSE)
        }

        DefaultAssay(cell_subset) <- "ATAC"
        Idents(cell_subset) <- "Sample"

        successful_plots <- 0

        # Method 1: Plot differential ATAC peaks directly
        print("   Creating peak-focused coverage plots...")

        for(i in 1:min(10, length(peak_regions))) {
            peak_region <- peak_regions[i]

            # Clean peak region name
            peak_clean <- gsub("[^A-Za-z0-9-]", "_", peak_region)

            tryCatch({
                # Create coverage plot for the peak region
                coverage_plot <- CoveragePlot(
                    object = cell_subset,
                    region = peak_region,
                    extend.upstream = 10000,
                    extend.downstream = 10000,
                    annotation = TRUE,
                    peaks = TRUE
                )

                plot_filename <- paste0(
                    "Peak_Coverage_", peak_clean, "_",
                    comparison_info$group1, "vs", comparison_info$group2, "_", cell_type, ".png"
                )

                ggsave(
                    filename = file.path(output_dir, plot_filename),
                    plot = coverage_plot,
                    width = 14,
                    height = 8,
                    dpi = 300
                )

                successful_plots <- successful_plots + 1

                if(i %% 3 == 0) {
                    print(paste(i, "peak plots created..."))
                }

            }, error = function(e) {
                print(paste("Peak plot failed for", peak_region, ":", e$message))
            })
        }

        # Method 2: Create combined peak + gene plots for top genes
        print("Creating gene-centered plots with peak highlights...")

        for(i in 1:min(5, nrow(top_genes))) {
            gene_name <- top_genes$gene[i]

            tryCatch({
                # Find peaks near this gene
                annotation <- Annotation(cell_subset)
                gene_coords <- annotation[mcols(annotation)$gene_name == gene_name &
                                              mcols(annotation)$type == "gene"]

                if(length(gene_coords) > 0) {
                    # Create coverage plot centered on gene
                    gene_coverage_plot <- CoveragePlot(
                        object = cell_subset,
                        region = gene_name,
                        features = gene_name,
                        expression.assay = "RNA",
                        extend.upstream = 25000,
                        extend.downstream = 25000,
                        annotation = TRUE,
                        peaks = TRUE,
                        heights = c(3, 1)
                    )

                    plot_filename_gene <- paste0(
                        "Gene_Peak_Coverage_", gene_name, "_",
                        comparison_info$group1, "vs", comparison_info$group2, "_", cell_type, ".png"
                    )

                    ggsave(
                        filename = file.path(output_dir, plot_filename_gene),
                        plot = gene_coverage_plot,
                        width = 14,
                        height = 10,
                        dpi = 300
                    )

                    successful_plots <- successful_plots + 1
                }

            }, error = function(e) {
                print(paste("Gene-peak plot failed for", gene_name, ":", e$message))
            })
        }

        print(paste("Created", successful_plots, "total plots"))
        return(successful_plots)

    }, error = function(e) {
        print(paste("Error processing files:", e$message))
        return(NULL)
    })
}

# Match gene and peak files and process them
process_matched_files <- function(seurat_obj, gene_files, peak_files, comparisons, output_dir) {

    total_plots <- 0
    processed_pairs <- 0

    for(gene_file in gene_files) {
        # Extract file info
        gene_file_name <- basename(gene_file)
        parts <- strsplit(gsub(".csv", "", gene_file_name), "_")[[1]]
        cell_type <- parts[1]
        comparison_name <- paste(parts[-1], collapse = "_")

        # Skip if comparison not recognized
        if(!comparison_name %in% names(comparisons)) {
            print(paste("Skipping unknown comparison:", comparison_name))
            next
        }

        # Find matching peak file
        peak_file_pattern <- paste0("Sig_results_", cell_type, "_", comparison_name, ".csv")
        matching_peak_files <- peak_files[grepl(peak_file_pattern, peak_files)]

        if(length(matching_peak_files) == 0) {
            print(paste("  No matching peak file for:", cell_type, comparison_name))
            next
        }

        peak_file <- matching_peak_files[1]  # Take first match

        # Process this gene-peak pair
        result <- create_atac_peak_coverage_plots(
            seurat_obj = seurat_obj,
            gene_file = gene_file,
            peak_file = peak_file,
            comparison_info = comparisons[[comparison_name]],
            output_dir = output_dir
        )

        if(!is.null(result) && result > 0) {
            total_plots <- total_plots + result
            processed_pairs <- processed_pairs + 1
        }
    }

    return(list(processed = processed_pairs, total_plots = total_plots))
}

# Run the analysis
atac_results <- process_matched_files(
    seurat_obj = combined_object,
    gene_files = sig_gene_files,
    peak_files = sig_peak_files,
    comparisons = comparisons,
    output_dir = output_dir_atac
)

# Create summary for ATAC peaks
create_atac_summary <- function(output_dir) {
    # Count different types of plots
    peak_plots <- list.files(output_dir, pattern = "^Peak_Coverage_.*\\.png$")
    gene_peak_plots <- list.files(output_dir, pattern = "^Gene_Peak_Coverage_.*\\.png$")

    summary_df <- data.frame(
        Plot_Type = c("Peak-Focused", "Gene-Peak Combined", "Total"),
        Count = c(length(peak_plots), length(gene_peak_plots), length(peak_plots) + length(gene_peak_plots))
    )

    write.csv(summary_df, file.path(output_dir, "atac_coverage_summary.csv"), row.names = FALSE)

    print("ATAC Coverage Plot Types:")
    print(summary_df)
}

create_atac_summary(output_dir_atac)
