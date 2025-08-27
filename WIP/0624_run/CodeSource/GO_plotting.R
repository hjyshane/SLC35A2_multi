count_gene_go_terms <- function(input_dir, output_dir_plot, output_dir_csv, output_file) {

    # Read all GO result files
    go_files <- list.files(input_dir, pattern = ".csv$", full.names = TRUE)

    # Create data frames to store results for each condition
    all_counts <- data.frame()
    all_counts_top <- data.frame()

    for (file in go_files) {
        # Extract condition name from filename
        condition <- gsub(".csv", "", basename(file))

        # Read GO results
        go_data <- read_csv(file)

        # Skip files with no data
        if (nrow(go_data) == 0) {
            message(sprintf("Skipping file %s: No data", file))
            next
        }

        # Ensure geneID column exists
        if (!"geneID" %in% colnames(go_data)) {
            message(sprintf("Skipping file %s: Missing 'geneID' column", file))
            next
        }

        # Ensure geneID is character type
        go_data$geneID <- as.character(go_data$geneID)

        # Select top 20 GO terms by adjusted p-value
        top_terms <- go_data %>%
            arrange(p.adjust) %>%
            slice_min(order_by = p.adjust, n = min(20, nrow(go_data)))

        # Convert geneID column from merged string format to long format
        gene_go_pairs <- go_data %>%
            dplyr::select(ID, geneID) %>%
            separate_rows(geneID, sep = "/")  # Expands "Tuba1a/Tubb5" to multiple rows

        # Count occurrences of each gene across GO terms
        gene_counts <- table(gene_go_pairs$geneID)
        go_id_map <- gene_go_pairs  # Store gene-to-GO term mapping

        # Repeat for top 20 GO terms
        gene_go_pairs_top <- top_terms %>%
            dplyr::select(ID, geneID) %>%
            separate_rows(geneID, sep = "/")

        gene_counts_top <- table(gene_go_pairs_top$geneID)
        go_id_map_top <- gene_go_pairs_top

        # Convert counts to data frames
        count_df <- data.frame(
            condition = condition,
            gene = names(gene_counts),
            go_term_count = as.numeric(gene_counts)
        ) %>% merge(go_id_map, by.x = "gene", by.y = "geneID", all.x = TRUE)

        count_df_top <- data.frame(
            condition = condition,
            gene = names(gene_counts_top),
            go_term_count = as.numeric(gene_counts_top)
        ) %>% merge(go_id_map_top, by.x = "gene", by.y = "geneID", all.x = TRUE)

        # Append results to main data frames
        all_counts <- bind_rows(all_counts, count_df)
        all_counts_top <- bind_rows(all_counts_top, count_df_top)

        # Analyze gene frequency for bar plots and gene reports
        gene_freq_df <- data.frame(gene = names(gene_counts), frequency = as.numeric(gene_counts)) %>%
            arrange(desc(frequency))

        # Create bar plot for top 30 genes
        plot <- ggplot(head(gene_freq_df, 30), aes(x = reorder(gene, frequency), y = frequency)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            geom_text(aes(label = frequency), hjust = -0.2, size = 3.5) +  # Add count labels
            coord_flip() +
            labs(title = paste("Top 30 Genes by GO Term Frequency -", condition),
                 x = "Gene",
                 y = "Number of GO Terms") +
            theme_minimal() +
            theme(plot.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10))

        # Save bar plot
        ggsave(filename = file.path(output_dir_plot, paste0("gene_frequency_", condition, ".png")),
               plot = plot, width = 10, height = 8, dpi = 300, bg = "white")

        # Create detailed gene report
        gene_report <- data.frame()

        for (gene in unique(gene_go_pairs$geneID)) {
            # Extract all GO terms for the gene
            gene_terms <- go_data %>%
                filter(grepl(paste0("\\b", gene, "\\b"), geneID)) %>%
                dplyr::select(ID, Description, p.adjust)

            # Create report entry
            report_entry <- data.frame(
                gene = gene,
                frequency = gene_counts[gene],
                go_terms = paste(
                    paste0(
                        gene_terms$Description, " (", gene_terms$ID,
                        ", p=", formatC(gene_terms$p.adjust, format = "e", digits = 2), ")"
                    ),
                    collapse = "; "
                )
            )
            gene_report <- bind_rows(gene_report, report_entry)
        }

        # Sort gene report by frequency
        gene_report <- gene_report %>% arrange(desc(frequency))

        # Save gene report
        write_csv(gene_report, file.path(output_dir_csv, paste0("gene_report_", condition, ".csv")))
    }

    # Write overall results to CSV
    write_csv(all_counts, file.path(input_dir, "gene_go_term_counts.csv"))
    write_csv(all_counts_top, file.path(input_dir, "gene_go_term_counts_top.csv"))

    # Combine all gene reports into one file
    file_list <- list.files(path = output_dir_csv, pattern = "*.csv", full.names = TRUE)

    combined_df <- lapply(file_list, function(file) {
        df <- read_csv(file)
        df$Source_File <- basename(file)  # Track source file
        return(df)
    }) %>% bind_rows()

    # Save combined results
    write_csv(combined_df, output_file)

    # Return results as a list
    return(list(all_counts = all_counts, all_counts_top = all_counts_top))
}

# Function to count gene occurrences and create a gene report for the top 20 GO terms
count_gene_top20_with_report <- function(input_dir, output_dir_plot, output_dir_csv, output_file) {

    # Read all GO result files
    go_files <- list.files(input_dir, pattern = ".csv$", full.names = TRUE)

    # Create a data frame to store results
    all_counts_top <- data.frame()

    for (file in go_files) {
        # Extract condition name from filename
        condition <- gsub(".csv", "", basename(file))

        # Read GO results
        go_data <- read_csv(file)

        # Skip files with no data
        if (nrow(go_data) == 0) {
            message(sprintf("Skipping file %s: No data", file))
            next
        }

        # Ensure 'geneID' column exists
        if (!"geneID" %in% colnames(go_data)) {
            message(sprintf("Skipping file %s: Missing 'geneID' column", file))
            next
        }

        # Ensure geneID is character type
        go_data$geneID <- as.character(go_data$geneID)

        # Select the top 20 GO terms by adjusted p-value
        top_terms <- go_data %>%
            arrange(p.adjust) %>%
            slice_min(order_by = p.adjust, n = min(20, nrow(go_data)))

        # Convert geneID column from merged string format to long format (only for top 20 GO terms)
        gene_go_pairs_top <- top_terms %>%
            dplyr::select(ID, Description, p.adjust, geneID) %>%
            separate_rows(geneID, sep = "/")  # Expands "Tuba1a/Tubb5" into separate rows

        # Count occurrences of each gene across the top 20 GO terms
        gene_counts_top <- table(gene_go_pairs_top$geneID)

        # Convert counts to a data frame
        count_df_top <- data.frame(
            condition = condition,
            gene = names(gene_counts_top),
            go_term_count = as.numeric(gene_counts_top)
        )

        # Append results to main data frame
        all_counts_top <- bind_rows(all_counts_top, count_df_top)

        # Analyze gene frequency for bar plots
        gene_freq_df <- count_df_top %>% arrange(desc(go_term_count))

        # Create bar plot for top 30 genes
        plot <- ggplot(head(gene_freq_df, 30), aes(x = reorder(gene, go_term_count), y = go_term_count)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            geom_text(aes(label = go_term_count), hjust = -0.2, size = 3.5) +  # Add count labels
            coord_flip() +
            labs(title = paste("Top 30 Genes in Top 20 GO Terms -", condition),
                 x = "Gene",
                 y = "Number of GO Terms") +
            theme_minimal() +
            theme(plot.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10))

        # Save bar plot
        ggsave(filename = file.path(output_dir_plot, paste0("gene_frequency_top20_", condition, ".png")),
               plot = plot, width = 10, height = 8, dpi = 300, bg = "white")

        # Generate detailed gene report for top 20 GO terms
        gene_report <- data.frame()

        for (gene in unique(gene_go_pairs_top$geneID)) {
            # Extract GO terms associated with the gene
            gene_terms <- gene_go_pairs_top %>%
                filter(geneID == gene) %>%
                dplyr::select(ID, Description, p.adjust)

            # Create report entry
            report_entry <- data.frame(
                gene = gene,
                frequency = gene_counts_top[gene],
                go_terms = paste(
                    paste0(
                        gene_terms$Description, " (", gene_terms$ID,
                        ", p=", formatC(gene_terms$p.adjust, format = "e", digits = 2), ")"
                    ),
                    collapse = "; "
                )
            )
            gene_report <- bind_rows(gene_report, report_entry)
        }

        # Sort gene report by frequency
        gene_report <- gene_report %>% arrange(desc(frequency))

        # Save gene report
        write_csv(gene_report, file.path(output_dir_csv, paste0("gene_report_top20_", condition, ".csv")))
    }

    # Write overall top gene counts to CSV
    write_csv(all_counts_top, file.path(output_dir_csv, "gene_go_term_counts_top20.csv"))

    # Return as list
    return(all_counts_top)
}

go_csv <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/GO/csv/"
output_dir_gf_csv <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/GO/Gene_frequency/summary/"
output_dir_gf_plot <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/GO/Gene_frequency/plots/"
output_file <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/GO/Gene_frequency/summary/gene_combined.csv"
output_file_top20 <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/GO/Gene_frequency/summary/gene_combined_top20.csv"

# Create necessary directories
dir.create(output_dir_gf_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_gf_csv, recursive = TRUE, showWarnings = FALSE)

# Run the function
count_gene_go_terms(input_dir = go_csv,
                    output_dir_plot = output_dir_gf_plot,
                    output_dir_csv = output_dir_gf_csv,
                    output_file = output_file)

# Run the function
count_gene_top20_with_report(input_dir = go_csv,
                             output_dir_plot = output_dir_gf_plot,
                             output_dir_csv = output_dir_gf_csv,
                             output_file = output_file_top20)

output_dir_base <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
csv_files <- list.files(bg_dir, pattern = ".csv$", full.names = TRUE)

all_counts_top20 <- read_csv("GO/Gene_frequency/summary/gene_go_term_counts_top20.csv")

# loop
for (csv_file in csv_files) {
    # Extract cell type and comparison name from file name
    file_name <- basename(csv_file)
    parts <- strsplit(gsub("DE_results_|.csv", "", file_name), "_")[[1]]
    cell_type <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    # Replace colons with underscores in comparison name
    safe_comparison_name <- gsub("[:]", "_", comparison_name)
    print(paste("Creating volcano plot for", cell_type, "-", comparison_name))

    condition_name <- gsub(".csv", "", file_name)[[1]]
    selected_genes <- all_counts_top20 %>%
        filter(condition == condition_name) %>%
        pull(gene) %>%
        unique()

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
                                        title = paste('Volcano plot -', cell_type),
                                        subtitle = safe_comparison_name,
                                        pCutoff = 0.05,
                                        FCcutoff = 0.2,
                                        pointSize = 1.0,
                                        labSize = 3.0,
                                        selectLab = selected_genes)

        # Save volcano plot with safe filename
        plot_file_name <- paste0("Volcano_plot_", cell_type, "_", safe_comparison_name, ".png")
        ggsave(paste0("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/GO/plot/VolcanoPlots_GOgene/", plot_file_name), volcano_plot, width = 12, height = 10)
        print(paste("Saved volcano plot:", plot_file_name))
    }}
