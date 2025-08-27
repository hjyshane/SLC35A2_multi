# Define DEG directory
deg_dir <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/DEG/sig_csv/"
deg_files <- list.files(deg_dir, pattern = ".csv", full.names = TRUE)

# Define output base
output_dir_base <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/fGSEA"

# Create output dirs
for (cat in c("fGSEA_GOMF", "fGSEA_GOBP", "fGSEA_KEGG", "fGSEA_TFT", "fGSEA_REACTOME")) {
    dir.create(file.path(output_dir_base, cat, "csv"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir_base, cat, "plot"), recursive = TRUE, showWarnings = FALSE)
}

# Load pathway databases
pathways_MF <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF")
pathways_BP <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
pathways_KG <- msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
pathways_RE <- msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
# pathways_TFT <- msigdbr::msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT")

# Format into named lists
list_MF <- split(pathways_MF$gene_symbol, pathways_MF$gs_name)
list_BP <- split(pathways_BP$gene_symbol, pathways_BP$gs_name)
list_KG <- split(pathways_KG$gene_symbol, pathways_KG$gs_name)
list_RE <- split(pathways_RE$gene_symbol, pathways_RE$gs_name)
# list_TFT <- split(pathways_TFT$gene_symbol, pathways_TFT$gs_name)

# Loop through DEG files
for (file in deg_files) {
    name <- gsub(".csv", "", basename(file))
    df <- read_csv(file)

    # Validate column existence
    if (!all(c("gene", "avg_log2FC") %in% colnames(df))) {
        warning(paste("Skipping", name, "- missing gene or avg_log2FC columns"))
        next
    }

    df <- df %>%
        arrange(desc(avg_log2FC)) %>%
        dplyr::select(gene, avg_log2FC)

    if (nrow(df) < 5) {
        message(sprintf("Skipping %s: only %d genes (minimum 5 required)", name, nrow(df)))
        next
    }

    ranked <- setNames(df$avg_log2FC, df$gene)

    # Run fGSEA for each category
    results <- list(
        GOBP = fgsea::fgseaMultilevel(list_BP, ranked, minSize = 5, maxSize = 500, nPerm = 1000) %>% arrange(padj),
        GOMF = fgsea::fgseaMultilevel(list_MF, ranked, minSize = 5, maxSize = 500, nPerm = 1000) %>% arrange(padj),
        KEGG = fgsea::fgseaMultilevel(list_KG, ranked, minSize = 5, maxSize = 500, nPerm = 1000) %>% arrange(padj),
        # TFT = fgseaMultilevel(list_TFT, ranked, minSize = 5, maxSize = 500, nPerm = 1000) %>% arrange(padj),
        REACTOME =  fgsea::fgseaMultilevel(list_RE, ranked, minSize = 5, maxSize = 500, nPerm = 1000) %>% arrange(padj)
    )

    # Save CSV outputs
    for (cat in names(results)) {
        out_path <- file.path(output_dir_base, paste0("fGSEA_", cat), "csv", paste0("fGSEA_", name, ".csv"))
        data.table::fwrite(results[[cat]], file = out_path)
    }

    # Plot significant results (padj < 0.05) – top 10 only
    for (cat in names(results)) {
        top <- results[[cat]] %>%
            filter(padj < 0.05) %>%
            head(10)

        if (nrow(top) == 0) {
            message(sprintf("No significant %s pathways for %s", cat, name))
            next
        }

        # Clean up names
        top$pathway <- top$pathway %>%
            str_remove("^GOBP_|^GOMF_|^KEGG_|^REACTOME|") %>%
            str_replace_all("_", " ")
        top$pathway <- factor(top$pathway, levels = top$pathway[order(top$NES)])

        # Plot
        p <- ggplot(top, aes(x = NES, y = pathway, size = size, color = -log10(padj))) +
            geom_point() +
            scale_color_gradient(low = "#0513d8", high = "#f517d3") +
            labs(title = paste("fGSEA", cat, "Enrichment -", name),
                 x = "Normalized Enrichment Score",
                 y = "Pathway",
                 size = "Size",
                 color = "padj") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 9))

        # Save plot
        plot_file <- file.path(output_dir_base, paste0("fGSEA_", cat), "plot", paste0("fGSEA_", name, ".png"))
        ggsave(plot_file, p, width = 12, height = 10, dpi = 300)
    }
}
