# read this function with
# source('path/to/this/file/16_find_dg.R')

#' Run differential expression analysis for specific cell types and group comparisons
#'
#' This function performs DEG analysis using either Wilcox or MAST test,
#' generates volcano plots, and saves background and significant results per comparison.
#'
#' @param seurat_obj Seurat object containing the data.
#' @param cell_type_meta Column name in `meta.data` that holds cell type labels (string).
#' @param comparison_meta Column name in `meta.data` that defines grouping for comparison (string).
#' @param cell_type Cell type to subset and analyze.
#' @param comparisons A list of named comparisons (e.g., list(PTZvsSAL = list(name = "PTZvsSAL", group1 = "PTZ", group2 = "SAL"))).
#' @param minimum_cnt Minimum number of cells required per group (default = 3).
#' @param mode Test method: "wilcox", "MAST", or "LR"(for ATAC) (default = "wilcox").
#' @param p_value Adjusted p-value threshold for significance (default = 0.05).
#' @param fc_value Minimum absolute log2 fold change to call a gene significant (default = 0.2).
#' @param save Logical; whether to save results (default = TRUE).
#' @param bg_dir Directory to save full background DEG tables.
#' @param sig_dir Directory to save filtered significant DEG tables.
#' @param plot_dir Directory to save volcano plots.
#'
#' @return None. Results are optionally saved to CSVs and PNGs.
#' @export

run_dg <- function(
    seurat_obj,
    cell_type_meta,
    comparison_meta,
    cell_types,
    comparisons,
    minimum_cnt = 10,
    mode = "MAST", 
    p_value = 0.05,
    fc_value = 0.2,
    save = TRUE,
    bg_dir = NULL,
    sig_dir = NULL,
    plot_dir = NULL) {
  # Save check
  if (save) {if (is.null(bg_dir)) {stop("You must provide 'bg_dir' when save = TRUE.")}
    if (!dir.exists(bg_dir)) {dir.create(bg_dir, recursive = TRUE)}
  }
  if (save) {if (is.null(sig_dir)) {stop("You must provide 'sig_dir' when save = TRUE.")}
    if (!dir.exists(sig_dir)) {dir.create(sig_dir, recursive = TRUE)}
  }
  if (save) {if (is.null(plot_dir)) {stop("You must provide 'plot_dir' when save = TRUE.")}
    if (!dir.exists(plot_dir)) {dir.create(plot_dir, recursive = TRUE)}
  }
  for (ct in cell_types) {
    # For each cell type, loop through each comparison group
    for (comparison in comparisons) {
      # Subset for cell type
      
      # Set idents
      Seurat::Idents(seurat_obj) <- cell_type_meta
      
      # subsetting cell type in each condition
      cell_type_data <- subset(seurat_obj,
                               subset = cell_type == ct & Sample %in% c(comparison$group1, comparison$group2)) # cell_type is metadata column name for cell type, Sample is for group
      
      if (length(cell_type_meta >= minimum_cnt)) {
      # Set identities
      Seurat::DefaultAssay(cell_type_data) <- "RNA"
      Seurat::Idents(cell_type_data) <- "Sample"
      
      # if (mode != "LR") {
      # DE analysis
      bg_results <- Seurat::FindMarkers(
        object = cell_type_data,
        assay = "RNA",
        ident.1 = comparison$group1,
        ident.2 = comparison$group2,
        logfc.threshold = 0.0,
        min.pct = 0.25,
        test.use = mode,
        min.cells.group = 10,
        latent.vars = "nFeature_RNA"
      )
      
      # Add metadata info
      bg_results$gene <- rownames(bg_results)
      bg_results$cell_type <- ct
      bg_results$comparison <- comparison$name
      
      # Add regulation info
      bg_results <- bg_results %>%
        mutate(Regulation = case_when(
          avg_log2FC > 0 ~ "Up",
          avg_log2FC < 0 ~ "Down",
          TRUE ~ "NoChange"))

      bg_results <- bg_results %>%
        mutate(Viz = case_when(
          abs(avg_log2FC) >= fc_value & p_value <= 0.05  ~ "Sig",
          abs(avg_log2FC) >= fc_value & p_value > 0.05 ~ "log2FC",
          abs(avg_log2FC) < fc_value & p_value <= 0.05 ~ "pValue",
          TRUE ~ "NoChange"))
      
      # Get significant genes - DEGs
      sig_results <- bg_results %>%
        filter(p_val < p_value & abs(avg_log2FC) >= fc_value)

      # Create VolcanoPlot
      if (nrow(bg_results) == 0) {
        message("There are no genes/peaks to plot")
        next
      } else {
         vplot1 <- EnhancedVolcano::EnhancedVolcano(
          bg_results,
          lab = bg_results$gene,
          x = 'avg_log2FC',
          y = 'p_val_adj',       
          title = paste0('MAST (combined p): ', comparison$group1, ' vs ', comparison$group2),
          subtitle = 'latent.vars = nFeature_RNA; RNA log-normalized',
          pCutoff = p_value,           
          FCcutoff = fc_value,         
          xlab = bquote(~Log[2]~ 'FC (Seurat/MAST)'),
          ylab = bquote(~-Log[10]~italic(P)),
          pointSize = 1.0,
          labSize = 3.0,
          legendPosition = 'right',
          drawConnectors = F,
          max.overlaps = 20)
         
       
         
         vplot2 <- ggplot2::ggplot(bg_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
           ggplot2::geom_point(aes(alpha = pmin(1, (pct.1 + pct.2)/2), color = Viz, size = 1)) +
           scale_color_manual(values = c("Sig" = "red", "log2FC" = "green", "pValue" = "blue", "NoChange" = "grey")) +
           ggplot2::geom_vline(xintercept = c(-fc_value, fc_value), linetype = 2) +
           ggplot2::geom_hline(yintercept = -log10(p_value), linetype = 2) +
           labs(x = "avg_log2FC", y = "-log10(p_val_adj)",
                title = paste0("MAST (combined p) : ", comparison$group1, " vs ", comparison$group2),
                subtitle = "latent.vars = nFeature_RNA; min.pct = 0.25, min.cells.group  = 10") +
           theme_minimal()
      }
      # Save
      if (save) {
        if (nrow(bg_results) == 0) {
          message("There are no genes/peaks to found")
          next
        } else {
          write.csv(bg_results, file = file.path(bg_dir, paste0(ct, "_", comparison$name, ".csv")))
          ggsave(plot = vplot1, file = file.path(plot_dir, paste0(ct, "_", comparison$name, "_EV.png")), width = 10, height = 12, dpi = 300)
          ggsave(plot = vplot2, file = file.path(plot_dir, paste0(ct, "_", comparison$name, "_gg.png")), width = 10, height = 12, dpi = 300)
          
                  }
        if (nrow(sig_results) == 0) {
          message("There are no deg genes/peaks to found")
          next
        } else {
          write.csv(sig_results, file = file.path(sig_dir, paste0(ct, "_", comparison$name, ".csv")))
        }
        
      }
      } else next
    }
  }
}

  

