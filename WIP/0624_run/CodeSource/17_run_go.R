# read this function with 
# source('path/to/this/file/17_run_go.R')


#' Run GO enrichment analysis with optional visualization (dotplot, cnetplot, RRvGO)
#'
#' This function loops over CSVs of significant genes (e.g., from DEG analysis), 
#' runs GO enrichment using `clusterProfiler`, and optionally saves:
#' - GO term tables
#' - dotplots (top terms)
#' - cnetplots
#' - RRvGO similarity plots
#'
#' @param seurat_obj A Seurat object (used for context, not directly used here).
#' @param input_sig Directory path containing significant gene CSV files.
#' @param orgdb Organism annotation package (default: org.Mm.eg.db).
#' @param ontology GO ontology: "BP", "MF", or "CC" (default: "BP").
#' @param bg_gene Vector of background gene symbols.
#' @param top_n Number of top GO terms to include in plots (default = 20).
#' @param save Logical. Whether to save results and plots (default: TRUE).
#' @param cnet Logical. Whether to generate cnetplots (default: FALSE).
#' @param rrvgo Logical. Whether to run RRvGO reduction (default: FALSE).
#' @param go_dir Directory to save GO term tables and dotplots.
#' @param cnet_dir Directory to save cnetplots (required if cnet = TRUE).
#' @param rrvgo_dir Directory to save RRvGO plots and CSVs (required if rrvgo = TRUE).
#'
#' @return None. Saves outputs to specified directories.
#' @export
run_go <- function(
    seurat_obj,
    input_sig,
    orgdb = org.Mm.eg.db,
    ontology = "BP",
    bg_gene,
    top_n = 20,
    save = TRUE,
    cnet = FALSE,
    rrvgo = FALSE,
    go_dir = NULL,
    cnet_dir = NULL,
    rrvgo_dir = NULL) {
  # Save check
  if (save && is.null(go_dir)) stop("You must provide 'go_dir' when save = TRUE.")
  if (cnet && is.null(cnet_dir)) stop("You must provide 'cnet_dir' when cnet = TRUE.")
  if (rrvgo && is.null(rrvgo_dir)) stop("You must provide 'rrvgo_dir' when rrvgo = TRUE.")
  if (save) dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  if (cnet) dir.create(cnet_dir, recursive = TRUE, showWarnings = FALSE)
  if (rrvgo) dir.create(rrvgo_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sig_file in list.files(path = input_sig, 
                             pattern = "\\.csv", 
                             full.names = TRUE)) {
    
    # Extract metadata from filename
    file_name <- basename(sig_file)
    parts <- strsplit(gsub(".csv", "", file_name), "_")[[1]]
    
    if (length(parts) == 4) {
      cell_type <-paste(parts[1:2], collapse = "_") 
      comparison <- paste(parts[3:4], collapse = "_")
    } else if (length(parts) == 3) {
      cell_type <-paste(parts[1], collapse = "_") 
      comparison <- paste(parts[-1], collapse = "_")
    } else {message("check file name")}
    
  
    # Read files
    sig_deg <- read_csv(sig_file)
    
    # Get gene
    sig_genes <- sig_deg$gene
    
    # Run GO analysis
    go_result <- clusterProfiler::enrichGO(
      gene = sig_genes,
      OrgDb = orgdb,
      keyType = "SYMBOL",
      ont = ontology,
      universe = bg_gene,
      pAdjustMethod = "bonferroni",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05)
    
    # Check if there is GO term found
    if (is.null(go_result) || 
        !("result" %in% slotNames(go_result)) || 
        nrow(go_result@result) == 0) {
      message("No GO terms found for ", sig_file)
      next
    } 
    
    # Get dataframe
    go_df <- as.data.frame(go_result)
    
    # Skip if go_df is empty
    if (nrow(go_df) == 0) {
      message("Empty GO dataframe for ", sig_file)
      next
    }
    
    # Compute additional matrics
    go_df <- go_df %>%
      mutate(
        log_p_adj = -log10(p.adjust),
        GeneRatio_numeric = sapply(strsplit(GeneRatio, '/'), function(x) as.numeric(x[1])/as.numeric(x[2])),
        BgRatio_numeric = sapply(strsplit(BgRatio, '/'), function(x) as.numeric(x[1])/as.numeric(x[2])),
        Fold_enrichment = GeneRatio_numeric / BgRatio_numeric,
        score = Fold_enrichment * log_p_adj,
        cell_type = cell_type,
        comparison = comparison)
    
    # Get top n terms
    top_terms <- go_df %>%
      arrange(p.adjust) %>%
      head(top_n)
    
    # DotPlot
    go_plot <- 
      ggplot(
        top_terms,
        aes(x = reorder(Description, log_p_adj),
            y = log_p_adj,
            fill = Fold_enrichment,
            size = Count)) +
      geom_point(shape = 21) +
      coord_flip() +
      scale_fill_viridis_c(option = 'magma') +
      theme_minimal() +
      labs(
        title = paste("Top go terms ", ontology, "_", cell_type, "_", comparison),
        x = "GO Term",
        y = "-log10(p-value)",
        fill = "Fold enrichment",
        size = "Gene count")
    
    # Save
    if (save) {
      dir.create(file.path(go_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(go_dir, "plot"), recursive = TRUE, showWarnings = FALSE)
      write.csv(go_df, file = file.path(go_dir, "csv",  paste0("GO_", ontology, "_", cell_type, "_",  comparison, ".csv")))
      ggsave(filename = file.path(go_dir, "plot", paste0("GOPlot_", ontology, "_",cell_type, "_",  comparison, ".png")), plot = go_plot, width = 12, height = 12, dpi = 300)
    }
    # CNET plot
    if (cnet) {
      # check to see if we have enough terms to make cnetplot
      if (nrow(top_terms) > 2) {
        cnet_plot <- enrichplot::cnetplot(
          go_result,
          showCategory = max(5, nrow(top_terms)),
          circular = TRUE,
          node_label = 'gene') +
          ggtitle(paste("GO cnetplot ", ontology, "_", cell_type, "_", comparison))
        
        ggsave(filename = file.path(cnet_dir, paste0("cnetplot_", ontology, "_", cell_type, "_",  comparison, ".png")), plot = cnet_plot, width = 12, height = 12, dpi = 300)
      } else {
        message("Not enough terms for CnetPlot ",ontology, "_", cell_type, '_', comparison)
      }
    }
    # rrvgo similarity plot
    if (rrvgo) {
      # check to see if we have enough terms to make rrvgo calculation
      if (nrow(go_df) > 1 && length(unique(go_df$ID)) >1) {
        similarity <- rrvgo::calculateSimMatrix(
          go_df$ID,
          orgdb = orgdb,
          ont = ontology,
          method = 'Rel')
        
        # get reduced terms
        if (!is.null(similarity) && length(similarity) > 0) {
          reduced <- rrvgo::reduceSimMatrix(
            similarity,
            scores = setNames(go_df$log_p_adj, go_df$ID),
            threshold = 0.5,
            orgdb = orgdb)
          
          # save
          if (!is.null(reduced) && length(reduced) > 0) {
            write.csv(reduced, file = file.path(rrvgo_dir, paste('rrvgo_', ontology, "_",cell_type, "_", comparison, '.csv')))
            png(file = file.path(rrvgo_dir, paste('rrvgo_',ontology, "_", cell_type, "_", comparison, '.png')), width = 12, height = 12, units = 'in', res = 300)
            dev.off()
              }
        } else {
          message("Not enough terms for rrvgo calculation ",ontology, "_", cell_type, '_', comparison)
        }
      }
    } 
  }
}

