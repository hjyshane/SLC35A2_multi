# read this function with 
# source('path/to/this/file/09_process_atac.R')

#' Process ATAC Seurat object: TF-IDF, LSI, clustering, markers, UMAP
#'
#' @param filtered_atac Filtered Seurat object (ATAC).
#' @param graph_name Vector of graph names (default: c("ATAC_nn", "ATAC_snn")).
#' @param umap_name Name for UMAP embedding (default: "ATAC_umap").
#' @param cluster_name Prefix for clustering result columns.
#' @param range_start Starting resolution (default: 0.4).
#' @param range_end Ending resolution (default: 1.2).
#' @param range_step Step size for resolution grid (default: 0.1).
#' @param save Logical. Whether to save outputs.
#' @param qsave_dir Directory to save .qs files (required if save = TRUE).
#' @param markers_dir Directory to save marker CSVs (required if save = TRUE).
#'
#' @return Processed Seurat object with dimensional reductions and clustering.
#' @export
process_atac <- function(
    filtered_atac,
    range_start = 0.4,
    range_end = 1.2,
    range_step = 0.1,
    save = TRUE,
    qsave_dir = NULL,
    markers_dir = NULL
    ) {
  # Save check
  if (save) {if (is.null(qsave_dir)) {stop("You must provide 'qsave_dir' when save = TRUE.")}
    if (is.null(markers_dir)) {stop("You must provide 'markers_dir' when save = TRUE.")}
    if (!dir.exists(qsave_dir)) {dir.create(qsave_dir, recursive = TRUE)}
    if (!dir.exists(markers_dir)) {dir.create(markers_dir, recursive = TRUE)}
  }
  
  
  # Core ATAC processing
  Seurat::DefaultAssay(filtered_atac) <- "ATAC"
  
  # Raw ATAC data is transformed using RunTFIDF and FindTopFeatures.
  #	Dimensionality reduction is performed using RunSVD, creating the lsi reduction.
  #	UMAP, clustering, and neighbors rely on the lsi reduction.
  
  processed_atac <- filtered_atac %>%
    # Normalizes ATAC-seq data using Term Frequency-Inverse Document Frequency (TF-IDF) - Dimensionality Reduction, which adjusts for the varying accessibility of regions across cells.
    # This step is crucial for focusing on biologically meaningful peaks.
    Signac::RunTFIDF(verbose = T) %>%
    # Selects the most variable peaks for downstream analysis, improving computational efficiency and accuracy.
    Signac::FindTopFeatures(
      min.cutoff = 10,
      verbose = T) %>%
    # Performs Latent Semantic Indexing (LSI), a dimensionality reduction technique tailored for sparse data like ATAC-seq.
    Signac::RunSVD(verbose = T) %>%
    # Builds k-nearest neighbor (kNN) and shared nearest neighbor (SNN) graphs based on LSI dimensions. The ATAC_snn graph is used for clustering.
    Seurat::FindNeighbors(
      reduction = "lsi", 
      dims = 2:30, 
      graph.name = c("ATAC_nn", "ATAC_snn"),
      verbose = T) %>%
    # Computes a UMAP embedding for visualization of clusters.Uses the LSI dimensions for input..
    Seurat::RunUMAP(
      reduction.name = 'ATAC_umap',
      reduction = 'lsi', 
      dims = 1:30,
      verbose = T)
  
  # Clustering
  cluster_range = seq(range_start, range_end, by = range_step)
  for (i in cluster_range) {
    processed_atac <- Seurat::FindClusters(
      processed_atac,
      graph.name = 'ATAC_snn',
      resolution = i,
      cluster.name = paste0('ATAC_cluster_', i))
  }
  
  # Find markers
  markers <- list()
  for (i in cluster_range) {
    cluster <- paste0('ATAC_cluster_', i)
    Seurat::Idents(processed_atac) <- processed_atac[[cluster]][,1]
    markers[[cluster]] <- Seurat::FindAllMarkers(
      processed_atac,
      assay = 'ATAC',
      slot = 'data',
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25)
  }
  # Save object
  if (save) {qs::qsave(processed_atac, file = file.path(qsave_dir, "processed_atac.qs"))
    message(' Saved processed ATAC object to: ', file.path(qsave_dir, "processed_atac.qs"))
    
    for (name in names(markers)) {
      write.csv(markers[[name]], file = file.path(markers_dir, paste0(name, "_markers.csv")))
      message('Saved processed markers to ', file.path(markers_dir, paste0(name, "_markers.csv")))
    }
  }
  return(processed_atac)
  }