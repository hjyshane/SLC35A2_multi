# read this function with 
# source('path/to/this/file/06_process_integrated.R')

#' Process integrated Seurat object (PCA, clustering, markers, UMAP/TSNE)
#'
#' @param integrated_obj Integrated Seurat object
#' @param assay Name of assay to use (default: "integrated")
#' @param graph_name Character vector with kNN and SNN graph names (default: c("RNA_nn", "RNA_snn"))
#' @param cluster_name Prefix for clustering results
#' @param umap_name Name to assign to UMAP embedding
#' @param range_start Start of resolution range for clustering
#' @param range_end End of resolution range for clustering
#' @param range_step Step size for resolution
#' @param save Logical, whether to save outputs
#' @param qsave_dir Directory to save processed object
#' @param markers_dir Directory to save marker CSVs
#' @return Processed Seurat object
#' @export
process_integrated <- function(
    integrated_rna,
    range_start = 0.4,
    range_end = 1.2,
    range_step = 0.1,
    save = TRUE,
    qsave_dir = NULL,
    markers_dir = NULL) {
  # Save check
  if (save) {if (is.null(qsave_dir)) {stop("You must provide 'qsave_dir' when save = TRUE.")}
    if (is.null(markers_dir)) {stop("You must provide 'markers_dir' when save = TRUE.")}
    if (!dir.exists(qsave_dir)) {dir.create(qsave_dir, recursive = TRUE)}
    if (!dir.exists(markers_dir)) {dir.create(markers_dir, recursive = TRUE)}
  }

  # Process
  Seurat::DefaultAssay(integrated_rna) <- 'integrated'
  
  integrated_rna <- integrated_rna %>%
    # PCA
    Seurat::RunPCA(
      assay = 'integrated', 
      verbose = T) %>% 
    #	Identify cell neighborhoods and perform clustering: 
    #	Compute the k-nearest neighbor graph based on the PCA results.
    #	kNN graph (RNA_nn): Represents direct neighbors of each cell.	A k-nearest neighbor graph 
    #	SNN graph (RNA_snn): Shared nearest neighbor graph, used for clustering. derived from the kNN graph
    Seurat::FindNeighbors(
      dims = 1:30,
      reduction = 'pca',
      graph.name = c('RNA_nn', 'RNA_snn'),
      verbose = T) %>%
    # Dimensional recudtion for visual
    Seurat::RunUMAP(
      reduction = 'pca',
      reduction.name = 'RNA_umap',
      dims = 1:30,
      verbose = T) %>%
    Seurat::RunTSNE(
      reduction = 'pca',
      dims = 1:30,
      verbose = T) %>%
    # When SCTransform is run, residuals of the gene expression are stored in the "scale.data" slot of the "SCT" assay.
    # PrepSCTFindMarkers ensures that the residual variance of the genes is appropriately scaled and stabilized before marker detection.
    # Run PrepSCTFindMarkers after clustering and before running FindAllMarkers.
    Seurat::PrepSCTFindMarkers(
      verbose = T)
  
  Seurat::DefaultAssay(integrated_rna) <- 'SCT'
  
  # Clustering
  cluster_range = seq(range_start, range_end, by = range_step)
  for (i in cluster_range) {
    integrated_rna <- Seurat::FindClusters(
      integrated_rna,
      graph.name = 'RNA_snn',
      resolution = i,
      cluster.name = paste0("RNA_cluster_", i))
  }
  
  # Find markers
  # The "SCT" assay contains normalized and variance-stabilized data suitable for identifying differentially expressed genes between clusters.
  # The "integrated" assay is batch-corrected for alignment and not designed for marker detection as it may suppress biological variation in favor of alignment.
  # only.pos = TRUE: Returns only positive markers (regions more accessible in the cluster).
  # min.pct: Minimum fraction of cells expressing the feature in the cluster.
  # logfc.threshold: Minimum log fold-change required to call a feature significant.
  markers <- list()
  for (i in cluster_range) {
    cluster <- paste0('RNA_cluster_', i)
    Seurat::Idents(integrated_rna) <- integrated_rna[[cluster]][,1]
    markers[[cluster]] <- Seurat::FindAllMarkers(
      integrated_rna,
      assay = 'SCT',
      slot = 'data',
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25)
  }
  
  preccesed_rna <- integrated_rna
  
  # Save object
  if (save) {qs::qsave(preccesed_rna, file = file.path(qsave_dir, "processed_rna.qs"))
  message('Saved processed integraed rna objec to ', file.path(qsave_dir, "processed_rna.qs"))
  
    for (name in names(markers)) {
    write.csv(markers[[name]], file = file.path(markers_dir, paste0(name, "_markers.csv")))
    message('Saved processed markers to ', file.path(markers_dir, paste0(name, "_markers.csv")))
    }
  }
  return(preccesed_rna)
  }
