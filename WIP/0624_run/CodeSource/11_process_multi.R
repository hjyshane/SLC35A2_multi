# read this function with 
# source('path/to/this/file/11_process_multi.R')

#' Process multimodal Seurat object with RNA + ATAC
#'
#' @param combined_obj Seurat object with RNA + ATAC
#' @param range_start Start of resolution range
#' @param range_end End of resolution range
#' @param range_step Step size for clustering resolutions
#' @param save Logical. Whether to save the annotated object.
#' @param save_dir Directory to save the annotated object if `save = TRUE`.
#' @param qsave_dir Path to save .qs file
#' @param markers_dir Path to save markers
#'
#' @return Processed multimodal Seurat object
#' @export
process_multi <- function(
    combined_obj,
    range_start = 0.1,
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
  cluster_range = seq(range_start, range_end, by = range_step)
  
  # Re-process ATAC
  Seurat::DefaultAssay(combined_obj) <- "ATAC"
  combined_obj <- combined_obj %>%
    # Normalize raw ATAC data
    Signac::RunTFIDF(verbose = T) %>%  
    # Select variable peaks
    Signac::FindTopFeatures(
      min.cutoff = 10,
      verbose = T) %>% 
    # Dimensionality reduction (creates "lsi")
    Signac::RunSVD(verbose = T) %>%
    Seurat::FindNeighbors(
      reduction = 'lsi', 
      dims = 2:30,
      graph.name =  c("ATAC_int_nn", "ATAC_int_snn"),
      verbose = T) %>%
    Seurat::RunUMAP(
      assay = 'ATAC',
      dims = 2:30,
      reduction = 'lsi',
      reduction.name = 'ATAC_int_umap',
      verbose = T)
  
  # Re-process RNA
  Seurat::DefaultAssay(combined_obj) <- "RNA"
  combined_obj <- combined_obj %>%
    Seurat::SCTransform(
      method = 'glmGamPoi',
      assay = 'RNA',
      vars.to.regress = c('percent_mt_rna'),
      vst.flavor = 'v2',
      verbose = T) %>%
    Seurat::RunPCA(
      assay = 'SCT',
      verbose = T) %>%
    Seurat::FindNeighbors(
      dims = 1:30,
      reduction = 'pca', 
      graph.name = c("RNA_int_nn", "RNA_int_snn"),
      verbose = T) %>%
    Seurat::RunUMAP(
      assay = 'SCT',
      dims = 1:30,
      reduction = 'pca',
      reduction.name = 'RNA_int_umap',
      verbose = T)
  
  # Multimodal neighbor
  # build a joint neighbor graph using both assays
  # The FindMultiModalNeighbors step focuses on integrating RNA and ATAC modalities into a shared neighbor graph but does not modify the residual variance or scaled data in the "SCT" assay.
  # integrates RNA (PCA reduction) and ATAC (LSI reduction) data into a shared nearest-neighbor graph.
  processed_comb <- combined_obj %>%
    Seurat::FindMultiModalNeighbors(
      reduction.list = list('pca', 'lsi'), # reduction.list: Specifies RNA and ATAC reductions ("pca" for RNA, "lsi" for ATAC).
      dims.list = list(1:30, 2:30), # dims.list: Defines the number of dimensions to use for each modality.
      knn.graph.name = "Multi_wknn", # Multimodal_wknn: Weighted k-nearest neighbor graph.
      snn.graph.name = "Multi_wsnn", # Multimodal_wsnn: Shared nearest neighbor (SNN) graph.
      weighted.nn.name = "Multi_w.nn", # Multimodal_weighted.nn: Final multi-modal neighbor graph.
      modality.weight.name = c('Multi_RNA.w', 'Multi_ATAC.w'), # modality.weight.name: Adds weights to indicate the relative contribution of RNA and ATAC.
      verbose = T) %>%
    # build a joint UMAP visualization
    # Creates a UMAP embedding based on the joint multi-modal neighbor graph ("Multimodal_weighted.nn").
    # nn.name: Uses the joint neighbor graph.
    # UMAP embedding for multi-modal neighbors does not rely on the data in the SCT assay but instead leverages the Multimodal_weighted.nn graph. 
    # Hence, the choice of assay (SCT or ATAC) does not impact the multi-modal UMAP embedding.
    Seurat::RunUMAP(
      nn.name = "Multi_w.nn",
      reduction.name = "Multi_UMAP",
      verbose = TRUE)
  
  # Clusters for RNA and ATAC
  for (i in cluster_range) {
    # RNA cluster
    processed_comb <- Seurat::FindClusters(
      processed_comb, 
      graph.name = "RNA_int_nn",
      resolution = i,
      cluster.name = paste0('RNA_int_cluster_', i),
      verbose = T)
    
    # ATAC cluster
    processed_comb <- Seurat::FindClusters(
      processed_comb,
      graph.name = "ATAC_int_nn",
      resolution = i,
      cluster.name = paste0('ATAC_int_cluster_', i),
      verbose = T)
    
    # Multimodal cluster
    processed_comb <- Seurat::FindClusters(
      processed_comb,
      graph.name = "Multi_wsnn",
      resolution = i,
      cluster.name = paste0('Multi_cluster_', i),
      verbose = T)
  }
  
  # Find markers
  processed_comb <- PrepSCTFindMarkers(processed_comb)
  markers <- list()
  for (i in cluster_range) {
    Seurat::Idents(processed_comb) <- processed_comb[[paste0('ATAC_int_cluster_', i)]][,1]
    markers[[paste0('ATAC_int_cluster_', i)]] <- Seurat::FindAllMarkers(
      processed_comb,
      assay = "ATAC",
      slot = 'data',
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = T)
    
    Seurat::Idents(processed_comb) <- processed_comb[[paste0('RNA_int_cluster_', i)]][,1]
    markers[[paste0('RNA_int_cluster_', i)]] <- Seurat::FindAllMarkers(
      processed_comb,
      assay = "SCT",
      slot = 'data',
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = T)   
    
    Seurat::Idents(processed_comb) <- processed_comb[[paste0('Multi_cluster_', i)]][,1]
    markers[[paste0('Multi_cluster_', i)]] <- Seurat::FindAllMarkers(
      processed_comb,
      assay = "SCT",
      slot = 'data',
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = T)
  }
  
  # Save object
  if (save) {qs::qsave(processed_comb, file = file.path(qsave_dir, "processed_comb.qs"))
    message('Saved combined object to: ', file.path(qsave_dir, "processed_comb.qs"), ".")
    for (name in names(markers)) {
      write.csv(markers[[name]], file = file.path(markers_dir, paste0(name, "_markers.csv")))
      message('Saved processed markers to ', file.path(markers_dir, paste0(name, "_markers.csv")))
    }
    
    
  }
  return(processed_comb)
}