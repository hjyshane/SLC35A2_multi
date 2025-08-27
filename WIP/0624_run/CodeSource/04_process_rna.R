# read this function with 
# source('path/to/this/file/04_process_rna.R')


#' Run SCTransform, PCA, UMAP, and TSNE on Seurat list
#' 
#' 
#' 
#' @param sobj_list A list of Seurat objects to process.
#' @param samples Vector of sample names.
#' @param dim Number of dimensions to use for PCA/UMAP/TSNE (default: 30).
#' @param reduction_to_use Dimensionality reduction method for UMAP/TSNE (default: "pca").
#' @param reduction_name_suffix Suffix to append to UMAP/TSNE names (default: "_rna").
#' @param save Logical. If TRUE, saves processed objects.
#' @param save_dir Directory to save .qs files (if save = TRUE).
#'
#' @return A list of processed Seurat objects.
#' @export

process_rna <- function(
    sobj_list,
    samples,
    save = TRUE,
    save_dir = NULL) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
    }
  
  # Process loop
  for (i in seq_along(samples)) {
    sample <- samples[i]
    sobj <- sobj_list[[sample]]
    message("Processing sample: ", sample)
    
    # Process
    sobj <- 
      sobj %>%
      # SCTransform
      #   Adjusts for sequencing depth (library size).Provides normalized counts while maintaining biological variability.
      #   Stabilizes the variance across genes, Identifies and prioritizes highly variable genes for analysis.
      # 	Improves signal-to-noise ratio by reducing the influence of technical noise.
      # 	Allows correction for technical artifacts or batch effects through regression.
      # vars.to.regress Variables to regress out (default: "percent_mt_rna") - Regress out unwanted sources of variation 
      # method Normalization method (default: "glmGamPoi")
      # return.only.var.genes Whether to keep only variable genes for downstream analysis
      # return.only.var.genes = FALSE to retain all genes in the normalized data. DEG and GO analyses typically require information from all genes, 
      # even those not classified as highly variable. but it will increase the computation time and memory usage. default is TRUE.
      Seurat::SCTransform(
        method = 'glmGamPoi', 
        assay = 'RNA', 
        vars.to.regress = c('percent_mt_rna'),
        verbose = T) %>%
      # PCA
      # Stores PCA results in the Seurat object under PTZ_1hr[["pca"]]
      #	Loadings: Contribution of each gene to each PC.
      #	Scores: Position of each cell in the PC space.
      # Eigenvalues: Variance captured by each PC.
      Seurat::RunPCA(
        assay = 'SCT',
        vervose = T) %>%
      # Compresses the dataset into fewer dimensions (principal components, or PCs) that capture the most variance in the data.
      # Identifies patterns of variability across genes and cells.
      # n.dims Number of dimensions to use (default: 30)
      # resolution Clustering resolution (default: 0.8)
      # assay to use SCT
      # reduction.name Name for UMAP reduction (default: "umap")
      # verbose Print progress
      Seurat::RunUMAP(
        dims = 1:30, 
        reduction = 'pca', 
        assay = 'SCT', 
        reduction.name = 'rna_umap',
        verbose = T) %>%
      Seurat::RunTSNE(
        dims = 1:30, 
        reduction = 'pca', 
        assay = 'SCT',
        verbose = T)
    
    sobj_list[[sample]] <- sobj
    
    # Save
    if (save) {qs::qsave(sobj, file = file.path(save_dir, paste0(sample, "_rna_processed.qs")))
    message("Saved processed object: ", file.path(save_dir, paste0(sample, "_rna_processed.qs")))
    }
    }
  return(sobj_list)
  }
  
