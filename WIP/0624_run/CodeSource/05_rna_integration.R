# read this function with 
# source('path/to/this/file/05_rna_integration.R')

#' Integrate SCTransformed Seurat objects
#'
#' This function performs SCTransform-based integration of a list of Seurat objects.
#'
#' @param sobj_list A named list of SCTransformed Seurat objects.
#' @param save Logical. Whether to save the integrated object (default: FALSE).
#' @param save_dir Directory to save the integrated object if `save = TRUE`.
#'
#' @return Integrated Seurat object.
#' @export
rna_integration <- function(
    sobj_list,
    save = TRUE,
    save_dir = NULL) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide save_dir when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  
  # Select features for integration
  # Identifies common features (genes) across datasets that are informative for integration. 
  # These features are used for anchor identification and subsequent data integration.
  features <- Seurat::SelectIntegrationFeatures(
    object.list = sobj_list,
    nfeatures = 3000)

  # Prepare for integration
  # Prepares the datasets for integration by ensuring that the same SCT model is applied to the selected features across datasets.
  sobj_list <- Seurat::PrepSCTIntegration(
    object.list = sobj_list,
    anchor.features = features)
  
  # Find integration anchors
  # Identifies anchors (shared cells or cell populations) across datasets based on the selected features.
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = sobj_list,
    normalization.method = 'SCT',
    anchor.features = features)

    # Find integration anchors
  # Identifies anchors (shared cells or cell populations) across datasets based on the selected features.
  integrated_rna <- Seurat::IntegrateData(
    anchorset = anchors,
    normalization.method = 'SCT')
  
  # Save object
  if (save) {qs::qsave(integrated_rna, file = file.path(save_dir, "integrated_rna.qs"))
  message('Saved integraed rna objec to ', file.path(save_dir, "integrated_rna.qs"))
  }
  return(integrated_rna)
  }