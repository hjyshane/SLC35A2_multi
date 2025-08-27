# read this function with 
# source('path/to/this/file/13_ref_mousebrain_annotation.R')

#' Annotate query object using mouse brain reference with SCT label transfer
#'
#' @param ref_obj Reference Seurat object with known annotations.
#' @param sobj Query Seurat object to be annotated.
#' @param save Logical. Whether to save the annotated object.
#' @param save_dir Directory to save the annotated object if `save = TRUE`.
#'
#' @return Annotated query object with predicted cell types.
#' @export
ref_mousebrain_annotation <- function(
    ref_obj,
    sobj,
    save = TRUE,
    save_dir = NULL) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  
  Seurat::DefaultAssay(sobj) <- "SCT"
  
  # Find anchors between sobj and ref_obj
  anchors <- Seurat::FindTransferAnchors(
    reference = ref_obj,
    query = sobj,
    dims = 1:30,
    normalization.method = 'SCT')
  
  # Predict cell types
  predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_obj$Description,   # or whichever column contains cell type labels
    dims = 1:30
  )
  
  # Add metadata
  sobj <- Seurat::AddMetaData(
    object = sobj,
    metadata = predictions)
  
  # Save object
  if (save) {qs::qsave(sobj, file = file.path(save_dir, "ref_annotated_obj.qs"))
    message('Saved reference annotated object to: ', file.path(save_dir, "ref_annotated_obj.qs"), ".")
  }
  return(sobj)
  }
  