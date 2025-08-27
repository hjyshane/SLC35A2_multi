# read this function with 
# source('path/to/this/file/10_atac_rna_merge.R')

#' Merge processed RNA and ATAC Seurat objects by common cells
#'
#' @param processed_rna Seurat object with RNA data (e.g., post-SCTransform)
#' @param processed_atac Seurat object with ATAC data (e.g., post-LSI)
#' @param save Logical. Whether to save the combined object
#' @param save_dir Directory to save the combined object if `save = TRUE`
#'
#' @return A Seurat object with RNA and ATAC assays combined.
#' @export

atac_rna_merge <- function(
    processed_rna,
    processed_atac,
    save = TRUE,
    save_dir) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  
  # Find common cell find
  common_cells <- intersect(
    colnames(processed_rna),
    colnames(processed_atac))
  if (length(common_cells) == 0) {
    stop("No overlapping cells found between RNA and ATAC objects.")
  }
  
  # Get RNA and ATAC data from common cells
  processed_rna <- processed_rna[, common_cells]
  processed_atac <- processed_atac[, common_cells]
  
  # Create copy of RNA object
  combined_obj <- processed_rna
  
  # Create ATAC data in combined object
  combined_obj[["ATAC"]] <- processed_atac[["ATAC"]]
  
  # Extract metadata and add missing meta data columns from ATAC
  # Get  meta data of each object 
  atac_metadata <- processed_atac@meta.data
  rna_metadata <- processed_rna@meta.data
  
  # Identify ATAC-specific meta data
  # We crated combined_object based on RNA so it does not have ATAC meta data
  atac_specific <- setdiff(
    colnames(atac_metadata),
    colnames(rna_metadata))
  
  # Add ATAC specific metadata to combined object
  for (missing in atac_specific) {
    combined_obj[[missing]] <- atac_metadata[[missing]]
  }
  
  # Save object
  if (save) {qs::qsave(combined_obj, file = file.path(save_dir, "combined_obj.qs"))
    message('Saved combined object to: ', file.path(save_dir, "combined_obj.qs"), ".")
    }
  return(combined_obj)
  }
