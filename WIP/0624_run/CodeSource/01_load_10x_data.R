# read this function with 
# source('path/to/this/file/01_load_10x_data.R')


#' Load multiple RNA expression matrices from 10X .h5 files
#'
#' This function loads RNA count matrices (e.g., "Gene Expression") from multiple
#' 10X Genomics HDF5 files. Optionally, it can save each matrix using `qs::qsave()`.
#'
#' @param input_dir Character. Directory containing sample folders.
#' @param samples Character vector. Names of sample folders to load. Each folder must contain
#'        `outs/filtered_feature_bc_matrix.h5`.
#' @param assay Character. Assay name to extract from each H5 file (default: "Gene Expression").
#' @param save Logical. If TRUE, saves each matrix using `qs::qsave()` (default: FALSE).
#' @param save_dir Character or NULL. Output directory for saved `.qs` files (one per sample). Required if `save = TRUE`.
#'
#' @return Named list of sparse matrices (class `dgCMatrix`), one per sample.
#'
#' @examples
#' samples <- c("PTZ_1hr", "PTZ_24hr")
#' rna_list <- load_rna_from_h5(input_dir = "~/project/data", samples = samples)
#'
#' @export

load_10x_data <- 
  function(
    input_dir, 
    assay = 'Gene Expression', 
    samples, 
    save = TRUE, 
    save_dir = NULL){
    # Save check
    if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
      if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
      }
    
    # Create empty list to store returned object
    rna_list <- list()
    
    # Read each sample's h5 files and store into rna_list[[sample]].
    message("Loading ", length(samples), " samples.")
    for (sample in samples){
      # Set file path
      h5_path <- file.path(input_dir, sample, 'outs', 'filtered_feature_bc_matrix.h5')
      
      if (file.exists(h5_path)) {
        message("Loading sample: ", sample, " from ", h5_path)
        # Read h5 and Store in rna_list
        h5 <- Seurat::Read10X_h5(h5_path)[[assay]]
        rna_list[[sample]] <- h5
        
        # Save if save is enabled
        if (save){
          qs::qsave(h5, file = file.path(save_dir, paste0(sample, "_h5.qs")))
          message("Saved object: ", sample, " to ", file.path(save_dir, paste0(sample, "_h5.qs")), ".")
          }
        } else {
          message("File not found for sample: ", sample, "\nMissing in ", h5_path)
          next
        }
      }
    message("Finished loading ", length(rna_list), " of ", length(samples), " samples.")
    return(rna_list)
    }

  
