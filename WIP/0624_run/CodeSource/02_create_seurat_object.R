# read this function with 
# source('path/to/this/file/02_create_seurat_object.R')

#' Create Seurat objects from a list of RNA count matrices
#'
#' This function takes a list of sparse RNA matrices and builds Seurat objects
#' for each sample. It uses barcode-level metadata from each sample's
#' `per_barcode_metrics.csv`. Optionally, each Seurat object can be saved to disk.
#'
#' Create Seurat objects from 10x data
#' Parameters:
#' - counts: Gene expression matrix from 10X h5 file 
#' - assay: Set RNA as the default assay type
#' - min.cells: Only keep genes expressed in at least 3 cells
#' - min.features: Only keep cells with at least 200 genes expressed
#' - max.RNA - Upper threshold for RNA counts to filter doublets (default: 25000)
#' - min.RNA - Lower threshold for RNA counts to filter empty droplets (default: 1000)
#' - meta.data: Import cell metadata from the per_barcode_metrics.csv file
#' - header = TRUE indicates the CSV has column names
#' - row.names = 1 indicates the first column is used as row names
#' 
#' @param rna_list Named list of RNA expression matrices (class `dgCMatrix`).
#' @param input_dir Character. Path containing sample subfolders (e.g., `sample/outs/...`).
#' @param project Character. Project name to assign to Seurat objects.
#' @param min.cell Integer. Minimum cells a gene must be expressed in to be retained. Default: 3.
#' @param min.feature Integer. Minimum genes a cell must express to be retained. Default: 200.
#' @param samples Character vector. Names of the samples in `rna_list`.
#' @param save Logical. If TRUE, saves each Seurat object using `qs::qsave()`. Default: FALSE.
#' @param save_dir Character or NULL. Path to directory where Seurat objects will be saved as `.qs`. Required if `save = TRUE`.
#'
#' @return Named list of Seurat objects (class `Seurat`).
#'
#' @examples
#' seurat_list <- create_seurat_object(
#'   rna_list = rna_list,
#'   input_dir = "~/project/data",
#'   samples = c("PTZ_1hr", "PTZ_24hr"),
#'   save = TRUE,
#'   save_dir = "~/project/objects"
#' )
#'
#' @export
create_seurat_object <- 
  function(rna_list, 
           input_dir,
           project = 'SeuratProject',
           min.cell = 3,
           min.feature = 200,
           samples, 
           save = TRUE, 
           save_dir = NULL) {
    # Save check
    if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
      if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
    }
    
    # Create empty seurat object list to store
    sobj_list <- list()
    
    for (i in seq_along(samples)) {
      message("Creating seurat Object for ", samples[i], ".")
      sample <- samples[i]
      
      # Set meta data path and read
      meta_path <- file.path(input_dir, sample, 'outs', 'per_barcode_metrics.csv')
      meta_data <- read.csv(meta_path, header =  TRUE, row.names = 1)
      
      # Create seurat object
      sobj <- Seurat::CreateSeuratObject(
        counts = rna_list[[i]],
        project = project,
        min.cells = min.cell,
        min.features = min.feature,
        meta.data = meta_data
        )
      
      # Add cell id as sample name
      sobj[['Sample']] <- sample
      sobj <- Seurat::RenameCells(sobj, add.cell.id = sample)
      
      # Store to list
      sobj_list[[sample]] <- sobj
      
      message("Seurat object is created for ", sample, '.')
      
      # Save if save is enabled
      if (save){qs::qsave(sobj, file = file.path(save_dir, paste0(sample, "_rna.qs")))
        message("Saved object: ", sample, " to ", file.path(save_dir, paste0(sample, "_rna.qs")), ".")
      }
    }
    message("Created ", length(rna_list), " Seurat objects.")
    return(sobj_list)
  }
