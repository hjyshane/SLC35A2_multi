# read this function with 
# source('path/to/this/file/16_filter_cluster.R')


#' Filter Seurat object by minimum cell count per group
#'
#' This function filters clusters (e.g., cell types) in a Seurat object that have fewer than 
#' `number_to_filter` cells per sample group. It saves and returns a subsetted Seurat object.
#'
#' @param seurat_obj A Seurat object to filter.
#' @param sample Character. Name of the metadata column indicating sample labels (e.g., `"Sample"`).
#' @param cell_type Character. Name of the metadata column indicating cluster or cell type (e.g., `"Cell_type_by_multimodal"`).
#' @param number_to_filter Integer. Minimum number of cells per group to retain a cluster. Default is 50.
#' @param save Logical. Whether to save the filtered object to disk. Default is TRUE.
#' @param save_dir Character. Directory to save the filtered object (`.qs` file). Required if `save = TRUE`.
#'
#' @return A filtered Seurat object containing only the clusters that meet the filtering criteria.
#'
#' @examples
#' filtered_50 <- filter_cluster(seurat_obj, sample = "Sample", cell_type = "Cell_type_by_multimodal", number_to_filter = 50, save_dir = "output/")
#'
#' @export
#' 
filter_cluster <- function(
    seurat_obj,
    sample,
    cell_type,
    number_to_filter = 50,
    save = TRUE,
    qsave_dir = NULL,
    csv_dir = NULL) {
  # Save check
  if (save) {if (is.null(qsave_dir)) {stop("You must provide 'qsave_dir' when save = TRUE.")}
    if (!dir.exists(qsave_dir)) {dir.create(qsave_dir, recursive = TRUE)}
  }
  if (save) {if (is.null(csv_dir)) {stop("You must provide 'csv_dir' when save = TRUE.")}
    if (!dir.exists(csv_dir)) {dir.create(csv_dir, recursive = TRUE)}
  }
  # Cell count data frame
  cell_counts <- table(seurat_obj[[cell_type]][,1], seurat_obj[[sample]][,1])
  cell_counts_df <- as.data.frame.matrix(cell_counts)

  # Add total per row
  cell_counts_df$Total <- rowSums(cell_counts_df)
  
  # Add cell type as column
  cell_counts_df <- tibble::rownames_to_column(cell_counts_df, var = "Cell_Type")
  
  # Reorder columns: all sample columns + Total
  cell_counts_df <- cell_counts_df[, c("Cell_Type", setdiff(colnames(cell_counts_df), c("Cell_Type", "Total")), "Total")]
  
  # cluster
  clusters_to_keep <-  cell_counts_df$Cell_Type[apply(cell_counts_df[,-1], 1, function(x) min(x) >= number_to_filter)]
  
  # Extract metadata column as vector
  meta_vec <- FetchData(seurat_obj, vars = cell_type)[[1]]
  
  # Identify matching cell barcodes
  cells_to_keep <- colnames(seurat_obj)[meta_vec %in% clusters_to_keep]

  # Filter the Seurat object
  filtered <- subset(seurat_obj, cells = cells_to_keep)

  # Save
  if (save) {
    qs::qsave(filtered, file = file.path(qsave_dir, paste0("filtered_", number_to_filter, "_cells.qs")))
    message("Saved object to ", file.path(qsave_dir, paste0("filtered_", number_to_filter, "_cells.qs")), ".")
    write.csv(cell_counts_df, file = file.path(csv_dir, 'Cell_counts_table.csv'))
    message("Saved cell count table to ", file.path(csv_dir, 'Cell_counts_table.csv'), ".")
    
  }
  return(filtered)
  }
