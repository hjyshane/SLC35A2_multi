# read this function with 
# source('path/to/this/file/12_link_atac_rna.R')

#' Link ATAC peaks to gene activity using Signac
#'
#' @param processed_comb Seurat object with ATAC and RNA assays
#' @param genome_data BSgenome object (e.g., BSgenome.Mmusculus.UCSC.mm10)
#' @param save Logical. Whether to save the annotated object.
#' @param save_dir Directory to save the annotated object if `save = TRUE`.
#'
#' @return Seurat object with linked gene activity (ACTIVITY assay)
#' @export
link_peak_rna <- function(
    processed_comb,
    genome_data = BSgenome.Mmusculus.UCSC.mm10,
    save = TRUE,
    save_dir = NULL
    ) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  
  # ref genome library check
  if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)) {
    stop("The package 'BSgenome.Mmusculus.UCSC.mm10' is not installed.\nYou can install it with BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')")
  } else {
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
  }
  
  Seurat::DefaultAssay(processed_comb) <- "ATAC"
  
  # Gene activity matrix calculation and add to assay object
  gene_activity <- Signac::GeneActivity(processed_comb)
  processed_comb[['ACTIVITY']] <- Seurat::CreateAssayObject(counts = gene_activity)
  
  # Normalize
  processed_comb <- Seurat::NormalizeData(processed_comb, assay = 'ACTIVITY',
                                          verbose = TRUE)

  # Get sequence for peaks
  processed_comb <- Signac::RegionStats(object = processed_comb,
                                        assay = 'ATAC',
                                        genome = genome_data)
  
  # change ATAC assay annotation data to have changed column id
  # get assay data
  annotation_data <- Signac::Annotation(processed_comb[['ATAC']])
  
  # rename
  mcols(annotation_data)$tx_id <- mcols(annotation_data)$transcript_id
  mcols(annotation_data)$transcript_id <- NULL
  
  # Add back to ATAC assay
  Signac::Annotation(processed_comb) <- annotation_data
  
  # rename object
  linked_obj <- processed_comb
  
  # Save object
  if (save) {qs::qsave(processed_comb, file = file.path(save_dir, "linked_obj.qs"))
    message('Saved linked object to: ', file.path(save_dir, "linked_obj.qs"), ".")
  }
  return(linked_obj)
}