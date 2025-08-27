# read this function with 
# source('path/to/this/file/08_atac_qc.R')

#' Run ATAC QC filtering based on TSS, nucleosome signal, and total counts
#'
#' @param atac_obj Seurat object with ATAC assay.
#' @param TSS_cutoff Minimum TSS enrichment (default = 2).
#' @param NSS_cutoff Maximum nucleosome signal (default = 2).
#' @param nCount_ATAC_low Minimum fragment count (default = 1000).
#' @param nCount_ATAC_high Maximum fragment count (default = 100000).
#' @param save Logical. Whether to save the filtered object (default = FALSE).
#' @param save_dir Directory to save the filtered object if `save = TRUE`.
#'
#' @return Filtered Seurat object with added TSS/NSS metadata.
#' @export
atac_qc <- function(
    atac_obj,
    TSS_cutoff = 2,
    NSS_cutoff = 4,
    nCount_ATAC_low = 1000,
    nCount_ATAC_high = 100000,
    nFeature_ATAC_low = NULL,
    nFeature_ATAC_high = NULL,
    save = TRUE,
    save_dir = NULL
    ) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
    }
  
  DefaultAssay(atac_obj) <- "ATAC"
  
  # Calculate for QC
  filtered_atac <- atac_obj %>%
    # nucleosome_signal: Indicator of chromatin accessibility.
    Signac::NucleosomeSignal() %>%
    # TSS.enrichment: Indicator of transcriptional activity.
    Signac::TSSEnrichment(fast = FALSE)
  
  # Add to meta.data
  filtered_atac$TSS <- ifelse(filtered_atac$TSS.enrichment > TSS_cutoff, "High", "Low")
  filtered_atac$NSS <- ifelse(filtered_atac$nucleosome_signal > NSS_cutoff, "High", "Low")
  
  # Subset
  filtered_atac <- subset(
    filtered_atac,
    subset = 
      nCount_ATAC < nCount_ATAC_high &
      nCount_ATAC > nCount_ATAC_low &
      # nFeature_ATAC < nFeature_ATAC_high &
      # nFeature_ATAC > nFeature_ATAC_low &
      TSS.enrichment > TSS_cutoff &
      nucleosome_signal < NSS_cutoff)
  
  # Save object
  if (save) {qs::qsave(filtered_atac, file = file.path(save_dir, "filtered_atac_obj.qs"))
    message('Saved filtered ATAC object to ', file.path(save_dir, "filtered_atac_obj.qs"))
  }
  
  return(filtered_atac)
}
