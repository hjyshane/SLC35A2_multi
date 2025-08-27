# read this function with 
# source('path/to/this/file/03_rna_qc.R'

#' Filter Seurat objects by RNA count and mitochondrial content
#' PercentageFeatureSet calculates the fraction of counts attributed to mitochondrial genes for each cell
#' mt.pattern - Pattern to identify mitochondrial genes (default: "^mt-"), create mitocondria gene information.
#'
#' Applies basic quality control filters to a list of Seurat objects, including
#' thresholds on total RNA counts and percent mitochondrial RNA. Optionally saves
#' the filtered objects using `qs::qsave()`.
#' 
#' @param sobj_list A named list of Seurat objects to be filtered.
#' @param samples Character vector of sample names corresponding to the objects.
#' @param assay Assay to apply mitochondrial calculation to (default: "RNA").
#' @param rna_count_high Numeric. Max nCount_RNA threshold (default: 25000).
#' @param rna_count_low Numeric. Min nCount_RNA threshold (default: 1000).
#' @param mt_rna Numeric. Max mitochondrial percent threshold (default: 5).
#' @param save Logical. If TRUE, saves filtered objects to `save_dir` (default: FALSE).
#' @param save_dir_qs Directory where `.qs` files should be saved (if `save = TRUE`).
#' @param save_dir_qc Directory where qc files should be saved (if `save = TRUE`).
#'
#' @return A list of filtered Seurat objects.
#' 
#'@examples
#' filtered_list <- filter_seurat_by_qc(
#'   sobj_list = seurat_list,
#'   samples = c("PTZ_1hr", "SAL_24hr"),
#'   save = TRUE,
#'   save_dir = "~/project/filtered"
#' )
#'
#' @export
rna_qc <- function(
    sobj_list,
    samples,
    nCount_RNA_high = 25000,
    nCount_RNA_low = 1000,
    nFeature_RNA_low = 200,
    nFeature_RNA_high = 7000,
    mt_rna = 5,
    save = TRUE,
    save_dir_qs = NULL,
    save_dir_qc = NULL
    ) {
  # Save check
  if (save) {if (is.null(save_dir_qs)) {stop("You must provide 'save_dir_qs' when save = TRUE.")}
    if (is.null(save_dir_qc)) {stop("You must provide 'save_dir_qc' when save = TRUE.")}
    if (!dir.exists(save_dir_qs)) {dir.create(save_dir_qs, recursive = TRUE)}
    if (!dir.exists(save_dir_qc)) {dir.create(save_dir_qc, recursive = TRUE)}
    }
  
  # Create QC summary 
  qc_summary <- data.frame(
    sample = character(),
    before = integer(),
    after = integer(),
    removed = integer()
  )
  
  # Make sure length of samples and sobj_list is the same
  if (length(samples) != length(sobj_list)){
    stop("Sample length and objects list length are not the same. \nSample length: ", length(samples), "\nObject list length: ", length(sobj_list))
    } else {
      for (i in seq_along(sobj_list)) {
        sample <- samples[i]
        sobj <- sobj_list[[sample]]
        
        # Calculate mt rna ratio
        message("Filtering Seurat object for: ", sample)
        sobj <- Seurat::PercentageFeatureSet(
          sobj,
          pattern = '^mt-',
          col.name = 'percent_mt_rna',
          assay = "RNA") 
        
        sobj <-  sobj %>% subset(
          subset = 
            nCount_RNA < nCount_RNA_high &
            nCount_RNA > nCount_RNA_low &
            nFeature_RNA > nFeature_RNA_low & 
            nFeature_RNA < nFeature_RNA_high & 
            percent_mt_rna < mt_rna)
        
        # Basic report 
        before <- ncol(sobj_list[[sample]])
        after <- ncol(sobj)
        message("Filtered ", sample, ": ", before - after, " cells removed (", after, " kept)")
        
        # Add to table
        qc_summary <- rbind(
          qc_summary,
          data.frame(
            sample = sample,
            before = before,
            after = after,
            removed = before - after))
        
        # Store
        sobj_list[[sample]] <- sobj
        
        if (save) {qs::qsave(sobj, file = file.path(save_dir_qs, paste0(sample, "_rna_filtered.qs")))
          message("Saved object: ", sample, " to ", file.path(save_dir_qs, paste0(sample, "_rna_filtered.qs")), ".")
        }
      }
      
      if (save) {readr::write_csv(qc_summary, file = paste0(save_dir_qc, "QC_summary.csv"))
        message("QC summary saved to: ", file.path(save_dir_qc, paste0(sample, "_rna_filtered.qs")))
        }
      
      return(sobj_list)
    }
}
