# read this function with 
# source('path/to/this/file/08_create_atac.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0428_run/Scripts/CodeSource/08_merge_atac_id.R')

#' Create ATAC object from 10X multiome data
#'
#' Merges ATAC peak and fragment files across samples, adds reference genome annotations.
#'
#' @param input_dir Base path containing sample folders.
#' @param samples Character vector of sample folder names.
#' @param ref_genome_path Path to reference GTF file.
#' @param ref_genome_name Genome build name (e.g., "mm10" or "hg38").
#' @param save Logical. If TRUE, saves ATAC object and annotation.
#' @param save_dir Directory to save .qs outputs if save = TRUE.
#'
#' @return Merged and annotated ATAC Seurat object.
#' @export
#' 
create_atac <- function(
    input_dir,
    samples,
    ref_genome_path,
    ref_genome_name,
    save = TRUE,
    save_dir = NULL) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  
  # Set up ATAC file paths
  peak_beds <- character(length(samples))
  frag_paths <- character(length(samples))
  
  for (i in seq_along(samples)){
    sample <- samples[i]
    peak_beds[i] <- file.path(input_dir, sample, "outs/atac_peaks.bed")
    frag_paths[i] <- file.path(input_dir, sample, 'outs/atac_fragments.tsv.gz')
  }
  
  # Using rrrSingleCellUtils library
  # peak_beds List of paths to peak bed files
  # min_peak_width Minimum width of peaks to include
  # max_peak_width Maximum width of peaks to include  
  # frag_paths List of paths to fragment files
  # cell_ids Sample IDs for labeling
  # n_regions_simul Number of regions to process simultaneously
  # threads Number of parallel threads to use
  # return Merged Seurat object with ATAC data
  atac_obj <- merge_atac_id(
    peak_beds = peak_beds, 
     frag_paths = frag_paths, 
     cell_ids = samples, 
     n_regions_simul = 100000, 
     threads = 3)
  
  # Assemble ref genome for annotation
  ref_genome <- rtracklayer::import(ref_genome_path)
  GenomeInfoDb::genome(ref_genome) <- ref_genome_name
  ref_genome$gene_biotype <- ref_genome$gene_type
  annotation <- ref_genome
  
  # Add annotation
  Seurat::DefaultAssay(atac_obj) <- "ATAC"
  Signac::Annotation(atac_obj) <- annotation

  # Save
  if (save) {
    qs::qsave(atac_obj, file = file.path(save_dir, "atac_obj.qs"))
    message("Saved merged atac object to ", file.path(save_dir, "atac_obj.qs"), ".")
    qs::qsave(annotation, file = file.path(save_dir, "atac_ref_annotation.qs"))
    message("Saved reference annotation for atac object to ", file.path(save_dir, "atac_obj.qs"), ".")
  }
  return(atac_obj)
  }

