# read this function with 
# source('path/to/this/file/08_merge_atac_id.R')

# function from Matt Cannon ---rrrsingglecellepxp library
# peak_beds List of paths to peak bed files
# min_peak_width Minimum width of peaks to include
# max_peak_width Maximum width of peaks to include  
# frag_paths List of paths to fragment files
# cell_ids Sample IDs for labeling
# n_regions_simul Number of regions to process simultaneously
# threads Number of parallel threads to use
# return Merged Seurat object with ATAC data
merge_atac_id <- function(peak_beds,
                          min_peak_width = 20,
                          max_peak_width = 10000,
                          frag_paths,
                          cell_ids,
                          n_regions_simul = 2000,
                          threads = 1) {
  message("Make sure that paths and cell_ids are in the same order")
  message("peak_beds: ", paste(peak_beds, sep = " "))
  message("frag_paths: ", paste(frag_paths, sep = " "))
  message("cell_ids: ", paste(cell_ids, sep = " "))
  
  # read peaks into a dataframe and reduce them
  reduced_peaks <- readr::read_tsv(peak_beds,
                                   comment = "#",
                                   col_names = c("chr", "start", "end"),
                                   show_col_types = FALSE) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%
    Signac::reduce()
  
  # filter out peaks that are too small or too large
  reduced_peaks <- reduced_peaks[
    reduced_peaks@ranges@width >= min_peak_width &
      reduced_peaks@ranges@width <= max_peak_width
  ]
  
  # Function to get fragment files, count and make Seurat object
  make_chrom_seurat <- function(i) {
    frags <- Signac::CreateFragmentObject(path = frag_paths[[i]])
    
    counts <- Signac::FeatureMatrix(fragments = frags,
                                    features = reduced_peaks,
                                    process_n = n_regions_simul)
    
    obj <- Signac::CreateChromatinAssay(counts = counts,
                                        fragments = frags) %>%
      SeuratObject::CreateSeuratObject(assay = "ATAC")
    
    return(obj)
  }
  
  # Use apply to make the Seurat objects and then merge them
  obj_list <-
    parallel::mclapply(seq_len(length(frag_paths)),
                       make_chrom_seurat,
                       mc.cores = threads)
  
  seurat_obj <- merge(x = obj_list[[1]],
                      y = obj_list[-1],
                      add.cell.ids = cell_ids)
  
  return(seurat_obj)
}