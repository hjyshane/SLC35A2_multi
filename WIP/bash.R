# Parallelization for process
library(future)
options(future.globals.maxSize = 128000 * 1024 ^ 2)  # Increase memory limit
library(parallel)  # Base R parallelization package

# Single-cell analysis (scRNA-seq & scATAC-seq)
library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(sctransform)    # Normalization for scRNA-seq
library(DoubletFinder)  # Detects doublets in single-cell datasets
library(scCustomize)    # Enhances Seurat visualization and analysis
library(SingleCellExperiment)

# Data manipulation & utilities
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(data.table)  # High-speed alternative to data frames
library(Matrix)      # Sparse matrix operations (important for single-cell data)
library(jsonlite)

# Visualization & Plotting
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots
library(cowplot)          # Additional ggplot utilities
library(ggrepel)          # Prevents text overlap in plots
library(ggvenn)           # Venn diagrams
library(ggVennDiagram)    # Another Venn diagram package
library(VennDiagram)      # Classic Venn diagrams
library(ggpubr)           # Publication-ready plots
library(ComplexHeatmap)   # Advanced heatmaps
library(circlize)         # Circular visualizations
library(pheatmap)         # Heatmaps
library(RColorBrewer)     # Color palettes
library(fields)           # Spatial plots and interpolation
library(enrichplot)       # Visualization for enrichment analysis
library(eulerr)           # Euler diagrams
library(gridExtra)        # Arrange multiple plots
library(corrplot)         # Correlation heatmaps
library(viridis)          # Colorblind-friendly palettes
library(scales)           # scale plot
library(EnhancedVolcano)
library(RColorBrewer)

# Clustering & Dimensionality Reduction
library(clustree)  # Clustering trees for resolution selection
library(igraph)    # Graph-based clustering

# Differential Expression & Statistical Analysis
library(DESeq2)   # Differential expression analysis
library(MAST)
library(rstatix)  # Statistical functions for ggplot2
library(ROCR)     # ROC curves for model evaluation

# Gene Ontology & Pathway Analysis
library(clusterProfiler)  # Functional enrichment analysis
library(org.Mm.eg.db)     # Mouse genome annotations
library(DOSE)             # Disease ontology enrichment
library(rrvgo)            # Reduce GO terms redundancy
library(treemap)          # Tree visualization for GO terms
library(AnnotationDbi)    # Database interface
library(GO.db)            # Gene Ontology database
library(ontologyIndex)    # General ontology functions
library(fgsea)            # fGSEA analysis
library(msigdbr)          # ref for fGSEA

# Pseudotime trajectory analysis
library(slingshot)
library(tradeSeq)

# Genomics & Motif Analysis
library(BSgenome.Mmusculus.UCSC.mm10)  # Mouse genome
library(rtracklayer)                   # Import/export genomic data
library(JASPAR2020)                    # Transcription factor motifs
library(TFBSTools)                      # Analyze transcription factor binding sites

# File Handling & Misc
library(openxlsx)  # Read/write Excel files
dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")  # Load HDF5 library
library(hdf5r)     # HDF5 file format for large datasets

# Interactive Applications
library(shiny)  # Web-based interactive applications

set.seed(42)

# Output folder path set
output_dir_base <-"~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test"

# Set base input folder directory files.
input_dir <- "~/PTZ_ATAC_scRNA_072024/File"

dir.create(output_dir_base, recursive = TRUE, showWarnings = FALSE)

# Set sample ID
Samples <- c("PTZ_1hr", "PTZ_24hr", "SAL_1hr", "SAL_24hr")

integrated_rna <- readRDS('~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/RDS_mid/integrated_rna_processed.rds')

DefaultAssay(integrated_rna) <- "SCT"

marker_results <- list()
for (i in seq(0.4, 1.2, by = 0.1)) {
  Idents(integrated_rna) <- integrated_rna[[paste0('RNA_cluster_', i)]][,1]
  marker_results[[paste0("RNA_cluster_", i)]] <- FindAllMarkers(integrated_rna,
                                                                assay = "SCT",
                                                                slot = "data",
                                                                only.pos = TRUE,
                                                                min.pct = 0.25,
                                                                logfc.threshold = 0.25)
}

for (name in names(marker_results)) {
  write.csv(marker_results[[name]], file = paste0("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/csv/", name, "_markers.csv"), row.names = FALSE)
}

message('integration done')

peak_beds <- c(
  file.path(input_dir, Samples[1], "outs/atac_peaks.bed"),
  file.path(input_dir, Samples[2], "outs/atac_peaks.bed"),
  file.path(input_dir, Samples[3], "outs/atac_peaks.bed"),
  file.path(input_dir, Samples[4], "outs/atac_peaks.bed"))

frag_paths <- c(
  file.path(input_dir, Samples[1], "outs/atac_fragments.tsv.gz"),
  file.path(input_dir, Samples[2], "outs/atac_fragments.tsv.gz"),
  file.path(input_dir, Samples[3], "outs/atac_fragments.tsv.gz"),
  file.path(input_dir, Samples[4], "outs/atac_fragments.tsv.gz"))

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

merged_atac <- merge_atac_id(peak_beds = peak_beds, frag_paths = frag_paths, cell_ids = Samples, n_regions_simul = 100000, threads = 3)
saveRDS(merged_atac,file.path(output_dir_base, "RDS_mid", "merged_atac.rds"))



mm10 <-import("/igm/apps/10X_chromium/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz")

genome(mm10) <- "mm10"
seqlevelsStyle(mm10) <- "UCSC"
mm10$gene_biotype <- mm10$gene_type
annotations <- mm10


DefaultAssay(merged_atac) <- "ATAC"
Annotation(merged_atac) <- annotations

# Claer workspace
rm("i", "name","mm10", "annotations", "peak_beds", "frag_paths")

saveRDS(annotation, file.path(output_dir_base, "RDS_mid", "atac_annotation.rds"))
saveRDS(merged_atac,file.path(output_dir_base, "RDS_mid", "merged_atac_annotated.rds"))

merged_atac <- NucleosomeSignal(merged_atac)
merged_atac <- TSSEnrichment(merged_atac, fast = FALSE)
merged_atac$high_tss <- ifelse(merged_atac$TSS.enrichment > 2, "High", "Low")
merged_atac$nucleosome_group <- ifelse(merged_atac$nucleosome_signal > 4,paste0("NS > ", 4), paste0("NS < ", 4))
merged_atac <- subset(merged_atac, subset = nCount_ATAC < 100000 & nCount_ATAC > 1000 & nucleosome_signal < 2 & TSS.enrichment > 2)
saveRDS(merged_atac,file.path(output_dir_base, "RDS_mid", "merged_atac_annotated_filtered.rds"))

DefaultAssay(merged_atac) <- "ATAC"
merged_atac <- RunTFIDF(merged_atac)
merged_atac <- FindTopFeatures(merged_atac, min.cutoff = 10) # Select features
merged_atac <- RunSVD(merged_atac) # Generate 'lsi' reduction
merged_atac <- FindNeighbors(merged_atac, reduction = "lsi", dims = 2:30, graph.name = c("ATAC_nn", "ATAC_snn"), verbose = T)
merged_atac <- RunUMAP(merged_atac, reduction.name = "ATAC_umap", reduction = "lsi", dims = 1:30, verbose = T)

# Cluster
for (i in seq(0.4, 1.2, by = 0.1)) {
  merged_atac <- FindClusters(merged_atac, graph.name = "ATAC_snn", resolution = i, verbose = T, cluster.name = paste0("ATAC_cluster_", i))
}
saveRDS(merged_atac,file.path(output_dir_base, "RDS_mid", "merged_atac_annotated_filtered_processed.rds"))
DefaultAssay(merged_atac) <- "ATAC"

marker_results_atac <- list()
for (i in seq(0.4, 1.2, by = 0.1)) {
  Idents(merged_atac) <- merged_atac[[paste0('ATAC_cluster_', i)]][,1]
  marker_results_atac[[paste0("ATAC_cluster_", i)]] <- FindAllMarkers(merged_atac,
                                                                      assay = "ATAC",
                                                                      slot = "data",
                                                                      only.pos = TRUE,
                                                                      min.pct = 0.25,
                                                                      logfc.threshold = 0.25)
}

for (name in names(marker_results_atac)) {
  write.csv(marker_results_atac[[name]], file = paste0("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/csv/", name, "_markers.csv"), row.names = FALSE)
}
message('atac done')

common_cells <- intersect(colnames(integrated_rna), colnames(merged_atac))
integrated_rna <- integrated_rna[, common_cells]
merged_atac <- merged_atac[, common_cells]
integrated_rna[["ATAC"]] <- merged_atac[["ATAC"]]
atac_meta <- merged_atac@meta.data
rna_meta <- integrated_rna@meta.data
atac_specific_cols <-setdiff(colnames(atac_meta), colnames(rna_meta))
for (col in atac_specific_cols) {
  integrated_rna[[col]] <- atac_meta[[col]]}
combined_merged <- integrated_rna
rm("integrated_rna", "merged_atac", "atac_meta", "rna_meta", "common_cells", "atac_specific_cols")

saveRDS(combined_merged, file.path(output_dir_base, "RDS_mid", "combined_merged.rds"))


summary_stats <- combined_merged@meta.data %>%
  dplyr::group_by(orig.ident, condition, timepoint) %>%
  dplyr::summarize(Cell_Count = dplyr::n(), ) %>% # Count the number of rows (cells) in each group
  dplyr::ungroup() %>%
  dplyr::rename(Group = orig.ident, Condition = condition, Timepoint = timepoint)

# save
write.csv(summary_stats, file.path(output_dir_base, "csv", "RNA_ATAC_merge_summary.csv"), row.names = FALSE)

DefaultAssay(combined_merged) <- "ATAC"
combined_merged <- RunTFIDF(combined_merged)  # Normalize raw ATAC data
combined_merged <- FindTopFeatures(combined_merged, min.cutoff = 10) # Select variable peaks
combined_merged <- RunSVD(combined_merged) # Dimensionality reduction (creates "lsi")
saveRDS(combined_merged, file.path(output_dir_base, "RDS_mid", "combined_merged_processed.rds"))

DefaultAssay(combined_merged) <- "RNA"
combined_merged <- SCTransform(combined_merged, method = "glmGamPoi", assay = "RNA",  vars.to.regress = c("percent_mt_rna"), vst.flavor = "v2", verbose = TRUE)
combined_merged <- RunPCA(combined_merged, verbose = TRUE) # Dimensionality reduction (creates "pca")
combined_merged <- FindNeighbors(combined_merged, dims = 1:30, assay = "SCT", reduction = "pca",  graph.name = c("RNA_integrated_nn", "RNA_integrated_snn"), verbose = T)
combined_merged <- RunUMAP(combined_merged, reduction.name = "RNA_integrated_umap",  assay = "SCT", reduction = "pca", dims = 1:30, verbose = T)
combined_merged <- FindNeighbors(combined_merged, dims = 2:30,  assay = "ATAC", reduction = "lsi",  graph.name = c("ATAC_integrated_nn", "ATAC_integrated_snn"), verbose = T)
combined_merged <- RunUMAP(combined_merged, reduction.name = "ATAC_integrated_umap",  assay = "ATAC", reduction = "lsi", dims = 2:30, verbose = T)

for (i in seq(0.4, 1.2, by = 0.1)) {
  combined_merged <- FindClusters(combined_merged, graph.name = "RNA_integrated_snn", resolution = i, verbose = T, cluster.name = paste0("RNA_integrated_cluster_", i))
  combined_merged <- FindClusters(combined_merged, graph.name = "ATAC_integrated_snn", resolution = i, verbose = T, cluster.name = paste0("ATAC_integrated_cluster_", i))
}

saveRDS(combined_merged, file.path(output_dir_base, "RDS_mid", "combined_merged_processed.rds"))

combined_merged <- PrepSCTFindMarkers(combined_merged)

marker_results_combined <- list()
for (i in seq(0.4, 1.2, by = 0.1)) {
  Idents(combined_merged) <- combined_merged[[paste0('ATAC_integrated_cluster_', i)]][,1]
  marker_results_combined[[paste0("RNA_integrated_cluster_", i)]] <- FindAllMarkers(combined_merged,
                                                                                    assay = "ATAC",
                                                                                    slot = "data",
                                                                                    only.pos = TRUE,
                                                                                    min.pct = 0.25,
                                                                                    logfc.threshold = 0.25,
                                                                                    graph.name = "RNA_integrated_snn")

  Idents(combined_merged) <- combined_merged[[paste0('RNA_integrated_cluster_', i)]][,1]
  marker_results_combined[[paste0("ATAC_integrated_cluster_", i)]] <- FindAllMarkers(combined_merged,
                                                                                     assay = "SCT",
                                                                                     slot = "data",
                                                                                     only.pos = TRUE,
                                                                                     min.pct = 0.25,
                                                                                     logfc.threshold = 0.25,
                                                                                     graph.name = "ATAC_integrated_snn")
}
for (name in names(marker_results_combined)) {
    write.csv(marker_results_combined[[name]], file = paste0("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/csv/", name, "_markers.csv"), row.names = FALSE)
}

combined_merged <- FindMultiModalNeighbors(object = combined_merged,
                                         reduction.list = list("pca", "lsi"), # reduction.list: Specifies RNA and ATAC reductions ("pca" for RNA, "lsi" for ATAC).
                                         dims.list = list(1:30, 2:30), # dims.list: Defines the number of dimensions to use for each modality.
                                         knn.graph.name = "Multimodal_wknn", # Multimodal_wknn: Weighted k-nearest neighbor graph.
                                         snn.graph.name = "Multimodal_wsnn", # Multimodal_wsnn: Shared nearest neighbor (SNN) graph.
                                         weighted.nn.name = "Multimodal_weighted.nn", # Multimodal_weighted.nn: Final multi-modal neighbor graph.
                                         modality.weight.name = c("Multimodal_RNA_weight", "Multimodal_ATAC_weight"), # modality.weight.name: Adds weights to indicate the relative contribution of RNA and ATAC.
                                         verbose = TRUE)
saveRDS(combined_merged, file.path(output_dir_base, "RDS_mid", "combined_merged_processed_multimodal.rds"))


combined_merged <- RunUMAP(object = combined_merged,
                           nn.name = "Multimodal_weighted.nn", # nn.name: Uses the joint neighbor graph.
                           # UMAP embedding for multi-modal neighbors does not rely on the data in the SCT assay
                           # but instead leverages the Multimodal_weighted.nn graph.
                           # Hence, the choice of assay (SCT or ATAC) does not impact the multi-modal UMAP embedding.
                           verbose = TRUE,
                           reduction.name = "Multi_umap")


for (i in seq(0.4, 1.2, by = 0.1)) {
  combined_merged <- FindClusters(combined_merged, graph.name = "Multimodal_wsnn", resolution = i, verbose = T, cluster.name = paste0("Multimodal_cluster_", i))
}

saveRDS(combined_merged, file.path(output_dir_base, "RDS_mid", "combined_merged_processed_multimodal_processed.rds"))

marker_results_multi <- list()
for (i in seq(0.4, 1.2, by = 0.1)) {
  Idents(combined_merged) <- combined_merged[[paste0('Multimodal_cluster_', i)]][,1]
  marker_results_multi[[paste0("Multimodal_cluster_", i)]] <- FindAllMarkers(combined_merged,
                                                                             assay = "SCT",
                                                                             slot = "data",
                                                                             only.pos = TRUE,
                                                                             min.pct = 0.25,
                                                                             logfc.threshold = 0.25)
}
for (name in names(marker_results_multi)) {
  write.csv(marker_results_multi[[name]], file = paste0("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/csv/", name, "_markers.csv"), row.names = FALSE)
}

message('multi done')

DefaultAssay(combined_merged) <- "ATAC"
gene.activities <- GeneActivity(combined_merged)
combined_merged[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
combined_merged <- NormalizeData(combined_merged, assay = "ACTIVITY")
DefaultAssay(combined_merged) <- "ATAC"
genome <- BSgenome.Mmusculus.UCSC.mm10  # Replace with appropriate genome if not mouse
combined_merged <- RegionStats(object = combined_merged, assay = "ATAC", genome = genome)
annotation_data <- Annotation(combined_merged[["ATAC"]])
mcols(annotation_data)$tx_id <- mcols(annotation_data)$transcript_id
mcols(annotation_data)$transcript_id <- NULL  # Remove the old column if desired
Annotation(combined_merged[["ATAC"]]) <- annotation_data
saveRDS(combined_merged, file.path(output_dir_base, "RDS_mid", "combined_merged_processed_multimodal_processed_linked.rds"))


message('done')