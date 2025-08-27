# Library load
library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(rrrSingleCellUtils)
dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")  # Load HDF5 library
library(hdf5r)
library(qs)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnhancedVolcano)
library(clusterProfiler)
library(rrvgo)
library(enrichplot)


# Load custom function scripts
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/01_load_10x_data.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/02_create_seurat_object.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/03_rna_qc.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/04_process_rna.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/05_rna_integration.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/06_process_integrated.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/07_create_atac.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/08_merge_atac_id.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/08_atac_qc.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/09_process_atac.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/10_atac_rna_merge.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/11_process_multi.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/12_link_atac_rna.R')
source('~/PTZ_ATAC_scRNA_072024/WIP/0609_run/CodeSource/16_filter_cluster.R')


# set seed
set.seed(42)

# Directory setup
input_dir <- '~/PTZ_ATAC_scRNA_072024/File'
save_dir <- '~/PTZ_ATAC_scRNA_072024/WIP/0624_run'
qsave_dir <- file.path(save_dir, "qsave")
qc_dir    <- file.path(save_dir, "QC")
csv_dir <- file.path(save_dir, 'CSV')

dir.create(save_dir, recursive = TRUE)
dir.create(qsave_dir, recursive = TRUE)
dir.create(qc_dir, recursive = TRUE)

# Set sample list
samples <- c("PTZ_1hr", "PTZ_24hr", "SAL_1hr", "SAL_24hr")


# Load H5
rna_list <- load_10x_data(
    input_dir,
    assay = 'Gene Expression',
    samples,
    save = T,
    save_dir = qsave_dir)

# Create seurat object
seurat_list <- create_seurat_object(
    rna_list,
    input_dir,
    samples,
    project = 'SeuratProject',
    min.cell = 3,
    min.feature = 200,
    save = T,
    save_dir = qsave_dir)

# QC
filtered_list <- rna_qc(
    seurat_list,
    samples,
    nCount_RNA_high = 25000,
    nCount_RNA_low = 1000,
    mt_rna = 5,
    save = T,
    save_dir_qs = qsave_dir,
    save_dir_qc = qc_dir)

# Process RNA
processed_list <- process_rna(
    filtered_list,
    samples,
    save = T,
    save_dir = qsave_dir)

# If you want to save as list
qs::qsave(rna_list, file = file.path(qsave_dir, 'H5_rna_list.qs'))
qs::qsave(seurat_list, file = file.path(qsave_dir, 'seurat_rna_list.qs'))
qs::qsave(filtered_list, file = file.path(qsave_dir, 'filtered_rna_list.qs'))
qs::qsave(processed_list, file = file.path(qsave_dir, 'processed_rna_list.qs'))


# Integrate
integrated_rna <- rna_integration(
    processed_list,
    save = T,
    save_dir = qsave_dir)


# Process integrated object
processed_rna <- process_integrated(
    integrated_rna,
    range_start = 0.4,
    range_end = 1.2,
    range_step = 0.1,
    save = T,
    qsave_dir = qsave_dir,
    markers_dir = csv_dir)

# Create ATAC object
atac_obj <- create_atac(
    input_dir,
    samples,
    ref_genome_path = "/igm/apps/10X_chromium/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz",
    ref_genome_name = "mm10",
    save = TRUE,
    save_dir = qsave_dir)

# ATAC QC
filtered_atac <- atac_qc(
    atac_obj,
    TSS_cutoff = 2,
    NSS_cutoff = 2,
    nCount_ATAC_low = 1000,
    nCount_ATAC_high = 100000,
    save = TRUE,
    save_dir = qsave_dir)

# Process ATAC
processed_atac <- process_atac(
    filtered_atac,
    range_start = 0.4,
    range_end = 1.2,
    range_step = 0.1,
    save = TRUE,
    qsave_dir = qsave_dir,
    markers_dir = csv_dir)

# Combine ATAC and RNA object
combined_object <- atac_rna_merge(
    processed_rna,
    processed_atac,
    save = TRUE,
    save_dir = qsave_dir)

# Process combined object and do MultiModalNeighbor process
processed_comb <- process_multi(
    combined_object,
    range_start = 0.1,
    range_end = 1.2,
    range_step = 0.1,
    save = TRUE,
    qsave_dir = qsave_dir,
    markers_dir = csv_dir)

# Link ATAC peaks and RNA gene expression
linked_obj <- link_peak_rna(
    processed_comb,
    genome_data = BSgenome.Mmusculus.UCSC.mm10,
    save = TRUE,
    save_dir = qsave_dir
)

combined_object <- qs::qread('~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/linked_obj.qs')

Seurat::Idents(combined_object) <- combined_object$Multi_cluster_0.4
DimPlot(combined_object, reduction = "Multi_UMAP", label = T) + NoLegend()

multimodal_cell_type_mapping <- c(
    "1" = "ExcitatoryNeuronsMatureDG", "12" = "ExcitatoryNeuronsMatureDG",

    "2" = "ExcitatoryNeuronsImmatureDG",

    "7" = "OPCs",

    "20" = "COPMFOL",

    "10" = "CGE",

    "22" = "Lamp5positive",

    "4" = "MGE",

    "17" = "LGE",

    "5" = "Astrocytes",

    "21" = "Microglia",

    "9" = "L23ExcitatoryNeurons", "13" = "L23ExcitatoryNeurons",

    "11" = "L45ExcitatoryNeurons",  "15" = "L45ExcitatoryNeurons", "16" = "L45ExcitatoryNeurons",

    "6" = "L56ExcitatoryNeurons", "23" = "L56ExcitatoryNeurons",

    "3" = "ExcitatoryNeuronsCA3", "18" = "ExcitatoryNeuronsCA3",

    "0" = "ExcitatoryNeuronsCA1", "8" = "ExcitatoryNeuronsCA1", "14" = "ExcitatoryNeuronsCA1",

    "19" = "NdnfRelnInhibitoryInterneurons")

# Assign the cell types to the metadata based on cluster IDs
combined_object@meta.data$cell_type <- multimodal_cell_type_mapping[as.character(combined_object$Multi_cluster_0.4)]

Idents(combined_object) <- "cell_type"  #"Multi_cluster_0.4" "Cell_type_by_multimodal
multi_anno <- DimPlot(combined_object, reduction = "Multi_UMAP", label = TRUE, repel = TRUE) + NoLegend()
combined_object[['orig.ident']] <- combined_object[['Sample']]

combined_object@meta.data <-
    combined_object@meta.data %>%
    mutate(
        condition = str_split_fixed(Sample, "_", 2)[, 1],
        timepoint = str_split_fixed(Sample, "_", 2)[, 2]
    )

qs::qsave(combined_object, file = file.path(qsave_dir, 'combined_cell_type.qs'))

