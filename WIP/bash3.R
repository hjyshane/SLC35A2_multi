# Parallelization for process
library(future)
options(future.globals.maxSize = 128000 * 1024 ^ 2)  # Increase memory limit
library(parallel)  # Base R parallelization package

# Single-cell analysis (scRNA-seq & scATAC-seq)
library(Seurat)         # scRNA-seq analysis
library(sctransform)    # Normalization for scRNA-seq
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate

library(qs)

set.seed(42)

script_50 <- readRDS("~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/RDS_mid/combined_merged_processed_multimodal_processed_linked_annotated_50_cells.rds")
DefaultAssay(script_50) <- "RNA"

script_50 <- Seurat::SCTransform(
    script_50,
    method = 'glmGamPoi',
    assay = 'RNA',
    vars.to.regress = c('percent_mt_rna'),
    vst.flavor = 'v2',
    return.only.var.genes = F,
    new.assay.name = "SCT",
    verbose = T)
saveRDS(script_50, file = "~/PTZ_ATAC_scRNA_072024/WIP/FullScript_test/RDS_mid/combined_merged_processed_multimodal_processed_linked_annotated_50_cells_sctall.rds")

module_50 <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/ModuleTest/qsave/filtered_50.qs")
DefaultAssay(module_50) <- "RNA"

module_50 <- Seurat::SCTransform(
    module_50,
    method = 'glmGamPoi',
    assay = 'RNA',
    vars.to.regress = c('percent_mt_rna'),
    vst.flavor = 'v2',
    return.only.var.genes = F,
    new.assay.name = "SCT",
    verbose = T)

qs::qsave(module_50, file = "~/PTZ_ATAC_scRNA_072024/WIP/ModuleTest/qsave/filtered_50_sctall.qs")

