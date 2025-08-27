# Cell typing Annotation
library(tidyverse)
library(Seurat)
library(qs)

# Load custom function scripts
source('~/PTZ_ATAC_scRNA_072024/WIP/0428_run/Scripts/CodeSource/13_ref_mousebrain_annotation.R')

# Directory setup
save_dir <- '~/PTZ_ATAC_scRNA_072024/WIP/Modularized'
qsave_dir <- file.path(save_dir, "qsave")

# Read object
ref_obj <-  readRDS("/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/Jesse_mousebrain_RDS//mouseatlas_cut_processed.rds")
combined_object <- qs::qread(file.path(qsave_dir, "linked_obj.qs"))

# From mousebrin.org
# Add reference cell type data
combined_object <- ref_mousebrain_annotation(
  ref_obj = ref_obj,
  sobj = combined_object,
  save = TRUE,
  save_dir = qsave_dir)

# Ref cell type prediction will be stored in seurat_obj$predicted.id.
Idents(combined_object) <- combined_object$predicted.id

# Visualize
DimPlot(combined_object, reduction = "Multi_UMAP", label = T, repel = T) + NoLegend()


# Save the list of cell types into a variable by getting unique cell types from the object
long_names <- unique(combined_object$predicted.id)

# Create a short names list- inspect long_names and shorten it manually.
short_names <- c(
  "Ex_Neu_CA1", 
  "Gran_NB_DG", 
  "Gran_Neu_DG",
  "HP_CTX", 
  "HP_CTX",
  "Ex_Neu_Cortex",
  "CGE",
  "DG_Glial",
  "In_HP_CTX",
  "OPC",
  "Astrocytes",
  "Ex_Neu_CA3", 
  "COP", 
  "HP", 
  "Neu_Progenitor",
  "Microglia",
  "MGE",
  "Mature_Oligo",
  "MFOL",
  "In_HP_CTX", 
  "Astrocytes",
  "In_HP",
  "CGE",
  "Int_HP_CTX",
  "Microglia_Act",
  "Vascular cells",
  "Peri_Macro",
  "NFOL",
  "Vascular cells", 
  "Ependymal", 
  "HP",
  "Macrophages", 
  "Pericytes",
  "Vascular cells")

# Name vector for mapping
name_map <- setNames(short_names, long_names)

# Function to replace long names to short names
shortened <- function(name) {
  ifelse(name %in% names(name_map), name_map[name], name)
}

# Apply function
combined_object$mouseref.short <- shortened(combined_object$predicted.id)

# Check UMAP
Idents(combined_object) <- "mouseref.short"
DimPlot(combined_object, reduction = "Multi_UMAP", label = T, repel = T) + NoLegend()

Idents(combined_object) <- "predicted.id"
DimPlot(combined_object, reduction = "Multi_UMAP", label = T, repel = T) + NoLegend()
