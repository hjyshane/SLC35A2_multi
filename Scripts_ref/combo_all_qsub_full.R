##### load libraries #####

dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")
library(hdf5r)
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(clustree)
library(plyr)
library(dplyr)
library(cowplot)
library(readr)
library(Matrix)
library(sctransform)
library(glmGamPoi)
library(ROCR)
library(parallel)
library(fields)
library(DoubletFinder)
library(clustree)


ptz1hr.counts <- Read10X_h5("/igm/projects/210907_Bedrosian_cellranger-arc-count/PTZ-1hr/outs/filtered_feature_bc_matrix.h5")
sal1hr.counts <- Read10X_h5("/igm/projects/210907_Bedrosian_cellranger-arc-count/PTZ-24hr/outs/filtered_feature_bc_matrix.h5")
ptz24hr.counts <- Read10X_h5("/igm/projects/210907_Bedrosian_cellranger-arc-count/SAL-1hr/outs/filtered_feature_bc_matrix.h5")
sal24hr.counts <- Read10X_h5("/igm/projects/210907_Bedrosian_cellranger-arc-count/SAL-24hr/outs/filtered_feature_bc_matrix.h5")

ptz1hr.fragpath <- '/igm/projects/210907_Bedrosian_cellranger-arc-count/PTZ-1hr/outs/atac_fragments.tsv.gz'
sal1hr.fragpath <- '/igm/projects/210907_Bedrosian_cellranger-arc-count/PTZ-24hr/outs/atac_fragments.tsv.gz'
ptz24hr.fragpath <- '/igm/projects/210907_Bedrosian_cellranger-arc-count/SAL-1hr/outs/atac_fragments.tsv.gz'
sal24hr.fragpath <- '/igm/projects/210907_Bedrosian_cellranger-arc-count/SAL-24hr/outs/atac_fragments.tsv.gz'

annotation <- readRDS("/igm/home/hxy008/mouse_annotation.rds")

ptz1hr.metadata <- read.csv(
  file = "/igm/projects/210907_Bedrosian_cellranger-arc-count/PTZ-1hr/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

sal1hr.metadata <- read.csv(
  file = "/igm/projects/210907_Bedrosian_cellranger-arc-count/PTZ-24hr/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

ptz24hr.metadata <- read.csv(
  file = "/igm/projects/210907_Bedrosian_cellranger-arc-count/SAL-1hr/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

sal24hr.metadata <- read.csv(
  file = "/igm/projects/210907_Bedrosian_cellranger-arc-count/SAL-24hr/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

ptz1hr <- CreateSeuratObject(
  counts = ptz1hr.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = ptz1hr.metadata
)

sal1hr <- CreateSeuratObject(
  counts = sal1hr.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = sal1hr.metadata
)

ptz24hr <- CreateSeuratObject(
  counts = ptz24hr.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = ptz24hr.metadata
)

sal24hr <- CreateSeuratObject(
  counts = sal24hr.counts$`Gene Expression`,
  assay = "RNA",
  meta.data = sal24hr.metadata
)

ptz1hr[["ATAC"]] <- CreateChromatinAssay(
  counts = ptz1hr.counts$Peaks,
  sep = c(":", "-"),
  fragments = ptz1hr.fragpath,
  annotation = annotation
)

sal1hr[["ATAC"]] <- CreateChromatinAssay(
  counts = sal1hr.counts$Peaks,
  sep = c(":", "-"),
  fragments = sal1hr.fragpath,
  annotation = annotation
)

ptz24hr[["ATAC"]] <- CreateChromatinAssay(
  counts = ptz24hr.counts$Peaks,
  sep = c(":", "-"),
  fragments = ptz24hr.fragpath,
  annotation = annotation
)

sal24hr[["ATAC"]] <- CreateChromatinAssay(
  counts = sal24hr.counts$Peaks,
  sep = c(":", "-"),
  fragments = sal24hr.fragpath,
  annotation = annotation
)

DefaultAssay(ptz1hr) <- "ATAC"
DefaultAssay(sal1hr) <- "ATAC"
DefaultAssay(ptz24hr) <- "ATAC"
DefaultAssay(sal24hr) <- "ATAC"

ptz1hr <- NucleosomeSignal(ptz1hr)
sal1hr <- NucleosomeSignal(sal1hr)
ptz24hr <- NucleosomeSignal(ptz24hr)
sal24hr <- NucleosomeSignal(sal24hr)

ptz1hr <- TSSEnrichment(ptz1hr, fast = F)
sal1hr <- TSSEnrichment(sal1hr, fast = F)
ptz24hr <- TSSEnrichment(ptz24hr, fast = F)
sal24hr <- TSSEnrichment(sal24hr, fast = F)

ptz1hr$pct_reads_in_peaks <- ptz1hr$atac_peak_region_fragments / ptz1hr$atac_fragments * 100
sal1hr$pct_reads_in_peaks <- sal1hr$atac_peak_region_fragments / sal1hr$atac_fragments * 100
ptz24hr$pct_reads_in_peaks <- ptz24hr$atac_peak_region_fragments / ptz24hr$atac_fragments * 100
sal24hr$pct_reads_in_peaks <- sal24hr$atac_peak_region_fragments / sal24hr$atac_fragments * 100

ptz1hr[["percent.mt"]] <- PercentageFeatureSet(ptz1hr, pattern = "^mt-")
sal1hr[["percent.mt"]] <- PercentageFeatureSet(sal1hr, pattern = "^mt-")
ptz24hr[["percent.mt"]] <- PercentageFeatureSet(ptz24hr, pattern = "^mt-")
sal24hr[["percent.mt"]] <- PercentageFeatureSet(sal24hr, pattern = "^mt-")

ptz1hr <- subset(
  x = ptz1hr,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2 &
    percent.mt < 5
)

sal1hr <- subset(
  x = sal1hr,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2 &
    percent.mt < 5
)

ptz24hr <- subset(
  x = ptz24hr,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2 &
    percent.mt < 5
)

sal24hr <- subset(
  x = sal24hr,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2 &
    percent.mt < 5
)

DefaultAssay(ptz1hr) <- "RNA"
DefaultAssay(sal1hr) <- "RNA"
DefaultAssay(ptz24hr) <- "RNA"
DefaultAssay(sal24hr) <- "RNA"

ptz1hr <- SCTransform(ptz1hr, method = "glmGamPoi",verbose = T)
sal1hr <- SCTransform(sal1hr, method = "glmGamPoi",verbose = T)
ptz24hr <- SCTransform(ptz24hr, method = "glmGamPoi",verbose = T)
sal24hr <- SCTransform(sal24hr, method = "glmGamPoi",verbose = T)

ptz1hr <- RunPCA(ptz1hr,verbose = T)
sal1hr <- RunPCA(sal1hr,verbose = T)
ptz24hr <- RunPCA(ptz24hr,verbose = T)
sal24hr <- RunPCA(sal24hr,verbose = T)

ptz1hr <- RunUMAP(ptz1hr, dims = 1:30, verbose = T)
sal1hr <- RunUMAP(sal1hr, dims = 1:30, verbose = T)
ptz24hr <- RunUMAP(ptz24hr, dims = 1:30, verbose = T)
sal24hr <- RunUMAP(sal24hr, dims = 1:30, verbose = T)

ptz1hr <- RunTSNE(ptz1hr, dims = 1:30, verbose = T)
sal1hr <- RunTSNE(sal1hr, dims = 1:30, verbose = T)
ptz24hr <- RunTSNE(ptz24hr, dims = 1:30, verbose = T)
sal24hr <- RunTSNE(sal24hr, dims = 1:30, verbose = T)

sweep.res.ptz1hr <- paramSweep_v3(ptz1hr, PCs = 1:30, sct = T)
sweep.res.sal1hr <- paramSweep_v3(sal1hr, PCs = 1:30, sct = T)
sweep.res.ptz24hr <- paramSweep_v3(ptz24hr, PCs = 1:30, sct = T)
sweep.res.sal24hr <- paramSweep_v3(sal24hr, PCs = 1:30, sct = T)

sweep.stats.ptz1hr<- summarizeSweep(sweep.res.ptz1hr, GT = FALSE)
sweep.stats.sal1hr<- summarizeSweep(sweep.res.sal1hr, GT = FALSE)
sweep.stats.ptz24hr<- summarizeSweep(sweep.res.ptz1hr, GT = FALSE)
sweep.stats.sal24hr<- summarizeSweep(sweep.res.sal1hr, GT = FALSE)

bcmvn_ptz1hr <- find.pK(sweep.stats.ptz1hr)
bcmvn_sal1hr <- find.pK(sweep.stats.sal1hr)
bcmvn_ptz24hr <- find.pK(sweep.stats.ptz24hr)
bcmvn_sal24hr <- find.pK(sweep.stats.sal24hr)

ptz1hrpk <- as.numeric(as.vector(top_n(bcmvn_ptz1hr,1,BCmetric)[["pK"]]))
sal1hrpk <- as.numeric(as.vector(top_n(bcmvn_sal1hr,1,BCmetric)[["pK"]]))
ptz24hrpk <- as.numeric(as.vector(top_n(bcmvn_ptz24hr,1,BCmetric)[["pK"]]))
sal24hrpk <- as.numeric(as.vector(top_n(bcmvn_sal24hr,1,BCmetric)[["pK"]]))

nExp_poi.ptz1hr <- round(0.075*nrow(ptz1hr@meta.data))
nExp_poi.sal1hr <- round(0.075*nrow(sal1hr@meta.data))
nExp_poi.ptz24hr <- round(0.075*nrow(ptz24hr@meta.data))
nExp_poi.sal24hr <- round(0.075*nrow(sal24hr@meta.data))

ptz1hr <- doubletFinder_v3(ptz1hr, PCs = 1:30, pN = 0.25, pK = ptz1hrpk, nExp = nExp_poi.ptz1hr, reuse.pANN = FALSE, sct = T)
sal1hr <- doubletFinder_v3(sal1hr, PCs = 1:30, pN = 0.25, pK = sal1hrpk, nExp = nExp_poi.sal1hr, reuse.pANN = FALSE, sct = T)
ptz24hr <- doubletFinder_v3(ptz24hr, PCs = 1:30, pN = 0.25, pK = ptz24hrpk, nExp = nExp_poi.ptz24hr, reuse.pANN = FALSE, sct = T)
sal24hr <- doubletFinder_v3(sal24hr, PCs = 1:30, pN = 0.25, pK = sal24hrpk, nExp = nExp_poi.sal24hr, reuse.pANN = FALSE, sct = T)

Idents(ptz1hr) <- ptz1hr@meta.data[[paste0("DF.classifications_0.25_",ptz1hrpk,"_",nExp_poi.ptz1hr)]]
Idents(sal1hr) <- sal1hr@meta.data[[paste0("DF.classifications_0.25_",sal1hrpk,"_",nExp_poi.sal1hr)]]
Idents(ptz24hr) <- ptz24hr@meta.data[[paste0("DF.classifications_0.25_",ptz24hrpk,"_",nExp_poi.ptz24hr)]]
Idents(sal24hr) <- sal24hr@meta.data[[paste0("DF.classifications_0.25_",sal24hrpk,"_",nExp_poi.sal24hr)]]

ptz1hr <- subset(ptz1hr, idents = "Singlet")
sal1hr <- subset(sal1hr, idents = "Singlet")
ptz24hr <- subset(ptz24hr, idents = "Singlet")
sal24hr <- subset(sal24hr, idents = "Singlet")

ptz1hr$sample <- "PTZ_1hr"
sal1hr$sample <- "Saline_1hr"
ptz24hr$sample <- "PTZ_24hr"
sal24hr$sample <- "Saline_24hr"

ptz1hr$time <- "1hr"
sal1hr$time <- "1hr"
ptz24hr$time <- "24hr"
sal24hr$time <- "24hr"

ptz1hr$treat <- "PTZ"
sal1hr$treat <- "Saline"
ptz24hr$treat <- "PTZ"
sal24hr$treat <- "Saline"


DefaultAssay(sal1hr) <- "SCT"
DefaultAssay(sal24hr) <- "SCT"
DefaultAssay(ptz1hr) <- "SCT"
DefaultAssay(ptz24hr) <- "SCT"


combo.list <- list("PTZ1Hour" = ptz1hr,"Saline1Hour" = sal1hr,"PTZ24Hour" = ptz24hr,"Saline24Hour" = sal24hr)

features <- SelectIntegrationFeatures(object.list = combo.list, nfeatures = 3000)

combo.list <- PrepSCTIntegration(object.list = combo.list, anchor.features = features)

combo.anchors <- FindIntegrationAnchors(object.list = combo.list, normalization.method = "SCT", anchor.features = features)

combo.all <- IntegrateData(anchorset = combo.anchors, normalization.method = "SCT")

saveRDS(combo.all, file = 'data/saved_objects/combo_all_df.rds')

#####Normalize ATAC-seq data#####

DefaultAssay(combo.all) <- "ATAC"

combo.all <- RunTFIDF(combo.all)
combo.all <- FindTopFeatures(combo.all, min.cutoff = 5)
combo.all <- RunSVD(combo.all)

#####Dim reduction and clustering#####

DefaultAssay(combo.all) <- "integrated"

# Run PCA

combo.all <- RunPCA(combo.all)

# Let's find neighbors, clusters, markers and identify our clusters

combo.all <- FindNeighbors(combo.all, dims = 1:30, verbose = T)

combo.all <- RunUMAP(combo.all, dims = 1:30, verbose = T)

combo.all <- FindClusters(combo.all,resolution = 0.2, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.3, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.4, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.5, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.6, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.7, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.8, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 0.9, verbose = T)
combo.all <- FindClusters(combo.all,resolution = 1, verbose = T)

saveRDS(combo.all, file = 'data/saved_objects/combo_all_df_clustered.rds')

