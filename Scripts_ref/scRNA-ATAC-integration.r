# Load required libraries
library(Seurat)
library(Signac)
library(GenomicFeatures)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(BSgenome.Mmusculus.UCSC.mm10)  # Mouse genome
library(EnsDb.Mmusculus.v79)  # Mouse gene annotations

# Load custom GTF file (if you're using one)
custom_gtf <- rtracklayer::import("path/to/your/custom_mouse.gtf")
gene.coords <- genes(custom_gtf)
seqlevelsStyle(gene.coords) <- "UCSC"
genome(gene.coords) <- "mm10"  # Mouse genome build

# If not using a custom GTF, you can use the built-in mouse annotations:
# gene.coords <- genes(EnsDb.Mmusculus.v79)
# seqlevelsStyle(gene.coords) <- "UCSC"

# Process RNA data for PTZ1hr
rna_PTZ1hr <- CreateSeuratObject(counts = counts_rna_PTZ1hr, project = "PTZ1hr", min.cells = 3, min.features = 200)
rna_PTZ1hr[["percent.mt"]] <- PercentageFeatureSet(rna_PTZ1hr, pattern = "^mt-")
rna_PTZ1hr[["percent.rb"]] <- PercentageFeatureSet(rna_PTZ1hr, pattern = "^Rp[sl]")

# Plot QC metrics for PTZ1hr
p1 <- VlnPlot(rna_PTZ1hr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
p2 <- FeatureScatter(rna_PTZ1hr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p1 + p2)

# Filter cells for PTZ1hr
rna_PTZ1hr <- subset(rna_PTZ1hr, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15 & percent.rb > 1 & percent.rb < 50)

# Normalize and scale data for PTZ1hr
rna_PTZ1hr <- NormalizeData(rna_PTZ1hr)
rna_PTZ1hr <- FindVariableFeatures(rna_PTZ1hr)
rna_PTZ1hr <- ScaleData(rna_PTZ1hr)

# Perform PCA and UMAP for PTZ1hr
rna_PTZ1hr <- RunPCA(rna_PTZ1hr)
rna_PTZ1hr <- RunUMAP(rna_PTZ1hr, dims = 1:30)

# Doublet detection for PTZ1hr
sweep.res.list_PTZ1hr <- paramSweep_v3(rna_PTZ1hr, PCs = 1:30, sct = FALSE)
sweep.stats_PTZ1hr <- summarizeSweep(sweep.res.list_PTZ1hr, GT = FALSE)
bcmvn_PTZ1hr <- find.pK(sweep.stats_PTZ1hr)
optimal_pk_PTZ1hr <- bcmvn_PTZ1hr$pK[which.max(bcmvn_PTZ1hr$BCmetric)] %>% as.character() %>% as.numeric()
nExp_poi_PTZ1hr <- round(0.075*nrow(rna_PTZ1hr@meta.data))
rna_PTZ1hr <- doubletFinder_v3(rna_PTZ1hr, PCs = 1:30, pN = 0.25, pK = optimal_pk_PTZ1hr, nExp = nExp_poi_PTZ1hr, reuse.pANN = FALSE, sct = FALSE)
rna_PTZ1hr <- subset(rna_PTZ1hr, cells = rownames(rna_PTZ1hr@meta.data)[rna_PTZ1hr@meta.data[, grep("DF.classifications", colnames(rna_PTZ1hr@meta.data))] == "Singlet"])

# Repeat the above steps for PTZ24hr, SAL1hr, and SAL24hr
# [Code for PTZ24hr]
# [Code for SAL1hr]
# [Code for SAL24hr]

# Process ATAC data for PTZ1hr
chrom_assay_PTZ1hr <- CreateChromatinAssay(
  counts = counts_atac_PTZ1hr,
  sep = c(":", "-"),
  fragments = fragments_PTZ1hr,
  annotation = gene.coords
)

atac_PTZ1hr <- CreateSeuratObject(
  counts = chrom_assay_PTZ1hr,
  assay = "ATAC",
  project = "PTZ1hr"
)

atac_PTZ1hr <- AddModelMatrix(
  object = atac_PTZ1hr,
  model = ~ 0 + atac_PTZ1hr$orig.ident
)

atac_PTZ1hr <- NucleosomeSignal(atac_PTZ1hr)
atac_PTZ1hr <- TSSEnrichment(atac_PTZ1hr, fast = FALSE)

atac_PTZ1hr$blacklist_ratio <- atac_PTZ1hr$ATAC_blacklist_region_fragments / atac_PTZ1hr$ATAC_peak_region_fragments
atac_PTZ1hr$frac_reads_in_peaks <- atac_PTZ1hr$ATAC_peak_region_fragments / atac_PTZ1hr$ATAC_fragments

# Plot QC metrics for ATAC PTZ1hr
p3 <- VlnPlot(
  object = atac_PTZ1hr,
  features = c("nCount_ATAC", "TSS.enrichment", "blacklist_ratio", "nucleosome_signal", "frac_reads_in_peaks"),
  pt.size = 0.1,
  ncol = 5
)
print(p3)

# Filter cells for ATAC PTZ1hr
atac_PTZ1hr <- subset(
  x = atac_PTZ1hr,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

# Normalize and reduce dimensions for ATAC PTZ1hr
atac_PTZ1hr <- RunTFIDF(atac_PTZ1hr)
atac_PTZ1hr <- FindTopFeatures(atac_PTZ1hr, min.cutoff = 'q0')
atac_PTZ1hr <- RunSVD(atac_PTZ1hr)
atac_PTZ1hr <- RunUMAP(atac_PTZ1hr, reduction = 'lsi', dims = 2:30)

# Repeat the above ATAC processing steps for PTZ24hr, SAL1hr, and SAL24hr
# [Code for PTZ24hr]
# [Code for SAL1hr]
# [Code for SAL24hr]

# Integrate RNA data
rna_list <- list(rna_PTZ1hr, rna_PTZ24hr, rna_SAL1hr, rna_SAL24hr)
rna_features <- SelectIntegrationFeatures(object.list = rna_list)
rna_anchors <- FindIntegrationAnchors(object.list = rna_list, anchor.features = rna_features)
rna_integrated <- IntegrateData(anchorset = rna_anchors)

# Process integrated RNA data
DefaultAssay(rna_integrated) <- "integrated"
rna_integrated <- ScaleData(rna_integrated)
rna_integrated <- RunPCA(rna_integrated)
rna_integrated <- RunUMAP(rna_integrated, dims = 1:30)
rna_integrated <- FindNeighbors(rna_integrated, dims = 1:30)
rna_integrated <- FindClusters(rna_integrated, resolution = 0.5)

# Integrate ATAC data
atac_list <- list(atac_PTZ1hr, atac_PTZ24hr, atac_SAL1hr, atac_SAL24hr)
atac_features <- SelectIntegrationFeatures(object.list = atac_list, assay = "ATAC")
atac_anchors <- FindIntegrationAnchors(object.list = atac_list, anchor.features = atac_features, assay = "ATAC")
atac_integrated <- IntegrateData(anchorset = atac_anchors)

# Process integrated ATAC data
DefaultAssay(atac_integrated) <- "integrated"
atac_integrated <- RunTFIDF(atac_integrated)
atac_integrated <- FindTopFeatures(atac_integrated, min.cutoff = 'q0')
atac_integrated <- RunSVD(atac_integrated)
atac_integrated <- RunUMAP(atac_integrated, reduction = 'lsi', dims = 2:30)

# Integrate RNA and ATAC data
combined <- CreateSeuratObject(
  counts = rna_integrated[["RNA"]]$counts,
  assay = "RNA"
)

# Add ATAC assay
combined[["ATAC"]] <- atac_integrated[["ATAC"]]

# Normalize RNA data
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

# Dimension reduction and clustering
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# Visualize results
p1 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p3 <- FeaturePlot(combined, features = c("Your_Mouse_Gene_Of_Interest1", "Your_Mouse_Gene_Of_Interest2"))

print(p1 + p2 + p3)

# Save the combined object
saveRDS(combined, file = "combined_PTZ_SAL_rna_atac_object.rds")

# Differential expression analysis
DefaultAssay(combined) <- "RNA"
markers_PTZ_vs_SAL <- FindMarkers(combined, ident.1 = c("PTZ1hr", "PTZ24hr"), ident.2 = c("SAL1hr", "SAL24hr"), min.pct = 0.25, logfc.threshold = 0.25)
markers_1hr_vs_24hr <- FindMarkers(combined, ident.1 = c("PTZ1hr", "SAL1hr"), ident.2 = c("PTZ24hr", "SAL24hr"), min.pct = 0.25, logfc.threshold = 0.25)

# Save marker genes
write.csv(markers_PTZ_vs_SAL, file = "markers_PTZ_vs_SAL.csv")
write.csv(markers_1hr_vs_24hr, file = "markers_1hr_vs_24hr.csv")

# Differential accessibility analysis
DefaultAssay(combined) <- "ATAC"
da_regions_PTZ_vs_SAL <- FindMarkers(combined, ident.1 = c("PTZ1hr", "PTZ24hr"), ident.2 = c("SAL1hr", "SAL24hr"), min.pct = 0.25, logfc.threshold = 0.25)
da_regions_1hr_vs_24hr <- FindMarkers(combined, ident.1 = c("PTZ1hr", "SAL1hr"), ident.2 = c("PTZ24hr", "SAL24hr"), min.pct = 0.25, logfc.threshold = 0.25)

# Save differentially accessible regions
write.csv(da_regions_PTZ_vs_SAL, file = "da_regions_PTZ_vs_SAL.csv")
write.csv(da_regions_1hr_vs_24hr, file = "da_regions_1hr_vs_24hr.csv")
