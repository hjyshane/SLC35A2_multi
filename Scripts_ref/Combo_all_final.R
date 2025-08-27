library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(tidyverse)
library(clustree)
library(plyr)
library(dplyr)
library(Nebulosa)
library(cowplot)

##### Import combined data set #####

combo.all <- readRDS("/igm/projects/Singlecell_analysis_projects/PTZ_mice/saved_objects/combo_all_df_clustered.rds")

# import reference
mb.cut <- readRDS("/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File/mouseatlas_cut_hippocampus.rds")

# clone query object for further analysis
combo.all.typed <- combo.all

##### ref and query #####

DefaultAssay(mb.cut) <- "SCT"

DefaultAssay(combo.all.typed) <- "SCT"

all.anchors <- FindTransferAnchors(reference = mb.cut, 
                                   query = combo.all.typed,
                                   # recompute.residuals = FALSE,
                                   normalization.method = 'SCT',
                                   dims = 1:30)

predictions <- TransferData(anchorset = all.anchors, 
                            refdata = mb.cut$Description,
                            dims = 1:30)
combo.all.typed <- AddMetaData(combo.all.typed, metadata = predictions)
combo.all.typed$Description <-  combo.all.typed$predicted.id

predictions <- TransferData(anchorset = all.anchors, 
                            refdata = mb.cut$TaxonomyRank1,
                            dims = 1:30)
combo.all.typed <- AddMetaData(combo.all.typed, metadata = predictions)
combo.all.typed$TaxonomyRank1 <-  combo.all.typed$predicted.id

predictions <- TransferData(anchorset = all.anchors, 
                            refdata = mb.cut$TaxonomyRank2,
                            dims = 1:30)
combo.all.typed <- AddMetaData(combo.all.typed, metadata = predictions)
combo.all.typed$TaxonomyRank2 <-  combo.all.typed$predicted.id

predictions <- TransferData(anchorset = all.anchors, 
                            refdata = mb.cut$TaxonomyRank3,
                            dims = 1:30)
combo.all.typed <- AddMetaData(combo.all.typed, metadata = predictions)
combo.all.typed$TaxonomyRank3 <-  combo.all.typed$predicted.id

predictions <- TransferData(anchorset = all.anchors, 
                            refdata = mb.cut$TaxonomyRank4,
                            dims = 1:30)
combo.all.typed <- AddMetaData(combo.all.typed, metadata = predictions)
combo.all.typed$TaxonomyRank4 <-  combo.all.typed$predicted.id



Idents(combo.all.typed) <- combo.all.typed$TaxonomyRank1
Idents(combo.all.typed) <- combo.all.typed$TaxonomyRank2
Idents(combo.all.typed) <- combo.all.typed$TaxonomyRank3
Idents(combo.all.typed) <- combo.all.typed$TaxonomyRank4
Idents(combo.all.typed) <- combo.all.typed$Description
# Verify the changes
cell_type_table <- table(combo.all.typed$predicted.id)
write.csv(cell_type_table, file = "/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/File/Merge_way/Jesse_cell_type_table.csv", row.names = FALSE)


#####Rename Idents#####

combo.all.typed <- RenameIdents(combo.all.typed,
                                'Dentate gyrus granule neurons'='Granule neurons',
                                'Non-glutamatergic neuroblasts'='Granule neuroblasts',
                                'Perivascular macrophages'='Microglia',
                                'Oligodendrocyte precursor cells'='OPCs',
                                'Vascular and leptomeningeal cells'='Vascular cells',
                                'Vascular endothelial cells'='Vascular cells',
                                'Ependymal cells'='Vascular cells',
                                'Pericytes'='Vascular cells',
                                'Telencephalon inhibitory interneurons'='Inhibitory interneurons',
                                'Telencephalon projecting excitatory neurons'='Excitatory neurons',
                                'Di- and mesencephalon excitatory neurons'='Cajal–Retzius cells')


combo.all.typed$celltype_MH <- Idents(combo.all.typed)

DimPlot(combo.all.typed, label=T, repel = TRUE, reduction = "umap") + NoLegend()
DimPlot(combo.all.typed, reduction = "umap")
DimPlot(combo.all.typed, label=T, repel = TRUE, reduction = "umap") + NoLegend()
DimPlot(combo.all.typed,  reduction = "umap")

#####Multimodal Neighbors#####

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight

DefaultAssay(combo.all.typed) <- "integrated"

combo.all.typed <- FindMultiModalNeighbors(
  object = combo.all.typed,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

#use the generated wknn graph to find clusters
#combo.all <- FindClusters(combo.all,  resolution = 0.4 ,graph.name = "wknn")

#build umap visualization using the neighbors found
combo.all.typed <- RunUMAP(
  object = combo.all.typed,
  nn.name = "weighted.nn",
  assay = "SCT",
  verbose = TRUE
)

levels(Idents(combo.all.typed))

DimPlot(combo.all.typed, reduction = 'umap')

highlight <- WhichCells(combo.all.typed, idents = c("Neuronal intermidate progenitor cells"))

DimPlot(combo.all.typed, reduction = 'umap', cells.highlight = highlight)

levels(Idents(mb.cut))


##### features #####

DefaultAssay(combo.all.typed) <- "ACT"
DefaultAssay(combo.all.typed) <- "SCT"

FeaturePlot(combo.all.typed, features = c("Auts2"))
FeaturePlot(combo.all.typed, features = c("Ptprc"))
plot_density(combo.all.typed, features = c("Auts2"))
FeaturePlot(combo.all.typed, features = c("Tmem119"))

FeaturePlot(combo.all.typed, features = c("Sox2"))
plot_density(combo.all.typed, features = c("Sox2"),pal = 'magma')
FeaturePlot(combo.all.typed, features = c("Satb1"))
plot_density(combo.all.typed, features = c("Satb1"),pal = 'magma')


##### gene activity #####
DefaultAssay(combo.all.typed) <- "ATAC"

combo.all.typed <- RunTFIDF(combo.all.typed)
combo.all.typed <- FindTopFeatures(combo.all.typed, min.cutoff = 'q0')
combo.all.typed <- RunSVD(combo.all.typed)

gene.activities <- GeneActivity(combo.all.typed)

combo.all.typed[['ACT']] <- CreateAssayObject(counts = gene.activities)

combo.all.typed <- NormalizeData(
  object = combo.all.typed,
  assay = 'ACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(combo.all.typed$nCount_ACT)
)

DefaultAssay(combo.all.typed) <- "ACT"

saveRDS(combo.all.typed, file = "data/saved_objects/combo_MH.rds")

FeaturePlot(combo.all.typed, features = c("Auts2"))

"Ptprc" %in% rownames(combo.all.typed)
"Auts2" %in% rownames(combo.all.typed)

length(rownames(combo.all.typed))

##### coverage #####

DefaultAssay(combo.all.typed) <- "ATAC"

CoveragePlot(
  object = combo.all.typed,
  region = "chr5-131435683-132544343",
  extend.upstream = 5000,
  extend.downstream = 5000,
  idents = c("Granule neurons", "Granule neuroblasts", "OPCs")
)

FeaturePlot(combo.all.typed, features = c("chr5-131435684-132544343"))

##### peaks as features #####

t <- as.list(rownames(combo.all.typed))

grep("chr5-13254", rownames(combo.all.typed),value = T)

FeaturePlot(combo.all.typed, features = "chr5-132543014-132543925")

plot_density(combo.all.typed, features = "chr5-132543014-132543925")


##### Make them plots #####

p1 <- DimPlot(combo.all.typed, 
               label=T, 
              # label.box = T,
               repel = TRUE,
              # pt.size = .5,
              reduction = "umap") +labs(title = "Hippocampal Region Cell types") + theme(
                plot.title = element_text(hjust = 0.5),axis.line=element_blank(),axis.text.x=element_blank(),
                axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),plot.background=element_blank()) + NoLegend()


DefaultAssay(combo.all.typed) <- "SCT"

p2 <- plot_density(combo.all.typed,
             features = "Auts2",
             pal = "plasma",
             #label=T, 
             # label.box = T,
             #repel = TRUE,
             # pt.size = .5,
             reduction = "umap") +labs(title = "Auts2 Gene Expression") + theme(
               plot.title = element_text(hjust = 0.5),axis.line=element_blank(),axis.text.x=element_blank(),
                                         axis.text.y=element_blank(),axis.ticks=element_blank(),
                                         axis.title.x=element_blank(),
                                         axis.title.y=element_blank(),
                                         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                         panel.grid.minor=element_blank(),plot.background=element_blank()) + NoLegend()

DefaultAssay(combo.all.typed) <- "ATAC"

p3 <- plot_density(combo.all.typed,
             features = "chr5-132543014-132543925",
             pal = "magma",
             #label=T, 
             # label.box = T,
             #repel = TRUE,
             # pt.size = .5,
             reduction = "umap") +labs(title = "Auts2 ATAC") + theme(
               plot.title = element_text(hjust = 0.5),axis.line=element_blank(),axis.text.x=element_blank(),
               axis.text.y=element_blank(),axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),plot.background=element_blank()) + NoLegend()

plot_grid(p1,p2,p3,nrow = 1,align = "h")



