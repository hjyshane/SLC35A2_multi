#single-cell multiome analysis

library(Seurat)
library(Signac)
library(ggvenn)

combined <- readRDS("/igm/home/tab013/multiome/combined.rds")

DimPlot(combined, reduction = "umap", label = T) + NoLegend()
DimPlot(combined, reduction = "umap", label = F, split.by = "sample") + NoLegend()

combined <- PrepSCTFindMarkers(combined)
markers <- FindAllMarkers(combined)
write.csv(markers, "/igm/home/tab013/multiome/combined_markers.csv")

combined.markers <- read.csv("/igm/home/tab013/multiome/combined_markers.csv")
top5 <- combined.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
DoHeatmap(combined, features = top5$gene, label = F) + theme(axis.text.y = element_text(size = 6))


DefaultAssay(combined) <- 'SCT'
FeaturePlot(combined, reduction = "umap", features = c("Calb", "Cacng5", "Fibcd1", "Dkk3"))
FeaturePlot(combined, reduction = "umap", features = c("Pdzd2", "Tox3", "Dcn", "Wfs1"))
FeaturePlot(combined, reduction = "umap", features = c("Cbin4", "Nr4a2", "Rxfp1", "Fermt1"))

FeaturePlot(combined, reduction = "umap", features = c("Mag", "Pdgfra"))
FeaturePlot(combined, reduction = "umap", features = c("Gad1", "Lhx6", "Adarb2"))
FeaturePlot(combined, reduction = "umap", features = c("Slc17a7", "Tle4", "Prox1", "Bcl11b"))
FeaturePlot(combined, reduction = "umap", features = c("Slc1a2", "P2ry12", "Reln"))
FeaturePlot(combined, reduction = "umap", features = c("Pax6", "Eomes"))


combined@meta.data$celltype_sample <- paste(combined@meta.data$celltypes, combined@meta.data$sample, sep = "_")
DotPlot(combined, features = c("Fos", "Jun", "Egr1", "Arc"), group.by = "sample")

sub <- subset(combined, celltypes == "unknown", invert = T)
DotPlot(sub, features = c("Fos", "Jun", "Egr1", "Arc"), group.by = "celltype_sample")

# # # # # # # # # # ## Differential Expression ## # # # # # # # # # # #
#Generate DEG between conditions for each individual cluster -- creates a single file 

##Append each cluster with sample
{
  DefaultAssay(combined) <- "SCT"
  combined $celltype.sample <- paste(Idents(combined), combined $sample, sep = "_")
  combined $celltype <- Idents(combined)
  Idents(combined) <- "celltype.sample"
  #saveRDS(combined, "/igm/home/tab013/combined_sample.rds")
}

print(unique(combined$celltype.sample))

#Define clusters to pass in ident.1 and ident.2 
{
  sample.class <- levels(combined @active.ident)
  PTZ_1hr <- sort(sample.class[grep("_PTZ_1hr" , sample.class)])
  SAL_1hr <- sort(sample.class[grep("_Saline_1hr", sample.class)])
}

#Check each of these lists and make sure they are the same items and in the same order. If not, it might mean a particular cluster is not present in one of the samples and needs to be removed. Items must match up between the two lists for this to work.

PTZ_1hr
SAL_1hr
#SAL_1hr <- SAL_1hr[-c(22)]

#Open new data frames to store temp individual and binded outputs
{
  response <- data.frame()
  all.response <- data.frame()
}

#set working directory to save all output files to
setwd("/igm/home/tab013/multiome/")

#Run loop to pass above defined clusters as variables ident.1 and ident.2
for(i in 1:length(PTZ_1hr)){
  skip_to_next <- FALSE
  response <-  FindMarkers(combined, assay = "SCT", slot = "data", ident.1 = PTZ_1hr[i], ident.2 = SAL_1hr[i], min.pct = 0.25, logfc.threshold = 0.0, return.threshold = 1.01)
  response[,"Cluster Comparison"] <- paste(PTZ_1hr[i],"vs", SAL_1hr[i], sep = " ", collapse = NULL) 
  response[, "original_gene_names"] <- row.names(x = response)
  tryCatch(
    if( i == 1 )
    {
      all.response <- response
    }
    else
    {
      all.response <- rbind(all.response, response)
    }
    , error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next) { next } 
  
  write.csv(all.response, "PTZ_1hr_vs_SAL_1hr.csv")
}

response.sig <- subset(all.response, p_val_adj < 0.05)
table(response.sig$`Cluster Comparison`)

# # # # # # # # # # ## Differential Expression ## # # # # # # # # # # #
#Generate DEG between conditions for each individual cluster -- creates a single file 

##Append each cluster with sample
{
  DefaultAssay(combined) <- "SCT"
  combined $celltype.sample <- paste(Idents(combined), combined $sample, sep = "_")
  combined $celltype <- Idents(combined)
  Idents(combined) <- "celltype.sample"
  #saveRDS(combined, "/igm/home/tab013/combined_sample.rds")
}

print(unique(combined$celltype.sample))

#Define clusters to pass in ident.1 and ident.2 
{
  sample.class <- levels(combined @active.ident)
  PTZ_1hr <- sort(sample.class[grep("_PTZ_1hr" , sample.class)])
  PTZ_24hr <- sort(sample.class[grep("_PTZ_24hr", sample.class)])
}

#Check each of these lists and make sure they are the same items and in the same order. If not, it might mean a particular cluster is not present in one of the samples and needs to be removed. Items must match up between the two lists for this to work.

PTZ_1hr
PTZ_24hr
#SAL_1hr <- SAL_1hr[-c(22)]

#Open new data frames to store temp individual and binded outputs
{
  response <- data.frame()
  all.response <- data.frame()
}

#set working directory to save all output files to
setwd("/igm/home/tab013/multiome/")

#Run loop to pass above defined clusters as variables ident.1 and ident.2
for(i in 1:length(PTZ_1hr)){
  skip_to_next <- FALSE
  response <-  FindMarkers(combined, assay = "SCT", slot = "data", ident.1 = PTZ_1hr[i], ident.2 = PTZ_24hr[i], min.pct = 0.25, logfc.threshold = 0.0, return.threshold = 1.01)
  response[,"Cluster Comparison"] <- paste(PTZ_1hr[i],"vs", PTZ_24hr[i], sep = " ", collapse = NULL) 
  response[, "original_gene_names"] <- row.names(x = response)
  tryCatch(
    if( i == 1 )
    {
      all.response <- response
    }
    else
    {
      all.response <- rbind(all.response, response)
    }
    , error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next) { next } 
  
  write.csv(all.response, "PTZ_1hr_vs_PTZ_24hr.csv")
}

response.sig <- subset(all.response, p_val_adj < 0.05 & avg_log2FC > abs(0.25))
table(response.sig$`Cluster Comparison`)

# # # # # # # # # # ## Differential Expression ## # # # # # # # # # # #
#Generate DEG between conditions for each individual cluster -- creates a single file 

##Append each cluster with sample
{
  DefaultAssay(combined) <- "SCT"
  combined $celltype.sample <- paste(Idents(combined), combined $sample, sep = "_")
  combined $celltype <- Idents(combined)
  Idents(combined) <- "celltype.sample"
  #saveRDS(combined, "/igm/home/tab013/combined_sample.rds")
}

print(unique(combined$celltype.sample))

#Define clusters to pass in ident.1 and ident.2 
{
  sample.class <- levels(combined @active.ident)
  PTZ_24hr <- sort(sample.class[grep("_PTZ_24hr" , sample.class)])
  SAL_24hr <- sort(sample.class[grep("_Saline_24hr", sample.class)])
}

#Check each of these lists and make sure they are the same items and in the same order. If not, it might mean a particular cluster is not present in one of the samples and needs to be removed. Items must match up between the two lists for this to work.

PTZ_24hr
SAL_24hr
#SAL_1hr <- SAL_1hr[-c(22)]

#Open new data frames to store temp individual and binded outputs
{
  response <- data.frame()
  all.response <- data.frame()
}

#set working directory to save all output files to
setwd("/igm/home/tab013/multiome/")

#Run loop to pass above defined clusters as variables ident.1 and ident.2
for(i in 1:length(PTZ_24hr)){
  skip_to_next <- FALSE
  response <-  FindMarkers(combined, assay = "SCT", slot = "data", ident.1 = PTZ_24hr[i], ident.2 = SAL_24hr[i], min.pct = 0.25, logfc.threshold = 0.0, return.threshold = 1.01)
  response[,"Cluster Comparison"] <- paste(PTZ_24hr[i],"vs", SAL_24hr[i], sep = " ", collapse = NULL) 
  response[, "original_gene_names"] <- row.names(x = response)
  tryCatch(
    if( i == 1 )
    {
      all.response <- response
    }
    else
    {
      all.response <- rbind(all.response, response)
    }
    , error = function(e) {skip_to_next <<- TRUE})
  
  if(skip_to_next) { next } 
  
  write.csv(all.response, "PTZ_24hr_vs_SAL_24hr.csv")
}

response.sig <- subset(all.response, p_val_adj < 0.05 & avg_log2FC > abs(0.25))
table(response.sig$`Cluster Comparison`)

#--------------------------------------------------------
ptz1h_sal1h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_SAL_1hr.csv")
ptz1h_sal1h <- subset(ptz1h_sal1h, p_val_adj < 0.05)
ptz1h_sal1h <- subset(ptz1h_sal1h, avg_log2FC > 0.25 | avg_log2FC < -0.25)
ptz1h_sal1h.degs <- as.data.frame(table(ptz1h_sal1h$Cluster.Comparison))

ptz24h_sal24h <- read.csv("/igm/home/tab013/multiome/PTZ_24hr_vs_SAL_24hr.csv")
ptz24h_sal24h <- subset(ptz24h_sal24h, p_val_adj < 0.05)
ptz24h_sal24h <- subset(ptz24h_sal24h, avg_log2FC > 0.25 | avg_log2FC < -0.25)
ptz24h_sal24h.degs <- as.data.frame(table(ptz24h_sal24h$Cluster.Comparison))

ptz1h_ptz24h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_PTZ_24hr.csv")
ptz1h_ptz24h <- subset(ptz1h_ptz24h, p_val_adj < 0.05)
ptz1h_ptz24h <- subset(ptz1h_ptz24h, avg_log2FC > 0.25 | avg_log2FC < -0.25)
ptz1h_ptz24h.degs <- as.data.frame(table(ptz1h_ptz24h$Cluster.Comparison))

#granule neurons
venn.data <-
  list(
    "PTZ1hr vs PTZ24hr" = ptz1h_ptz24h %>% filter(Cluster.Comparison == "Granule Neurons_PTZ_1hr vs Granule Neurons_PTZ_24hr") %>% pull(original_gene_names),
    "PTZ1hr vs SAL1hr" = ptz1h_sal1h %>% filter(Cluster.Comparison == "Granule Neurons_PTZ_1hr vs Granule Neurons_Saline_1hr") %>% pull(original_gene_names),
    "PTZ24hr vs SAL24hr" = ptz24h_sal24h %>% filter(Cluster.Comparison == "Granule Neurons_PTZ_24hr vs Granule Neurons_Saline_24hr") %>% pull(original_gene_names)
  )

ggvenn(venn.data)

#granule neuroblasts
venn.data <-
  list(
    "PTZ1hr vs PTZ24hr" = ptz1h_ptz24h %>% filter(Cluster.Comparison == "Granuale Neuroblasts_PTZ_1hr vs Granuale Neuroblasts_PTZ_24hr") %>% pull(original_gene_names),
    "PTZ1hr vs SAL1hr" = ptz1h_sal1h %>% filter(Cluster.Comparison == "Granuale Neuroblasts_PTZ_1hr vs Granuale Neuroblasts_Saline_1hr") %>% pull(original_gene_names),
    "PTZ24hr vs SAL24hr" = ptz24h_sal24h %>% filter(Cluster.Comparison == "Granuale Neuroblasts_PTZ_24hr vs Granuale Neuroblasts_Saline_24hr") %>% pull(original_gene_names)
  )

ggvenn(venn.data)

#excitatory neurons CA1
venn.data <-
  list(
    "PTZ1hr vs PTZ24hr" = ptz1h_ptz24h %>% filter(Cluster.Comparison == "Excitatory Neurons (CA1)_PTZ_1hr vs Excitatory Neurons (CA1)_PTZ_24hr") %>% pull(original_gene_names),
    "PTZ1hr vs SAL1hr" = ptz1h_sal1h %>% filter(Cluster.Comparison == "Excitatory Neurons (CA1)_PTZ_1hr vs Excitatory Neurons (CA1)_Saline_1hr") %>% pull(original_gene_names),
    "PTZ24hr vs SAL24hr" = ptz24h_sal24h %>% filter(Cluster.Comparison == "Excitatory Neurons (CA1)_PTZ_24hr vs Excitatory Neurons (CA1)_Saline_24hr") %>% pull(original_gene_names)
  )

ggvenn(venn.data)

#excitatory neurons 2
venn.data <-
  list(
    "PTZ1hr vs PTZ24hr" = ptz1h_ptz24h %>% filter(Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_PTZ_24hr") %>% pull(original_gene_names),
    "PTZ1hr vs SAL1hr" = ptz1h_sal1h %>% filter(Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_Saline_1hr") %>% pull(original_gene_names),
    "PTZ24hr vs SAL24hr" = ptz24h_sal24h %>% filter(Cluster.Comparison == "Excitatory Neurons 2_PTZ_24hr vs Excitatory Neurons 2_Saline_24hr") %>% pull(original_gene_names)
  )

ggvenn(venn.data)

#overlap table
ptz1h_ptz24h <- subset(ptz1h_ptz24h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_PTZ_24hr")
ptz1h_sal1h <- subset(ptz1h_sal1h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_Saline_1hr")
ptz24h_sal24h <- subset(ptz24h_sal24h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_24hr vs Excitatory Neurons 2_Saline_24hr")

overlap <- full_join(ptz1h_ptz24h, ptz1h_sal1h, by = "original_gene_names")
overlap <- full_join(overlap, ptz24h_sal24h, by = "original_gene_names")
write.csv(overlap, file = "/igm/home/tab013/multiome/excneuron2_degoverlap.csv")

##### DEG Plots #######
library(EnhancedVolcano)
ptz1h_sal1h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_SAL_1hr.csv")
ptz1h_sal1h_granuleneurons <- subset(ptz1h_sal1h, Cluster.Comparison == "Granule Neurons_PTZ_1hr vs Granule Neurons_Saline_1hr")
EnhancedVolcano(ptz1h_sal1h_granuleneurons,
                lab = ptz1h_sal1h_granuleneurons$original_gene_names,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                labSize = 4,
                FCcutoff = 0.25,
                pCutoff = .01,
                ylab = "-Log10(adjusted P-value)",
                title = "Granule Neurons",
                subtitle = "PTZ 1hr vs. Saline 1hr")

ptz1h_sal1h_excneurons2 <- subset(ptz1h_sal1h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_Saline_1hr")
EnhancedVolcano(ptz1h_sal1h_excneurons2,
                lab = ptz1h_sal1h_excneurons2$original_gene_names,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                labSize = 4,
                FCcutoff = 0.25,
                pCutoff = .01,
                ylab = "-Log10(adjusted P-value)",
                title = "Excitatory Neurons 2",
                subtitle = "PTZ 1hr vs. Saline 1hr")

ptz1h_sal1h_granuleneuroblasts <- subset(ptz1h_sal1h, Cluster.Comparison == "Granuale Neuroblasts_PTZ_1hr vs Granuale Neuroblasts_Saline_1hr")
EnhancedVolcano(ptz1h_sal1h_granuleneuroblasts,
                lab = ptz1h_sal1h_granuleneuroblasts$original_gene_names,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                labSize = 4,
                FCcutoff = 0.25,
                pCutoff = .01,
                ylab = "-Log10(adjusted P-value)",
                title = "Granule Neuroblasts",
                subtitle = "PTZ 1hr vs. Saline 1hr")
#24hr
ptz1h_ptz24h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_PTZ_24hr.csv")
ptz1h_ptz24h_granuleneuroblasts <- subset(ptz1h_ptz24h, Cluster.Comparison == "Granuale Neuroblasts_PTZ_1hr vs Granuale Neuroblasts_PTZ_24hr")
EnhancedVolcano(ptz1h_ptz24h_granuleneuroblasts,
                lab = ptz1h_ptz24h_granuleneuroblasts$original_gene_names,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                labSize = 4,
                FCcutoff = 0.25,
                pCutoff = .01,
                ylab = "-Log10(adjusted P-value)",
                title = "Granule Neuroblasts",
                subtitle = "PTZ 1hr vs. PTZ 24hr")

ptz1h_ptz24h_excneurons2 <- subset(ptz1h_ptz24h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_PTZ_24hr")
EnhancedVolcano(ptz1h_ptz24h_excneurons2,
                lab = ptz1h_ptz24h_excneurons2$original_gene_names,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                labSize = 4,
                FCcutoff = 0.25,
                pCutoff = .01,
                ylab = "-Log10(adjusted P-value)",
                title = "Excitatory Neurons 2",
                subtitle = "PTZ 1hr vs. PTZ 24hr")

##############################################################################################################
ptz1h_sal1h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_SAL_1hr.csv")
ptz1h_sal1h_granuleneurons_universe <- subset(ptz1h_sal1h, Cluster.Comparison == "Granule Neurons_PTZ_1hr vs Granule Neurons_Saline_1hr")

ptz1h_sal1h <- subset(ptz1h_sal1h, p_val_adj < 0.05)
ptz1h_sal1h <- subset(ptz1h_sal1h, avg_log2FC > 0.25 | avg_log2FC < -0.25)
ptz1h_sal1h_granuleneurons <- subset(ptz1h_sal1h, Cluster.Comparison == "Granule Neurons_PTZ_1hr vs Granule Neurons_Saline_1hr")

library(clusterProfiler)
library(enrichplot)
library(AnnotationHub)
go1 <- enrichGO(
  gene = ptz1h_sal1h_granuleneurons %>% pull(original_gene_names),
  universe      = ptz1h_sal1h_granuleneurons_universe %>% pull(original_gene_names),
  OrgDb         = "org.Mm.eg.db",
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.05
)
head(go1@result)
barplot(go1@result)

##
ptz1h_sal1h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_SAL_1hr.csv")
ptz1h_sal1h_excneurons2_universe <- subset(ptz1h_sal1h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_Saline_1hr")

ptz1h_sal1h <- subset(ptz1h_sal1h, p_val_adj < 0.05)
ptz1h_sal1h <- subset(ptz1h_sal1h, avg_log2FC > 0.25 | avg_log2FC < -0.25)
ptz1h_sal1h_excneurons2 <- subset(ptz1h_sal1h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_Saline_1hr")

library(clusterProfiler)
library(enrichplot)
library(AnnotationHub)
go1 <- enrichGO(
  gene = ptz1h_sal1h_excneurons2 %>% pull(original_gene_names),
  universe      = ptz1h_sal1h_excneurons2_universe %>% pull(original_gene_names),
  OrgDb         = "org.Mm.eg.db",
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.05
)
head(go1@result)

##
ptz1h_ptz24h <- read.csv("/igm/home/tab013/multiome/PTZ_1hr_vs_PTZ_24hr.csv")
ptz1h_ptz24h_excneurons2_universe <- subset(ptz1h_ptz24h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_PTZ_24hr")

ptz1h_ptz24h <- subset(ptz1h_ptz24h, p_val_adj < 0.05)
ptz1h_ptz24h <- subset(ptz1h_ptz24h, avg_log2FC > 0.25 | avg_log2FC < -0.25)
ptz1h_ptz24h_excneurons2 <- subset(ptz1h_ptz24h, Cluster.Comparison == "Excitatory Neurons 2_PTZ_1hr vs Excitatory Neurons 2_PTZ_24hr")

go1 <- enrichGO(
  gene = ptz1h_ptz24h_excneurons2 %>% pull(original_gene_names),
  universe      = ptz1h_ptz24h_excneurons2_universe %>% pull(original_gene_names),
  OrgDb         = "org.Mm.eg.db",
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "fdr",
  qvalueCutoff  = 0.05
)
head(go1@result)
barplot(go1)
