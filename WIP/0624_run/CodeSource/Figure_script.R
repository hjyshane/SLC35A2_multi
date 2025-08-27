base::library(Seurat)
base::library(SeuratExtend)
base::library(ggplot2)
base::library(tidyverse)
base::library(qs)
base::library(patchwork)
base::library(SeuratExtend)

# object with all celltypes
cobj <- qs::qread("./qsave/combined_cell_type.qs")
DefaultAssay(cobj) <- "RNA"
cobj <- Seurat::NormalizeData(cobj)

# object with filtered with 50 cells or more
cobj_2 <- qs::qread("./qsave/filtered_50_cells_atac_cleaned.qs")
cobj_2 <- NormalizeData(cobj_2)

#### Figure 1C/D UMAP ####
# shortend cluster name
short <- c("ExcitatoryNeuronsCA1" = "ExN_CA1",
           "ExcitatoryNeuronsCA3" = "ExN_CA3",
           "ExcitatoryNeuronsMatureDG" = "ExN_mDG",
           "ExcitatoryNeuronsImmatureDG" = "ExN_iDG",
           "L56ExcitatoryNeurons" = "ExN_5/6",
           "L45ExcitatoryNeurons" = "ExN_4/5",
           "L23ExcitatoryNeurons" = "ExN_2/3",
           "LGE" = "InN_LGE",
           "CGE" = "InN_CGE",
           "MGE" = "InN_MGE",
           "NdnfRelnInhibitoryInterneurons" = "InN_Ndnf",
           "Lamp5positive" = "InN_Lamp5",
           "Astrocytes" = "Astrocytes",
           "Microglia" = "Microglia",
           "COPMFOL" = "COP/MFOL",
           "OPCs" = "OPCs")
cobj@meta.data$short_ct <- short[as.character(cobj$cell_type)]

# RNA UMAP by cluster
umap_by_ct <- SeuratExtend::DimPlot2(cobj, 
                                     reduction = "RNA_int_umap", 
                                     features = "short_ct", 
                                     label = TRUE, 
                                     label.size = 5,
                                     repel = TRUE,
                                     theme = list(
                                         Seurat::NoLegend(), 
                                         Seurat::NoAxes(), 
                                         ggtitle("RNA UMAP by Cell Type")
                                         )) +  
    theme_umap_arrows(x_label = "UMAP_1",
                      y_label = "UMAP_2") 

ggplot2::ggsave(umap_by_ct, file = "./Plots/UMAPbyCellType.png", width = 6, height = 6, dpi = 300)

# Multimodal UMAP by cluster
umap_by_mtct <- SeuratExtend::DimPlot2(cobj, 
                                       reduction = "Multi_UMAP", 
                                       features = "short_ct", 
                                       label = TRUE,
                                       label.size = 5,
                                       repel = TRUE,
                                       theme = list(
                                           Seurat::NoLegend(), 
                                           Seurat::NoAxes(), 
                                           ggtitle("Multimodal UMAP by Cell Type")
                                           )) +  
    theme_umap_arrows(x_label = "UMAP_1",
                      y_label = "UMAP_2") 


ggplot2::ggsave(umap_by_mtct, file = "./Plots/MultiUMAPbyCellType.png", width = 6, height = 6, dpi = 300)


#### Figure 1E Cell type markers ####
# Marker gene list
genes_to_plot <- c(
    "Aqp4",  # Astrocytes
    "P2ry12",  # Microglia
    "Pdgfra",  # OPCs
    "Mbp",  # COPMFOL (Myelin-associated)
    "Dcn",  # CA1
    "Matn2",  # CA1
    "Grik4",  # CA3
    "Prox1",  # DG in general CA1-
    "Calb1",  # ExcitatoryNeuronsMatureDG
    "Dcx",  # ExcitatoryNeuronsImmatureDG
    "Satb2",  # cortical neuron
    "Cux2",  # L234ExcitatoryNeurons
    "Bcl11b",  # L5ExcitatoryNeurons CA1+
    # "Rorb", # L4/5ExcitatoryNeurons
    "Etv1",  # L5ExcitatoryNeurons
    "Tbr1",  # L56ExcitatoryNeurons
    "Foxp2",  # L6
    "Gad2", # inhibitory neurons
    "Vip",  # CGE-derived cells
    "Adarb2",
    "Meis2",  # LGE-derived interneurons
    "Lhx6",
    "Sst",  # MGE-derived interneurons
    "Ndnf",  # NdnfRelnInhibitoryInterneurons
    "Lamp5",  # Lamp5 Positive
    "Col1a1",  # Vascular Cells
    "Vim"     # Marker for progenitor or stem-like cells
)

# Order the colunm
custom_order <- c(
    "Astrocytes", "Microglia", "OPCs", "COPMFOL",
    "ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG","ExcitatoryNeuronsImmatureDG",
    "L23ExcitatoryNeurons","L45ExcitatoryNeurons", "L56ExcitatoryNeurons",
    "CGE","LGE", "MGE",  "NdnfRelnInhibitoryInterneurons","Lamp5positive")

# order column
cobj$cell_type <-
    base::factor(cobj$cell_type,
           levels = custom_order)


Seurat::DefaultAssay(cobj) <- "RNA"
Seurat::Idents(cobj) <- "cell_type"

# DotPlot
dot_plot <- Seurat::DotPlot(
    cobj,
    features = genes_to_plot,
    cols = c("darkgrey", "indianred1"),
    dot.scale = 8) +
    ggplot2::theme(
        axis.text.x =
            ggplot2::element_text(
                angle = 90,
                vjust = 0.8,
                hjust = 1,
                size = 10
            ),
        axis.text.y =
            ggplot2::element_text(
                size = 10
            ),
        axis.title.x =
            ggplot2::element_blank(),
        axis.title.y =
            ggplot2::element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
    ) +
    ggplot2::coord_flip()

ggplot2::ggsave(dot_plot, filename = "./Plots/MarkersDotPlot.png", width = 7, height = 7, dpi = 300)

#### Figure 1F Gene activity score ####
Seurat::DefaultAssay(cobj) <- "ACTIVITY"
Seurat::Idents(cobj) <- "cell_type"

# DotPlot
dot_plot <- Seurat::DotPlot(
    cobj,
    scale.min = 40,
    col.min = 0,
    features = genes_to_plot,
    cols = c("darkgrey", "indianred1"),
    dot.scale = 8) +
    ggplot2::theme(
        axis.text.x =
            ggplot2::element_text(
                angle = 90,
                vjust = 0.8,
                hjust = 1,
                size = 10
            ),
        axis.text.y =
            ggplot2::element_text(
                size = 10
            ),
        axis.title.x =
            ggplot2::element_blank(),
        axis.title.y =
            ggplot2::element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
    ) +
    ggplot2::coord_flip()

ggplot2::ggsave(dot_plot, filename = "./Plots/ActivityMarkersDotPlot.png", width = 7, height = 7, dpi = 300)


#### Figure 1G IEG VlnPlots/DotPlot ####
Seurat::DefaultAssay(cobj) <- "RNA"
Seurat::Idents(cobj) <- "cell_type"

IEGs <- c("Jun", "Fos", "Junb", "Fosb", "Egr1", "Homer1", "Snap25", "Nr4a3")
Seurat::Idents(cobj) <- "Sample"
vp <- Seurat::VlnPlot(cobj,features = "Fos", cols = colorspace::rainbow_hcl(4))

for (gene in IEGs) {
    Seurat::Idents(cobj) <- "Sample"
    vp <- Seurat::VlnPlot(cobj,features = gene, cols = colorspace::rainbow_hcl(4)) + Seurat::NoLegend()
    ggsave(vp, file = paste0("./Plots/", gene, "_Sample_Vlnplot.png"), width = 8, height = 6, dpi = 300)
    
    Seurat::Idents(cobj) <- "cell_type"
    vp2 <- Seurat::VlnPlot(cobj,features = gene, cols = colorspace::rainbow_hcl(16)) + Seurat::NoLegend()
    ggsave(vp2, file = paste0("./Plots/", gene, "_CellType_Vlnplot.png"), width = 8, height = 6, dpi = 300)
}

for (sample in unique(cobj$Sample)) {
    ob <-  subset(cobj, subset = Sample == sample)
    Seurat::Idents(ob) <- "cell_type"
    
    for (gene in IEGs) {
    vp <- Seurat::VlnPlot(ob,features = gene) + Seurat::NoLegend()
    ggsave(vp, file = paste0("./Plots/", gene, "_", sample, "_CellType_Vlnplot.png"), width = 8, height = 6, dpi = 300)
    }
}

for (ct in unique(cobj$cell_type)) {
    ob <-  subset(cobj, subset = cell_type == ct)
    Seurat::Idents(ob) <- "Sample"
    
    for (gene in IEGs) {
        vp <- Seurat::VlnPlot(ob,features = gene) + NoLegend()
        ggsave(vp, file = paste0("./Plots/", gene, "_", ct, "_CellType_Vlnplot.png"), width = 8, height = 6, dpi = 300)
    }
}

# All IEGs in one DotPlot
IEGdot <- SeuratExtend::DotPlot2(
    cobj,
    group.by = "Sample",
    features = c("Jun", "Fos", "Junb", "Fosb", "Egr1", "Homer1", "Snap25", "Nr4a3"),
    border = F,
    color_scheme = "Reds") +
    ggplot2::theme(
        axis.text.x =
        ggplot2::element_text(
            size = 10),
        axis.text.y =
            ggplot2::element_text(
                size = 10),
        axis.title.x =
            ggplot2::element_blank(),
        axis.title.y =
            ggplot2::element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
ggplot2::ggsave(IEGdot, filename = "./Plots/IEGdot_celltyep.png", width = 8, height = 6,dpi = 300)

#### Figure 2A DEG stacked barplot ####
library(patchwork)

deg_dir <- "./DEG/sig_csv"
deg_files <- list.files(deg_dir, pattern = ".csv")
full_deg <- data.frame()
for (file in deg_files) {
    deg <- read.csv(file = file.path(deg_dir, file), row.names = "X")
    rownames(deg) <- NULL
    full_deg <- rbind(full_deg, deg)
}

write.csv(file = "./DEG/full_sig_deg.csv", full_deg, row.names = F)

full_deg <- read.csv("./DEG/full_sig_deg.csv")
hp_deg <- full_deg %>%
    filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))
hp_deg$cell_type <- factor(hp_deg$cell_type, levels = c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))
hp_deg$Regulation <- factor(hp_deg$Regulation, levels = c("Up", "Down"))

summary_total <- hp_deg %>%
    group_by(cell_type, Regulation) %>%
    summarize(count = n()) %>%
    mutate(shortend = case_when(cell_type == "ExcitatoryNeuronsCA1" ~ "CA1",
                                cell_type == "ExcitatoryNeuronsCA3" ~ "CA3",
                                cell_type == "ExcitatoryNeuronsMatureDG" ~ "DG",
                                cell_type == "ExcitatoryNeuronsImmatureDG" ~ "Immature DG",
                                TRUE ~ NA))

summary_table <- hp_deg %>%
    group_by(cell_type, comparison, Regulation) %>%
    summarize(count = n()) %>%
    mutate(shortend = case_when(cell_type == "ExcitatoryNeuronsCA1" ~ "CA1",
                                cell_type == "ExcitatoryNeuronsCA3" ~ "CA3",
                                cell_type == "ExcitatoryNeuronsMatureDG" ~ "DG",
                                cell_type == "ExcitatoryNeuronsImmatureDG" ~ "Immature DG",
                                TRUE ~ NA))

degcountplot <- ggplot(data = summary_total, aes(x = reorder(shortend, -count), y = count, fill = Regulation)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=count),position = position_stack(vjust = 0.6), color = "black",size = 3) +
    theme(title = element_text(size = 11),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,  size = 10)) +
    labs(title = "Number Differentially expressed gene in Hippocampus",
         y = "Number of Genes")

# separate plot for each comparison
comparisons <- c("PTZvsSAL_1hr", "PTZvsSAL_24hr", "24hrvs1hr_PTZ")
regulation_colors <- c("Up" = "coral", "Down" = "darkturquoise")
plot_list <- list()
for(comp in comparisons) {
    plot_data <- summary_table %>% filter(comparison == comp)
    
    p <- ggplot(data = plot_data, aes(x = shortend, y = count, fill = Regulation)) +
        geom_bar(stat="identity") +
        geom_text(aes(label=count), position = position_stack(vjust = 0.6), color = "black", size = 3) +
        theme_minimal() +
        ylim(0, 650) +
        labs(title = comp, y = "Number of Genes", x = "") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
              legend.position = "none", panel.grid = element_blank())  # Remove individual legends
    
    plot_list[[comp]] <- p
}

# Combine with shared legend
combined_plot <- wrap_plots(plot_list, ncol = 3) +
    plot_layout(guides = "collect") &  # Collect legends
    theme(legend.position = "bottom")   # Position shared legend


ggsave(degcountplot, file = "./Plots/deg_count_plot.png", width = 6, height = 8, dpi = 300)
ggsave(combined_plot, file = "./Plots/deg_count_combined.png", width = 6, height = 8, dpi = 300)

#### Figure 2B Vann diagram ####
full_deg <- read_csv("./DEG/full_sig_deg.csv")
comparisons <- unique(full_deg$comparison)
hpct <- c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG")

get_gene <- function(data, comparisons, celltypes) {
    result <- list()
    for (comp in comparisons) {
        result[[comp]] <- list()
        for (ct in celltypes) {
            result[[comp]][[ct]] <- data %>% filter(comparison == comp, cell_type == ct) %>% pull(gene)
        }
    }
    return (result)
}

gene_list <- get_gene(full_deg, comparisons, hpct)

library(VennDiagram)
library(ggvenn)
library(RColorBrewer)

ptzsal1hr <- list(
    "CA1" = gene_list$PTZvsSAL_1hr$ExcitatoryNeuronsCA1,
    "CA3" = gene_list$PTZvsSAL_1hr$ExcitatoryNeuronsCA3,
    "DG" = gene_list$PTZvsSAL_1hr$ExcitatoryNeuronsMatureDG
    )

ptzsal1hrplot <- ggvenn(ptzsal1hr, fill_color = c("orange", "lightblue", "lightgreen"),
       stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "PTZ vs SAL 1hr")

ptzsal24hr <- list(
    "CA1" = gene_list$PTZvsSAL_24hr$ExcitatoryNeuronsCA1,
    "CA3" = gene_list$PTZvsSAL_24hr$ExcitatoryNeuronsCA3,
    "DG" = gene_list$PTZvsSAL_24hr$ExcitatoryNeuronsMatureDG
)

ptzsal24hrplot <- ggvenn(ptzsal24hr, fill_color = c("orange", "lightblue", "lightgreen"),
                        stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "PTZ vs SAL 24hr")

ptztemp <- list(
    "CA1" = gene_list$`24hrvs1hr_PTZ`$ExcitatoryNeuronsCA1,
    "CA3" = gene_list$`24hrvs1hr_PTZ`$ExcitatoryNeuronsCA3,
    "DG" = gene_list$`24hrvs1hr_PTZ`$ExcitatoryNeuronsMatureDG
)

ptztempplot <- ggvenn(ptztemp, fill_color = c("orange", "lightblue", "lightgreen"),
                        stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "24hr vs 1hr PTZ")

saltemp <- list(
    "CA1" = gene_list$`24hrvs1hr_SAL`$ExcitatoryNeuronsCA1,
    "CA3" = gene_list$`24hrvs1hr_SAL`$ExcitatoryNeuronsCA3,
    "DG" = gene_list$`24hrvs1hr_SAL`$ExcitatoryNeuronsMatureDG
)

saltempplot <- ggvenn(saltemp, fill_color = c("orange", "lightblue", "lightgreen"),
                        stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "24hr vs 1hr SAL")

ggsave(ptzsal1hrplot, file = "./Plots/venn_ptzsal1hr.png", width = 6, height = 6, dpi = 300)
ggsave(ptzsal24hrplot, file = "./Plots/venn_ptzsal24hr.png", width = 6, height = 6, dpi = 300)
ggsave(ptztempplot, file = "./Plots/venn_ptztemp.png", width = 6, height = 6, dpi = 300)
ggsave(saltempplot, file = "./Plots/venn_saltemp.png", width = 6, height = 6, dpi = 300)

#### Figure 2D Gene behavior bar plots for HP ####
gb <- read.csv("./DEG/gene_behavior.csv") %>%
    filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))

# Create function to ensure all patterns are present
make_complete_data <- function(data, cell_type_name) {
    data %>%
        complete(pattern = c("Increase then Baseline",                                                                        
                             "Decrease then Baseline", 
                             "Late Downregulation",
                             "Late Upregulation",
                             "Decrease then Reversal", 
                             "Increase then Reversal", 
                             "Sustained Decrease", 
                             "Sustained Increase"), fill = list(n = 0)) %>%
        mutate(pattern = factor(pattern, levels = c("Sustained Increase", "Increase then Reversal", "Increase then Baseline",                                                                        
                                                    "Sustained Decrease", "Decrease then Reversal", "Decrease then Baseline", 
                                                    "Late Upregulation", "Late Downregulation")))  # Set consistent order
}

ca1_gv <- gb %>% filter(cell_type == "ExcitatoryNeuronsCA1", pattern %in% c("Sustained Increase", "Increase then Reversal", "Increase then Baseline",                                                                        
                                                                            "Sustained Decrease", "Decrease then Reversal", "Decrease then Baseline", 
                                                                            "Late Upregulation", "Late Downregulation")) %>%
    count(pattern, sort = T) %>%
    make_complete_data("CA1")

ca1_gv$pattern <- factor(ca1_gv$pattern, levels = c("Increase then Baseline",                                                                        
                                                    "Decrease then Baseline", 
                                                    "Late Downregulation",
                                                    "Late Upregulation",
                                                    "Decrease then Reversal", 
                                                    "Increase then Reversal", 
                                                    "Sustained Decrease", 
                                                    "Sustained Increase"))

ca1 <- ggplot(ca1_gv, mapping = aes(x = pattern, y = n, fill = pattern)) +
    geom_bar(stat = "identity") +
    ylim(0, 150) +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
    labs(x = NULL, y = "Number of Genes", title = "Gene behavior in CA1")

ca3_gv <- gb %>% filter(cell_type == "ExcitatoryNeuronsCA3", pattern %in% c("Sustained Increase", "Increase then Reversal", "Increase then Baseline",                                                                        
                                                                            "Sustained Decrease", "Decrease then Reversal", "Decrease then Baseline", 
                                                                            "Late Upregulation", "Late Downregulation")) %>%
    count(pattern, sort = T) %>%
    make_complete_data("CA3")

ca3_gv$pattern <- factor(ca3_gv$pattern, levels = c("Increase then Baseline",                                                                        
                                                    "Decrease then Baseline", 
                                                    "Late Downregulation",
                                                    "Late Upregulation",
                                                    "Decrease then Reversal", 
                                                    "Increase then Reversal", 
                                                    "Sustained Decrease", 
                                                    "Sustained Increase"))

ca3<- ggplot(ca3_gv, mapping = aes(x = pattern, y = n, fill = pattern)) +
    geom_bar(stat = "identity") +
    NoLegend() +
    ylim(0, 150) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
    labs(x = NULL, y = "Number of Genes", title = "Gene behavior in CA3")

dg_gv <- gb %>% filter(cell_type == "ExcitatoryNeuronsMatureDG", pattern %in% c("Sustained Increase", "Increase then Reversal", "Increase then Baseline",                                                                        
                                                                            "Sustained Decrease", "Decrease then Reversal", "Decrease then Baseline", 
                                                                            "Late Upregulation", "Late Downregulation")) %>%
    count(pattern, sort = T) %>%
    make_complete_data("DG")

dg_gv$pattern <- factor(dg_gv$pattern, levels = c("Increase then Baseline",                                                                        
                                                    "Decrease then Baseline", 
                                                    "Late Downregulation",
                                                    "Late Upregulation",
                                                    "Decrease then Reversal", 
                                                    "Increase then Reversal", 
                                                    "Sustained Decrease", 
                                                    "Sustained Increase"))

dg <- ggplot(dg_gv, mapping = aes(x = pattern, y = n, fill = pattern)) +
    geom_bar(stat = "identity") +
    NoLegend() +
    ylim(0, 150) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
    labs(x = NULL, y = "Number of Genes", title = "Gene behavior in Mature DG")

idg_gv <- gb %>% filter(cell_type == "ExcitatoryNeuronsImmatureDG", pattern %in% c("Sustained Increase", "Increase then Reversal", "Increase then Baseline",                                                                        
                                                                            "Sustained Decrease", "Decrease then Reversal", "Decrease then Baseline", 
                                                                            "Late Upregulation", "Late Downregulation")) %>%
    count(pattern, sort = T) %>%
    make_complete_data("ImmatureDG")

idg_gv$pattern <- factor(idg_gv$pattern, levels = c("Increase then Baseline",                                                                        
                                                   "Decrease then Baseline", 
                                                   "Late Downregulation",
                                                   "Late Upregulation",
                                                   "Decrease then Reversal", 
                                                   "Increase then Reversal", 
                                                   "Sustained Decrease", 
                                                   "Sustained Increase"))

idg <- ggplot(idg_gv, mapping = aes(x = pattern, y = n, fill = pattern)) +
    geom_bar(stat = "identity") +
    NoLegend() +
    ylim(0, 150) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid = element_blank()) +
    labs(x = NULL, y = "Number of Genes", title = "Gene behavior in Immature DG")

ggsave(ca1, file = "./Plots/ca1_behavior.png", width = 4, height = 4, dpi = 300)
ggsave(ca3, file = "./Plots/ca3_behavior.png", width = 4, height = 4, dpi = 300)
ggsave(dg, file = "./Plots/dg_behavior.png", width = 4, height = 4, dpi = 300)
ggsave(idg, file = "./Plots/idg_behavior.png", width = 4, height = 4, dpi = 300)

#### Figure 2E Gene behavior line graph example genes ####
#### Figure 3A DAG stacked barplot ####
library(patchwork)

dag_dir <- "./DAG/Sig_ATAC/"
dag_files <- list.files(dag_dir, pattern = ".csv")
full_dag <- data.frame()
for (file in dag_files) {
    dag <- read.csv(file = file.path(dag_dir, file))
    full_dag <- rbind(full_dag, dag)
}

write.csv(file = "./DAG/full_sig_dag.csv", full_deg)

full_dag <- read.csv("./DAG/full_sig_dag.csv")
hp_dag <- full_dag %>%
    filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))
hp_dag$cell_type <- factor(hp_dag$cell_type, levels = c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))
hp_dag$Regulation <- factor(hp_dag$Regulation, levels = c("Up", "Down"))

summary_total <- hp_dag %>%
    group_by(cell_type, Regulation) %>%
    summarize(count = n()) %>%
    mutate(shortend = case_when(cell_type == "ExcitatoryNeuronsCA1" ~ "CA1",
                                cell_type == "ExcitatoryNeuronsCA3" ~ "CA3",
                                cell_type == "ExcitatoryNeuronsMatureDG" ~ "DG",
                                cell_type == "ExcitatoryNeuronsImmatureDG" ~ "Immature DG",
                                TRUE ~ NA))

summary_table <- hp_dag %>%
    group_by(cell_type, comparison, Regulation) %>%
    summarize(count = n()) %>%
    mutate(shortend = case_when(cell_type == "ExcitatoryNeuronsCA1" ~ "CA1",
                                cell_type == "ExcitatoryNeuronsCA3" ~ "CA3",
                                cell_type == "ExcitatoryNeuronsMatureDG" ~ "DG",
                                cell_type == "ExcitatoryNeuronsImmatureDG" ~ "Immature DG",
                                TRUE ~ NA))

degcountplot <- ggplot(data = summary_total, aes(x = shortend, y = count, fill = Regulation)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=count),position = position_stack(vjust = 0.6), color = "black",size = 3) +
    theme(title = element_text(size = 11),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,  size = 10)) +
    labs(title = "DAG in Hippocampus",
         y = "Number of Peaks")

# separate plot for each comparison
comparisons <- c("PTZvsSAL_1hr", "PTZvsSAL_24hr", "24hrvs1hr_PTZ")
regulation_colors <- c("Up" = "coral", "Down" = "darkturquoise")
plot_list <- list()
for (comp in comparisons) {
    plot_data <- summary_table %>% filter(comparison == comp)
    
    p <- ggplot(plot_data, aes(x = shortend, y = count, fill = Regulation)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = count), position = position_stack(vjust = 0.6), color = "black", size = 3) +
        scale_fill_manual(values = regulation_colors, drop = FALSE) +  # fix colors
        theme_minimal() +
        ylim(0, 100) +
        labs(title = comp, y = "Number of Peaks", x = "") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
            legend.position = "none",
            panel.grid = element_blank()
        )
    
    plot_list[[comp]] <- p
}

# Combine with shared legend
combined_plot <- wrap_plots(plot_list, ncol = 3) +
    plot_layout(guides = "collect") &  # Collect legends
    theme(legend.position = "bottom")   # Position shared legend


ggsave(degcountplot, file = "./Plots/dag_count_plot.png", width = 6, height = 8, dpi = 300)
ggsave(combined_plot, file = "./Plots/dag_count_combined.png", width = 6, height = 8, dpi = 300)


#### Figure 3B Vann diagram ATAC ####
dag_files <- list.files("./DAG/Sig_ATAC/", full.names = T)
full_dag <- data.frame()
for (file in dag_files) {
    df <- read.csv(file)
    if (nrow(df) != 0){
    full_dag <- bind_rows(full_dag, df)
    }
}
write.csv(full_dag, file = "./DAG/full_sig_dag.csv")

comparisons <- unique(full_dag$comparison)
hpct <- c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG")

get_peak <- function(data, comparisons, celltypes) {
    result <- list()
    for (comp in comparisons) {
        result[[comp]] <- list()
        for (ct in celltypes) {
            result[[comp]][[ct]] <- data %>% filter(comparison == comp, cell_type == ct) %>% pull(ranges) # change closestgene <-> ranges for which you pick
        }
    }
    return (result)
}

peak_list <- get_peak(full_dag, comparisons, hpct)

library(VennDiagram)
library(ggvenn)
library(RColorBrewer)

ptzsal1hr <- list(
    "CA1" = peak_list$PTZvsSAL_1hr$ExcitatoryNeuronsCA1,
    "CA3" = peak_list$PTZvsSAL_1hr$ExcitatoryNeuronsCA3,
    "DG" = peak_list$PTZvsSAL_1hr$ExcitatoryNeuronsMatureDG
)

ptzsal1hrplot <- ggvenn(ptzsal1hr, fill_color = c("orange", "lightblue", "lightgreen"),
                        stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "PTZ vs SAL 1hr")

ptzsal24hr <- list(
    "CA1" = peak_list$PTZvsSAL_24hr$ExcitatoryNeuronsCA1,
    "CA3" = peak_list$PTZvsSAL_24hr$ExcitatoryNeuronsCA3,
    "DG" = peak_list$PTZvsSAL_24hr$ExcitatoryNeuronsMatureDG
)

ptzsal24hrplot <- ggvenn(ptzsal24hr, fill_color = c("orange", "lightblue", "lightgreen"),
                         stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "PTZ vs SAL 24hr")

ptztemp <- list(
    "CA1" = peak_list$`24hrvs1hr_PTZ`$ExcitatoryNeuronsCA1,
    "CA3" = peak_list$`24hrvs1hr_PTZ`$ExcitatoryNeuronsCA3,
    "DG" = peak_list$`24hrvs1hr_PTZ`$ExcitatoryNeuronsMatureDG
)

ptztempplot <- ggvenn(ptztemp, fill_color = c("orange", "lightblue", "lightgreen"),
                      stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "24hr vs 1hr PTZ")

saltemp <- list(
    "CA1" = peak_list$`24hrvs1hr_SAL`$ExcitatoryNeuronsCA1,
    "CA3" = peak_list$`24hrvs1hr_SAL`$ExcitatoryNeuronsCA3,
    "DG" = peak_list$`24hrvs1hr_SAL`$ExcitatoryNeuronsMatureDG
)

saltempplot <- ggvenn(saltemp, fill_color = c("orange", "lightblue", "lightgreen"),
                      stroke_size = 0.5) + 
    theme(title = element_text(size = 20)) +
    labs(title = "24hr vs 1hr SAL")

ggsave(ptzsal1hrplot, file = "./Plots/venn_ptzsal1hr_ATAC_peak.png", width = 6, height = 6, dpi = 300)
ggsave(ptzsal24hrplot, file = "./Plots/venn_ptzsal24hr_ATAC_peak.png", width = 6, height = 6, dpi = 300)
ggsave(ptztempplot, file = "./Plots/venn_ptztemp_ATAC_peak.png", width = 6, height = 6, dpi = 300)
ggsave(saltempplot, file = "./Plots/venn_saltemp_ATAC_peak.png", width = 6, height = 6, dpi = 300)
#### Supp UMAP per sample ####
# RNA UMAP by sample group
umap_by_group <- SeuratExtend::DimPlot2(cobj, 
                                        reduction = "RNA_int_umap", 
                                        features = "Sample", 
                                        cols = c("PTZ_1hr" = "pink",
                                                 "PTZ_24hr" = "darkred",
                                                 "SAL_1hr" = "lightblue",
                                                 "SAL_24hr" = "darkblue"),
                                        theme = list(
                                            Seurat::NoAxes(), 
                                            ggtitle("RNA UMAP by Sample"),
                                        )) +  
    theme_umap_arrows(x_label = "UMAP_1",
                      y_label = "UMAP_2") 

ggplot2::ggsave(umap_by_group, file = "./Plots/UMAPbySample.png", width = 6, height = 6, dpi = 300)

#### Supp QC plots ####
# RNA QC plot
DefaultAssay(cobj) <- "RNA"
Idents(cobj) <- "Sample"
rnaqc <- VlnPlot2(cobj, 
                 features = c( "nCount_RNA", "nFeature_RNA", "percent_mt_rna"),
                 cols = colorspace::rainbow_hcl(4),
                 pt = F,
                 box = F) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(rnaqc, file = "./Plots/RNA_vlnplot.png", width = 12, height = 4)


# ATAC QC plot
DefaultAssay(cobj) <- "ATAC"
Idents(cobj) <- "Sample"
atacqc <- VlnPlot2(cobj, 
                  features = c("nCount_ATAC", "nFeature_ATAC", "blacklist_fraction",  "nucleosome_signal", "TSS.enrichment", "pct_reads_in_peaks"),
                  cols = colorspace::rainbow_hcl(4),
                  pt = F,
                  box = F) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(atacqc, file = "./Plots/ATAC_vlnplot.png", width = 12, height = 8)


allqc <- VlnPlot2(cobj, 
                   features = c("nCount_RNA", "nFeature_RNA", "percent_mt_rna", "nCount_ATAC", "nFeature_ATAC", "blacklist_fraction",  "nucleosome_signal", "TSS.enrichment", "pct_reads_in_peaks"),
                   cols = colorspace::rainbow_hcl(4),
                   pt = F,
                   box = F) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(allqc, file = "./Plots/QC_vlnplot.png", width = 12, height = 8)

#### Overlapped/Unique genes/peaks ####
full_deg <- read.csv("./DEG/full_sig_deg.csv")

hp_deg <- full_deg %>% filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))

comparisons <- unique(full_deg$comparison)
hpct <- c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG")

get_gene <- function(data, comparisons, celltypes) {
    result <- list()
    for (comp in comparisons) {
        result[[comp]] <- list()
            for (ct in celltypes) {
            result[[comp]][[ct]] <- data %>% filter(comparison == comp, cell_type == ct) %>% pull(gene)
        }
    }
    return (result)
}

gene_list <- get_gene(full_deg, comparisons, hpct)

library(VennDetail)

# retrieve from Figure 2B
ptzsal1hr_obj <- VennDetail::venndetail(ptzsal1hr)
ptzsal24hr_obj <- VennDetail::venndetail(ptzsal24hr)
ptztemp_obj <- VennDetail::venndetail(ptztemp)
saltemp_obj <- VennDetail::venndetail(saltemp)

ptzsal1hr_obj <- as.data.frame(ptzsal1hr_obj@result) %>%  rename(gene = Detail, venn_group = Subset) %>% mutate(comparison = "PTZvsSAL_1hr")
ptzsal24hr_obj <- as.data.frame(ptzsal24hr_obj@result) %>%  rename(gene = Detail, venn_group = Subset) %>%  mutate(comparison = "PTZvsSAL_24hr")
ptztemp_obj <- as.data.frame(ptztemp_obj@result) %>% rename(gene = Detail, venn_group = Subset) %>%  mutate(comparison = "24hrvs1hr_PTZ")
saltemp_obj <- as.data.frame(saltemp_obj@result) %>% rename(gene = Detail, venn_group = Subset) %>% mutate(comparison = "24hrvs1hr_SAL")

venn_list <- bind_rows(ptzsal1hr_obj, ptzsal24hr_obj, ptztemp_obj, saltemp_obj)
venn_list_deg <- venn_list %>%
    left_join(hp_deg, by = c("gene", "comparison"))

write.csv(venn_list_deg, file = "./DEG/hp_deg_venn.csv")

# ATAC part (pick closest gene or range)
full_dag <- read.csv("./DAG/full_sig_dag.csv")
hp_dag <- full_dag %>% filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG"))

comparisons <- unique(full_dag$comparison)
hpct <- c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG", "ExcitatoryNeuronsImmatureDG")

get_peak <- function(data, comparisons, celltypes) {
    result <- list()
    for (comp in comparisons) {
        result[[comp]] <- list()
        for (ct in celltypes) {
            result[[comp]][[ct]] <- data %>% filter(comparison == comp, cell_type == ct) %>% pull(ranges) # change closestgene <-> ranges for which you pick
        }
    }
    return (result)
}

peak_list <- get_peak(full_dag, comparisons, hpct)

# retrieve from Figure 3B
ptzsal1hr_obj <- VennDetail::venndetail(ptzsal1hr)
ptzsal24hr_obj <- VennDetail::venndetail(ptzsal24hr)
ptztemp_obj <- VennDetail::venndetail(ptztemp)
saltemp_obj <- VennDetail::venndetail(saltemp)

# change closestgene <-> ranges for which you pick
ptzsal1hr_obj <- as.data.frame(ptzsal1hr_obj@result) %>%  rename(closestgene = Detail, venn_group = Subset) %>% mutate(comparison = "PTZvsSAL_1hr")
ptzsal24hr_obj <- as.data.frame(ptzsal24hr_obj@result) %>%  rename(closestgene = Detail, venn_group = Subset) %>%  mutate(comparison = "PTZvsSAL_24hr")
ptztemp_obj <- as.data.frame(ptztemp_obj@result) %>% rename(closestgene = Detail, venn_group = Subset) %>%  mutate(comparison = "24hrvs1hr_PTZ")
saltemp_obj <- as.data.frame(saltemp_obj@result) %>% rename(closestgene = Detail, venn_group = Subset) %>% mutate(comparison = "24hrvs1hr_SAL")

venn_list <- bind_rows(ptzsal1hr_obj, ptzsal24hr_obj, ptztemp_obj, saltemp_obj)
venn_list_dag <- venn_list %>%
    left_join(hp_dag, by = c("closestgene", "comparison"))

write.csv(venn_list_dag, file = "./DAG/hp_dag_venn_gene.csv")

#### GO with gene behavior ####
gb <- read.csv("./DEG/gene_behavior.csv")
counts_matrix <- GetAssayData(cobj_2, assay = "RNA", layer = "counts")
bg_gene <- rownames(counts_matrix)[Matrix::rowSums(counts_matrix > 0) > 10]

dir.create("./GO/Gene_behavior/csv", recursive = T)

hp_gb <- gb %>%
    filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG"))

go_list <- list()

for (i in unique(hp_gb$pattern)) {
    gb_set <- filter(hp_gb, pattern == i)
    go <- clusterProfiler::enrichGO(
        gene = gb_set$gene,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        universe = bg_gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05
        )
    
    go_df <- go@result %>%
        mutate(
            log10_p_adjust = -log10(p.adjust),
            GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
            BgRatio_numeric = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
            fold_enrichment = GeneRatio_numeric / BgRatio_numeric,
            Pvalue_base_enrichment_score = fold_enrichment * -log10(p.adjust)
        )
    if (nrow(go) != 0){
        top_terms <- go_df %>%
            arrange(p.adjust) %>% head(20)
            
        go_plot <- ggplot(top_terms, aes(
            x = reorder(Description, p.adjust),  
            y = p.adjust,
            fill = fold_enrichment,
            size = Count)) +
            geom_point(shape = 21) +  
            coord_flip() +  
            viridis::scale_fill_viridis(option = "magma") +  
            theme_minimal() +
            theme(axis.text.y = element_text(size = 20),
                  title = element_text(size = 22),
                  panel.grid = element_blank(),
                  panel.background = element_blank(),
                  legend.position = "bottom",
                  legend.title = element_text(size = 15),
                  plot.margin = margin(t = 15, r = 15, b = 15, l = 15), # More bottom space for legend
                 ) +
            labs(
                title = paste(i),
                x = "",
                y = "p.adjust",
                fill = "Fold Enrichment",
                size = "Gene Count")
        
        ggsave(file.path("./GO/Gene_behavior/Plots", paste0("GO_dotplot_", i, ".png")), 
               go_plot, width = 18, height = 8, dpi = 300)
        write.csv(go, file = paste0("./GO/Gene_behavior/csv/", i, "_GO.csv"))
        go_list[[i]] <- go_df}
}

#### GO with Venn group ####
venn_list_deg <- read.csv("./DEG/hp_deg_venn.csv")
go_files <- list.files("./GO/csv", full.names = T)

# for (go_file in go_files) {
#     # Extract filename info
#     file_info <- gsub(".csv", "", basename(go_file))
#     parts <- strsplit(gsub(".csv", "", basename(go_file)), "_")[[1]]
#     ct <- parts[1]
#     comp <- paste(parts[-1], collapse = "_")
#     
#     # Read GO results
#     go_data <- read_csv(go_file) %>% mutate(cell_type = ct,
#                                             comparison = comp)
#     
#     write.csv(go_data, file = go_file)
# }
    
go_files <- go_files[setdiff(6:19, c(14, 15))]
go_data <- purrr::map_dfr(go_files, read.csv)
go_data[[1]] <- NULL
go_data[[1]] <- NULL

go_data_10 <- go_data %>%
    group_by(cell_type, comparison) %>%
    slice_head(n = 10) %>%
    ungroup()

go_data_unwrapped <- go_data_10 %>%
    mutate(gene = strsplit(geneID, "/")) %>%
    unnest(gene)

go_data_unwrapped$geneID <- NULL

venn_list_deg_go <- venn_list_deg %>%
    inner_join(go_data_unwrapped, by = c("gene", "comparison", "cell_type"))

go_plot_data <- venn_list_deg_go %>%
    group_by(comparison, cell_type, venn_group, Description) %>%
    summarise(
        logp = -log10(min(p.adjust, na.rm = TRUE)),
        fold = min(Fold_enrichment),
        Count = min(Count),
        .groups = "drop"
    )


ggplot(go_plot_data, aes(x = reorder(Description, -logp), y = logp, fill = venn_group)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(cell_type ~ comparison, scales = "free_y", space = "free") +
    coord_flip() +
    theme_minimal() +
    theme(
        strip.text = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        legend.position = "bottom"
    ) +
    labs(
        x = "GO term",
        y = expression(-log[10](adjusted~p)),
        fill = "Venn Group",
        title = "GO Enrichment per Venn Group, Cell Type, and Comparison"
    )


split_go <- go_plot_data %>%
    group_split(comparison, cell_type)

split_keys <- go_plot_data %>%
    group_keys(comparison, cell_type) %>%
    transmute(name = paste(comparison, cell_type, sep = "_")) %>%
    pull(name)

names(split_go) <- split_keys

plot_list <- map(split_names, function(name) {
    data <- split_go[[name]]
    ggplot(data, aes(x = reorder(Description, -logp), y = logp, fill = venn_group)) +
        geom_bar(stat = "identity", position = "stack") +
        coord_flip() +
        labs(
            title = name,
            x = "GO Term",
            y = expression(-log[10](adjusted~p)),
            fill = "Venn Group"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 10),
            legend.position = "bottom"
        )
})


library(scatterpie)

go_pie_data <- venn_list_deg_go %>%
    mutate(logp = -log10(p.adjust)) %>%
    dplyr::select(Description, comparison, cell_type, venn_group, Count, Fold_enrichment, logp) %>%
    group_by(Description, comparison, cell_type, venn_group) %>%
    summarise(
        Count = sum(Count),
        Fold_enrichment = mean(Fold_enrichment),
        logp = max(logp),
        .groups = "drop"
    )

go_pie_data <- go_pie_data %>%
    group_by(comparison, cell_type) %>%
    mutate(x = row_number()) %>%
    ungroup()

ggplot(go_pie_data, aes(
    x = logp,
    y = Description,
    size = Count,
    fill = Fold_enrichment,
    color = venn_group
)) +
    geom_point(shape = 21, stroke = 1) +
    facet_grid(cell_type ~ comparison, scales = "free_y", space = "free_y")

#### Gene behavior of GO genes ####
go_files <- list.files("./GO/csv/", pattern = ".csv", full.names = TRUE)[6:19]

gb <- read.csv("./DEG/gene_behavior.csv")
hp_gb <- gb %>%
    filter(cell_type %in% c("ExcitatoryNeuronsCA1", "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsMatureDG")) %>%
    mutate(pattern_simplified = case_when(pattern %in% c("Increase then Baseline", "Decrease then Baseline") ~ "Acture Response",
                                          pattern %in% c("Late Downregulation", "Late Upregulation") ~ "Delayed Response",
                                          pattern %in% c("Decrease then Reversal", "Increase then Reversal") ~ "Dynamic Response",
                                          pattern %in% c("Sustained Decrease", "Sustained Increase") ~ "Sustained Response",
                                          pattern %in% c("Increase in PTZ", "Decrease in PTZ") ~ "PTZ only",
                                          pattern %in% c("No change in PTZvsSAL") ~ "No respone",
                                          TRUE ~ "Other"))
hp_patterns <- unique(hp_gb$pattern)

results = list()

for (go_file in go_files) {
    # Read GO results
    go_data <- read_csv(go_file) %>% arrange(p.adjust) %>% head(20)

    # Extract filename info
    file_info <- gsub(".csv", "", basename(go_file))
    parts <- strsplit(gsub(".csv", "", basename(go_file)), "_")[[1]]
    ct <- parts[1]
    comp <- paste(parts[-1], collapse = "_")
    
    go_deg <- data.frame()
    
    for(i in 1:nrow(go_data)) {
    term_name <- go_data$Description[i]
    term_genes <- unlist(strsplit(go_data$geneID[i], "/"))
    
    # Find patterns for these genes
    selected_genes <- hp_gb %>%
            filter(cell_type == ct) %>%
            filter(gene %in% term_genes) %>%
            select(gene, cell_type,  pattern_simplified) %>%
            drop_na() %>%
            mutate(GO_term = term_name,
                   comparison = comp,
                   GO_pvalue = go_data$p.adjust[i],
                   file = file_info)
        
    go_deg <- bind_rows(go_deg, selected_genes)
    }
    
    go_plot_data <- go_deg %>%
    group_by(GO_term, pattern_simplified, GO_pvalue) %>%
    summarise(gene_count = n(), .groups = "drop") %>%
    mutate(logp = -log10(GO_pvalue))


    plot <- ggplot(data = go_plot_data, aes( x = reorder(GO_term, -logp), y = gene_count, fill = pattern_simplified)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip() + 
    theme_minimal() +
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          axis.text.y = element_text(size = 15),
          title = element_text(size = 15),
          panel.background = element_blank(),
          legend.title = element_text(size = 15),
          ) + # More bottom space for legend)
        labs(x = "",
         y = "Number of Genes",
         fill = "Pattern",
         title = paste(ct, comp),
         )
    
    ggsave(plot = plot, filename = paste0("./GO/Gene_behavior/Plots/", ct, "_", comp, ".png"), width = 16, height = 8, dpi = 300)
    
    if (nrow(go_deg) !=0){
    results[[file_info]] = go_deg}
    
}


#### Module Score ####
DefaultAssay(cobj) <- "RNA"
cobj <- NormalizeData(cobj)

# Get top 20 genes enriched 
ca1enriched <- FindMarkers(cobj, ident.1 = "ExcitatoryNeuronsCA1", verbose = FALSE) %>%
    arrange(-avg_log2FC) %>%
    rownames_to_column(var = "gene") %>%
    pull(gene) %>% 
    .[1:20]

cobj <- AddModuleScore(cobj,
                      features = list(ca1enriched),
                      name="CA1_enriched")

# Plot scores
ca1enrichedplot <- FeaturePlot(cobj,
            features = "CA1_enriched11", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

ca3enriched <- FindMarkers(cobj, ident.1 = "ExcitatoryNeuronsCA3", verbose = FALSE) %>%
    arrange(-avg_log2FC) %>%
    rownames_to_column(var = "gene") %>%
    pull(gene) %>% 
    .[1:20]

cobj <- AddModuleScore(cobj,
                       features = list(ca3enriched),
                       name="CA3_enriched")

# Plot scores
ca3enrichedplot <- FeaturePlot(cobj,
            features = "CA3_enriched1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

dgenriched <- FindMarkers(cobj, ident.1 = "ExcitatoryNeuronsMatureDG", verbose = FALSE) %>%
    arrange(-avg_log2FC) %>%
    rownames_to_column(var = "gene") %>%
    pull(gene) %>% 
    .[1:20]

cobj <- AddModuleScore(cobj,
                       features = list(dgenriched),
                       name="DG_enriched")

# Plot scores
dgenrichedplot <- FeaturePlot(cobj,
            features = "DG_enriched1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

idgenriched <- FindMarkers(cobj, ident.1 = "ExcitatoryNeuronsImmatureDG", verbose = FALSE) %>%
    arrange(-avg_log2FC) %>%
    rownames_to_column(var = "gene") %>%
    pull(gene) %>% 
    .[1:20]

cobj <- AddModuleScore(cobj,
                       features = list(idgenriched),
                       name="IDG_enriched")

# Plot scores
idgenrichedplot <- FeaturePlot(cobj,
            features = "IDG_enriched1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

ggsave(ca1enrichedplot, file = "./Plots/CA1_ModuleScore.png", width = 8, height = 8, dpi = 300)
ggsave(ca3enrichedplot, file = "./Plots/CA3_ModuleScore.png", width = 8, height = 8, dpi = 300)
ggsave(dgenrichedplot, file = "./Plots/DG_ModuleScore.png", width = 8, height = 8, dpi = 300)
ggsave(idgenrichedplot, file = "./Plots/IDG_ModuleScore.png", width = 8, height = 8, dpi = 300)
