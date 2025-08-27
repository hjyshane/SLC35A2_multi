#### GRN ####
library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(GENIE3)
library(data.table)
library(doParallel)
library(foreach)
library(igraph)

set.seed(42)

base_dir <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
setwd(base_dir)
cobj <- qs::qread("./qsave/filtered_50_cells_atac_cleaned.qs")

# DEG file
deg_dir <- "./DEG/sig_csv/"
if(!dir.exists(file.path(base_dir, "GRN"))) {dir.create(file.path(base_dir, "GRN"))}
  
deg_files <- list.files(deg_dir, pattern = ".csv", full.names = T)
neu_deg <- deg_files[9:36]

# Extract expression matrix from seurat
DefaultAssay(cobj) <- "RNA"
cobj <- NormalizeData(cobj, normalization.method = "LogNormalize", verbose = TRUE, slot = "data")
expr_matrix <- GetAssayData(cobj, assay = "RNA", layer = "data")

# TF list
# download.file("https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt", destfile = "./DEG/tflist.txt")
tflist <- read.table("./DEG/tflist.txt", col.names = "Gene") %>% as.data.frame()

for (file in neu_deg) {
  tryCatch({
  deg_list <- read.csv(file)
  message(basename(file))
  
  # extract gene list
  gene_list <- unique(deg_list$gene)

  # subset deg gene from matrix
  deg_matrix <- expr_matrix[rownames(expr_matrix) %in% gene_list, ] %>% as.matrix()

    # get tf that in deg list
  tf_filtered <- intersect(tflist$Gene, rownames(deg_matrix))

  # run GENIE3
  message("running GENIE3 ", basename(file))
  weightMarix <- GENIE3(deg_matrix, regulators = tf_filtered, nCores = 16)
  message("done GENIE3 ", basename(file))
  # melt
  melted <- reshape2::melt(weightMarix, value.name = "Score")
  melted <- dplyr::rename(melted, "target" = "Var1", "source" = "Var2", "weight" = "Score")
  
  # save
  write.csv(melted, file = file.path(base_dir, "GRN", basename(file)))
  message("saved ", basename(file))
  }, error = function(e) { 
    message("error in ", basename(file))}
  ) 
}


list <- list.files("./GRN/", full.names = T)
full <- data.frame()

for (file in list) {
  data <- read_csv(file) %>% mutate(file = basename(file))
  full <- rbind(full, data)
}

full <- full %>%
  arrange(desc(weight)) %>%
  filter(weight > 0.5)

write.csv(full, file = "./GRN/filtered_full.csv")



grn_graph <- igraph::make_full_graph(full)
igraph::plot.igraph(grn_graph)

grn_graph <- graph_from_data_frame(melted, directed = T)
igraph::plot.igraph(grn_graph)