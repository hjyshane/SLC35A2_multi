

#load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)


combo.1hr.typed$celltype <- Idents(combo.1hr.typed)

combo.1hr.typed$sample_type <- paste0(combo.1hr.typed$celltype, "_", combo.1hr.typed$sample)

#####DE between conditions (SCT)#####

DefaultAssay(combo.1hr.typed) <- "SCT"

Idents(combo.1hr.typed) <- combo.1hr.typed$sample_type

for (i in 1:length(levels(combo.1hr.typed$celltype))) {
  
  DEmarkers <-  FindMarkers(combo.1hr.typed, 
                            ident.1 = paste0(levels(combo.1hr.typed$celltype)[i],"_PTZ_1hr"),
                            ident.2 = paste0(levels(combo.1hr.typed$celltype)[i],"_Saline_1hr"),
                            min.pct = 0.3, 
                            logfc.threshold = 0.25,)
  
  DEmarkers.sig <- filter(DEmarkers, p_val_adj <= 0.05)
  
  DEmarkers.sig[6] <- rownames(DEmarkers.sig)
  
  names(DEmarkers.sig)[6] <- "gene"
  
  DEmarkers.sig[7] <- levels(combo.1hr.typed$celltype)[i]
  
  names(DEmarkers.sig)[7] <- "cell_type"
  
  
  if (nrow(DEmarkers.sig) > 0) {
    
    DElist <- DEmarkers.sig
    
    rownames(DElist) <- c()
    
  }
  
  if (i == 1) {
    
    DElist.all <- DElist
    
  } else {
    
    DElist.add <- DElist
    
    DElist.all <- bind_rows(DElist.all, DElist.add)
    
  }
  
}

saveRDS(DElist.all, "data/diff_lists/DE_genes_by_condition_1hr.rds")
write_excel_csv(DElist.all, "data/diff_lists/DE_genes_by_condition_1hr.csv")