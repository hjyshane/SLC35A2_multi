# Load additional required libraries
library(SingleR)
library(celldex)

# ... [Previous code for data processing and integration remains the same] ...

# After creating the combined object and performing clustering:

# Load the mouse brain atlas data
mouseatlas.cut <- readRDS("~/PTZ_ATAC_scRNA_072024/File/mouseatlas_cut_processed.rds")

# Prepare the reference dataset
ref_data <- mouseatlas.cut@meta.data
ref_labels <- mouseatlas.cut$Taxonomy_group  # Adjust this based on the level of annotation you want to use

# Perform annotation using SingleR
singler_results <- SingleR(test = GetAssayData(combined, assay = "RNA", slot = "data"),
                           ref = ref_data,
                           labels = ref_labels,
                           method = "cluster",
                           clusters = combined@meta.data$seurat_clusters)

# Add SingleR labels to the Seurat object
combined$singler_labels <- singler_results$labels[match(combined@meta.data$seurat_clusters, rownames(singler_results))]

# Visualize the annotation results
p4 <- DimPlot(combined, reduction = "umap", group.by = "singler_labels", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types Annotated by SingleR")

# Combine with previous visualizations
print(p1 + p2 + p3 + p4)

# You can also create a table to summarize the cell type compositions
cell_type_summary <- table(combined$singler_labels, combined$orig.ident)
print(cell_type_summary)

# Visualize cell type proportions
cell_type_props <- prop.table(cell_type_summary, margin = 2)
cell_type_props_df <- as.data.frame(cell_type_props)
colnames(cell_type_props_df) <- c("CellType", "Group", "Proportion")

p5 <- ggplot(cell_type_props_df, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Type Proportions Across Groups",
       x = "Group", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# Save the updated combined object with annotations
saveRDS(combined, file = "combined_PTZ_SAL_rna_atac_object_annotated.rds")

# Differential expression analysis with cell type annotations
DefaultAssay(combined) <- "RNA"

# Function to perform DE analysis for each cell type
de_by_celltype <- function(seurat_obj, cell_type, group1, group2) {
  seurat_obj_subset <- subset(seurat_obj, subset = singler_labels == cell_type)
  FindMarkers(seurat_obj_subset, 
              ident.1 = group1, 
              ident.2 = group2, 
              group.by = "orig.ident",
              min.pct = 0.25, 
              logfc.threshold = 0.25)
}

# Get unique cell types
cell_types <- unique(combined$singler_labels)

# Perform DE analysis for each cell type
de_results_PTZ_vs_SAL <- lapply(cell_types, function(ct) {
  de_by_celltype(combined, ct, c("PTZ1hr", "PTZ24hr"), c("SAL1hr", "SAL24hr"))
})
names(de_results_PTZ_vs_SAL) <- cell_types

de_results_1hr_vs_24hr <- lapply(cell_types, function(ct) {
  de_by_celltype(combined, ct, c("PTZ1hr", "SAL1hr"), c("PTZ24hr", "SAL24hr"))
})
names(de_results_1hr_vs_24hr) <- cell_types

# Save DE results
for (ct in cell_types) {
  write.csv(de_results_PTZ_vs_SAL[[ct]], file = paste0("de_PTZ_vs_SAL_", ct, ".csv"))
  write.csv(de_results_1hr_vs_24hr[[ct]], file = paste0("de_1hr_vs_24hr_", ct, ".csv"))
}

# Differential accessibility analysis with cell type annotations
DefaultAssay(combined) <- "ATAC"

# Function to perform DA analysis for each cell type
da_by_celltype <- function(seurat_obj, cell_type, group1, group2) {
  seurat_obj_subset <- subset(seurat_obj, subset = singler_labels == cell_type)
  FindMarkers(seurat_obj_subset, 
              ident.1 = group1, 
              ident.2 = group2, 
              group.by = "orig.ident",
              min.pct = 0.25, 
              logfc.threshold = 0.25)
}

# Perform DA analysis for each cell type
da_results_PTZ_vs_SAL <- lapply(cell_types, function(ct) {
  da_by_celltype(combined, ct, c("PTZ1hr", "PTZ24hr"), c("SAL1hr", "SAL24hr"))
})
names(da_results_PTZ_vs_SAL) <- cell_types

da_results_1hr_vs_24hr <- lapply(cell_types, function(ct) {
  da_by_celltype(combined, ct, c("PTZ1hr", "SAL1hr"), c("PTZ24hr", "SAL24hr"))
})
names(da_results_1hr_vs_24hr) <- cell_types

# Save DA results
for (ct in cell_types) {
  write.csv(da_results_PTZ_vs_SAL[[ct]], file = paste0("da_PTZ_vs_SAL_", ct, ".csv"))
  write.csv(da_results_1hr_vs_24hr[[ct]], file = paste0("da_1hr_vs_24hr_", ct, ".csv"))
}
