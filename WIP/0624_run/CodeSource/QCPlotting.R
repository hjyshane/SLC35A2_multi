#### directory ####
qc_dir <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/QC/"
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
  dir.create(file.path(qc_dir, "Plots"), recursive = TRUE)
  }

# read object
qcobject <- qs::qread("~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/combined_cell_type.qs")

#### QC fueature plots ####
mt_rna <- FeaturePlot(object = qcobject,
                      features = "percent_mt_rna",
                      reduction = "Multi_UMAP",  # Use UMAP or tSNE for dimensionality reduction
                      label = F, # Label clusters
                      cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                      pt.size = 0.2               # Adjust point size for better visibility
                      )

nCount <- FeaturePlot(object = qcobject,
                      features = "nCount_RNA",
                      reduction = "Multi_UMAP", # Use UMAP or tSNE for dimensionality reduction
                      label = F, # Label clusters
                      cols = c("lightgrey", "blue"),# Adjust colors to highlight high/low values
                      pt.size = 0.2               # Adjust point size for better visibility
                      )

nFeature <- FeaturePlot(object = qcobject,
                        features = "nFeature_RNA",
                        reduction = "Multi_UMAP", # Use UMAP or tSNE for dimensionality reduction
                        label = F, # Label clusters
                        cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                        pt.size = 0.2               # Adjust point size for better visibility
                        )

nCount_atac <- FeaturePlot(object = qcobject,
                           features = "nCount_ATAC",
                           reduction = "Multi_UMAP", # Use UMAP or tSNE for dimensionality reduction
                           label = F, # Label clusters
                           cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                           pt.size = 0.2               # Adjust point size for better visibility
                           )

nFeature_atac <- FeaturePlot(object = qcobject,
                           features = "nFeature_ATAC",
                           reduction = "Multi_UMAP", # Use UMAP or tSNE for dimensionality reduction
                           label = F, # Label clusters
                           cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                           pt.size = 0.2               # Adjust point size for better visibility
                           )

TSSscore <- FeaturePlot(object = qcobject,
                        features = "TSS.enrichment",
                        reduction = "Multi_UMAP", # Use UMAP or tSNE for dimensionality reduction
                        label = F, # Label clusters
                        cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                        pt.size = 0.2               # Adjust point size for better visibility
                        )


NSSenrichment <- FeaturePlot(object = qcobject,
                             features = "nucleosome_signal",
                             reduction = "Multi_UMAP", # Use UMAP or tSNE for dimensionality reduction
                             label = F, # Label clusters
                             cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                             pt.size = 0.2               # Adjust point size for better visibility
                             )

ATAC_QC <- DensityScatter(qcobject, 
                          x = 'nCount_ATAC', 
                          y = 'TSS.enrichment', 
                          log_x = TRUE, 
                          quantiles = TRUE
                          )
# combined
RNAqcplot <- mt_rna | nCount | nFeature
ATACqcplot <- nFeature_atac / nCount_atac | TSSscore / NSSenrichment | ATAC_QC

# save
ggsave(mt_rna, filename = file.path(qc_dir, "Plots", "mt_rna.png"), width = 8, height = 8, dpi = 300)
ggsave(nCount, filename = file.path(qc_dir, "Plots", "nCount.png"), width = 8, height = 8, dpi = 300)
ggsave(nFeature, filename = file.path(qc_dir, "Plots", "nFeature.png"), width = 8, height = 8, dpi = 300)
ggsave(nCount_atac, filename = file.path(qc_dir, "Plots", "nCount_ATAC.png"), width = 8, height = 8, dpi = 300)
ggsave(nFeature_atac, filename = file.path(qc_dir, "Plots", "nFeature_ATAC.png"), width = 8, height = 8, dpi = 300)
ggsave(TSSscore, filename = file.path(qc_dir, "Plots", "TSSscore.png"), width = 8, height = 8, dpi = 300)
ggsave(NSSenrichment, filename = file.path(qc_dir, "Plots", "NSSenrichment.png"), width = 8, height = 8, dpi = 300)
ggsave(ATAC_QC, filename = file.path(qc_dir, "Plots", "ATAC_QC.png"), width = 8, height = 8, dpi = 300)
ggsave(RNAqcplot, filename = file.path(qc_dir, "Plots", "RNAqcplot.png"), width = 16, height = 8, dpi = 300)
ggsave(ATACqcplot, filename = file.path(qc_dir, "Plots", "ATACqcplot.png"), width = 16, height = 8, dpi = 300)

#### matrics ####
qc_metrics <- c( "nCount_RNA", "nFeature_RNA", "percent_mt_rna", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment")
qc_thresholds <- list(c(1000, 25000), c(200, 7000), c(NA, 5), c(1000, 100000), c(500, 20000), c(2, NA))

# get cell counts for each matrics
qc_pass <- data.frame(
  RNA_count_pass = sum(qcobject$nCount_RNA >= 1000 & qcobject$nCount_RNA <= 25000),
  RNA_feature_pass = sum(qcobject$nFeature_RNA >= 200),
  MT_percent_pass = sum(qcobject$percent_mt_rna <= 5),
  ATAC_count_pass = sum(qcobject$nCount_ATAC >= 1000 & qcobject$nCount_ATAC <= 100000),
  ATAC_feature_pass = sum(qcobject$nFeature_ATAC >= 500 & qcobject$nFeature_ATAC <= 100000),
  TSS_enrichment_pass = sum(qcobject$TSS.enrichment >= 2))

# Calculate percentages
qc_pass_percent <- round(qc_pass / ncol(qcobject) * 100, 2)
qc_pass_table <- rbind(qc_pass, qc_pass_percent)

# Calculate percentages
qc_pass_percent <- round(qc_pass / ncol(qcobject) * 100, 2)
qc_pass_table <- rbind(qc_pass, qc_pass_percent)

summary <-
  data.frame(summary(qcobject@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent_mt_rna", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment")]))

write.csv(summary, file = paste0(qc_dir, "QC_matrics.csv"))

#### Distribution ####
p1 <- ggplot(qcobject@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mt_rna)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "RNA metrics",
       x = "UMI count",
       y = "Gene count",
       color = "Percent mitochondrial")

p2 <- ggplot(qcobject@meta.data, aes(x = orig.ident, y = nCount_RNA)) +
  geom_violin(aes(fill = orig.ident), show.legend = FALSE) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_minimal() +
  labs(title = "RNA count distribution",
       x = "Sample",
       y = "UMI count")

p3 <- ggplot(qcobject@meta.data, aes(x = orig.ident, y = percent_mt_rna)) +
  geom_violin(aes(fill = orig.ident), show.legend = FALSE) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_minimal() +
  labs(title = "Mitochondrial percentage",
       x = "Sample",
       y = "Percent mitochondrial")

# ATAC QC metrics
p4 <-ggplot(qcobject@meta.data, aes(x = nCount_ATAC, y = nFeature_ATAC, color = TSS.enrichment)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "ATAC metrics",
       x = "ATAC count",
       y = "Peak count",
       color = "TSS enrichment")

p5 <-ggplot(qcobject@meta.data, aes(x = orig.ident, y = nCount_ATAC)) +
  geom_violin(aes(fill = orig.ident), show.legend = FALSE) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_minimal() +
  labs(title = "ATAC count distribution",
       x = "Sample",
       y = "ATAC count")

p6 <-ggplot(qcobject@meta.data, aes(x = orig.ident, y = TSS.enrichment)) +
  geom_violin(aes(fill = orig.ident), show.legend = FALSE) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_minimal() +
  labs(title = "TSS enrichment",
       x = "Sample",
       y = "TSS enrichment score")

# Combine plots
combined_plot <- (p1 + p2 + p3) / (p4 + p5 + p6) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

# distribution plots
p <- ggplot(qcobject@meta.data, aes(x = .data[[metric]])) +
  geom_histogram(bins = 100, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = paste(metric, "distribution"),
       x = metric,
       y = "Count")

p <- p + geom_vline(xintercept = threshold_low, color = "red", linetype = "dashed")
p <- p + geom_vline(xintercept = threshold_high, color = "red", linetype = "dashed")