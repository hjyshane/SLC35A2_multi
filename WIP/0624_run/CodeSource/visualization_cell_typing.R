# Visualization for more cell typing using VlnPlot, FeaturePlot, Heatmap and DotPlot

# DimPlot 
make_dimplot <- function(
    object,
    reduction,
    groupby = NULL,
    save = TRUE,
    save_dir = NULL
    ) {
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  
  # 
  
  }