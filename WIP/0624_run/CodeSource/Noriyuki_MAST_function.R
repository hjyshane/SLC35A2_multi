test <- FindMarkersMAST(combined_object, ident.1 = "PTZ_1hr", ident.2 = "SAL_1hr", group.by = "Sample", subset.ident = "ExcitatoryNeuronsCA1")

FindMarkersMAST <-
  function(
    object,
    ident.1,
    ident.2,
    group.by,
    subset.ident,
    min.pct = 0.1,
    random.effect = character(0)
  ) {
    # subset cluster of interest
    object.subset <- subset(object, idents = subset.ident)
    
    # Set default assay to RNA and normalize
    DefaultAssay(object.subset) <- "RNA"
    object.subset <- NormalizeData(object.subset)
    
    # Extract counts from subset object and filter out genes expressed in less than 10% of cells
    counts <- as.matrix(object.subset@assays$RNA$data)
    counts <- counts[rowSums(counts > 0) / (ncol(counts)) >= min.pct,]
    
    # Prepare a list of genes to model
    features <- data.frame(primerid = rownames(counts))
    
    # Include important metadata
    levels <- c(ident.2, ident.1)
    metadata <- data.frame(
      condition = factor(
        object.subset@meta.data[[group.by]],
        levels = levels)
    )
    
    # If a random effect is specified, make sure to include it in the metadata
    if (length(random.effect)) {
      metadata$re <- factor(object.subset@meta.data[[random.effect]])
    }
    
    # Create single cell assay
    assay <- MAST::FromMatrix(
      exprsArray = counts,
      cData = metadata,
      fData = features
    )
    
    # Recalculate cellular detection rate and add to assay
    detections <- colSums(SummarizedExperiment::assay(assay) > 0)
    SummarizedExperiment::colData(assay)$cngeneson <- scale(detections)
    
    # Perform zero-inflated regression
    # If a random effect is specified, the glmer method is used instead of the default method
    if (length(random.effect)) {
      fit <-
        MAST::zlm(
          ~ condition + cngeneson + (1 | re),
          assay,
          method = 'glmer',
          ebayes = FALSE,
          strictConvergence = FALSE
        )
    } else {
      fit <- MAST::zlm(formula = ~ condition + cngeneson, sca = assay)
    }
    
    # Summarize experiment in data table
    experiment <- paste0('condition', levels[2])
    summary <- MAST::summary(fit, doLRT = experiment)
    table <- summary$datatable
    
    # convert p value to FDR
    hurdle <- merge(
      table[
        contrast == experiment & component == "H",
        .(primerid, `Pr(>Chisq)`)
      ], table[
        contrast == experiment & component == "logFC",
        .(primerid, coef)
      ], by = "primerid"
    )
    hurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
    hurdle <- stats::na.omit(as.data.frame(hurdle))
    
    return(hurdle)
  }
