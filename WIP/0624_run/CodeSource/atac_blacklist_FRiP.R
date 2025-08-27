# Blacklist Ratio and FRiP 

cobj <- qs::qread("./qsave/combined_cell_type.qs")
DefaultAssay(cobj) <- "RNA"
cobj <- NormalizeData(cobj)

DefaultAssay(cobj) <- "ATAC"

# get FRiP
# The percentage of total ATAC-seq fragments (per cell) that fall within called peak regions.
# High pct_reads_in_peaks suggests that the majority of the ATAC signal is concentrated in biologically relevant open chromatin (peaks), not random noise or background signal.
# FRiP = (reads in peaks) / (total reads/fragments)
cobj$pct_reads_in_peaks <- cobj$atac_peak_region_fragments / cobj$atac_fragments * 100

# add blacklist ratio
cobj$blacklist_fraction <- FractionCountsInRegion(cobj, assay = "ATAC", regions = blacklist_mm10)

