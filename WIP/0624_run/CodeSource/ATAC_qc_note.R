# 1. pct_reads_in_peaks (a.k.a. FRiP)
# 	Definition: % of total fragments that fall within called peak regions.
# 	Good values:
#   > 15% is minimum acceptable
# 	> 25–30% is high quality
# 	Use: Filter low-quality cells; reflects signal-to-noise ratio.
# 	Low values suggest** background noise or dying cells.
# 
# 2. nCount_ATAC
# 	Definition: Total number of fragments per cell.
# 	Good values (typical range in 10X multiome):
#   	> 3,000 – 10,000 is minimum
# 	> 20,000 – 50,000+ is strong, but watch for doublets
# 	Use:
#   	Filter low coverage cells (<3k–5k)
# 	Check for extreme outliers (potential multiplets)

# 3. nFeature_ATAC
# 	Definition: Number of peaks detected (nonzero) per cell.
# 	Good values:
#   	>1,000 – 2,000 peaks is baseline
# 	>5,000 = high quality, but >15,000 could be a doublet
# 	Use: Gauge chromatin accessibility breadth

# 4. TSS.enrichment
# 
# From TSSEnrichment() in Signac
# 	Definition: Enrichment of reads near transcription start sites (TSSs).
# 	Good values:
#   	> 2 = minimal acceptable
# 	> 4 = moderate
# 	> 6–8+ = high-quality cells
# 	Why it matters: High-quality nuclei show sharp peaks of accessibility near TSSs.
# 
# 
# 5. NucleosomeSignal (sometimes written nucleosome_signal)
# 	Definition: Ratio of mono- to di-nucleosome fragment counts.
# 	Good values:
#   	< 4 is typical
# 	< 2 is high-quality
# 	High values (>10): indicate degraded chromatin (bad cell integrity)
# 
# 
# 6. blacklist_fraction
# 	Definition: % of fragments falling into genomic blacklist regions (like centromeres, telomeres, repeat-heavy regions).
# 	Good values:
#   	< 0.05 (i.e. <5%) is preferred
# 	Why: High values suggest non-specific or artifactual signal