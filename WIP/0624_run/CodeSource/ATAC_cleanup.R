library(Signac)         # scATAC-seq analysis
library(Seurat)         # scRNA-seq analysis
library(tidyverse)   # tidyr, dplyr, readr, purr, stringr, ggplot2, tibble, forcats, ludridate
library(ggplot2)          # Core visualization package
library(patchwork)        # Arrange multiple ggplots
library(AnnotationDbi)    # Database interface
library(BSgenome.Mmusculus.UCSC.mm10)  # Mouse genome
library(rtracklayer)                   # Import/export genomic data
library(JASPAR2020)                    # Transcription factor motifs
library(TFBSTools)                      # Analyze transcription factor binding sites
combined_object <- qs::qread('~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/filtered_50_sctall.qs')
mm10 <-import("/igm/apps/10X_chromium/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz")


DefaultAssay(combined_object) <- "ATAC"
seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- seqlevelsStyle(Annotation(combined_object))

# Check  sequence levels contain meaningful peaks or annotations
seqinfo(Annotation(combined_object))
seqlevels(Annotation(combined_object))
table(seqnames(Annotation(combined_object)))


# Check for mismatches - most likely non-standard chr will pop up.
setdiff(seqlevels(Annotation(combined_object)), seqlevels(BSgenome.Mmusculus.UCSC.mm10))

# keep only standard chromosome, Filter Non-Matching Seqlevels
Annotation(combined_object) <- keepStandardChromosomes(Annotation(combined_object), pruning.mode = "coarse")

# Update the seqlengths in your annotations
seq_lengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
seqlengths(Annotation(combined_object)) <- seq_lengths[names(seqlengths(Annotation(combined_object)))]

# Check for mismatches
setdiff(seqlevels(Annotation(combined_object)), seqlevels(BSgenome.Mmusculus.UCSC.mm10))

# Check  sequence levels contain meaningful peaks or annotations
seqinfo(Annotation(combined_object))
seqlevels(Annotation(combined_object))
table(seqnames(Annotation(combined_object)))

# since I imported Grange object, and it is missing row names (it is not inherent)
# Add identifiers as metadata
mcols(Annotation(combined_object))$gene_id <- make.unique(as.character(mcols(Annotation(combined_object))$gene_id))
names(Annotation(combined_object)) <- mcols(Annotation(combined_object))$gene_id

# Verify the result
head(mcols(Annotation(combined_object))$gene_id)
head(names(Annotation(combined_object)))

# this is extra if needed a separate obj for ID.
# # Extract gene IDs as a vector
# gene_ids <- make.unique(as.character(mcols(Annotation(combined_object))$gene_id))
# # Convert GRanges to data.frame with row names
# annotation_df <- as.data.frame(Annotation(combined_object))
# rownames(annotation_df) <- gene_ids
# # Use gene_ids for other operations
# head(gene_ids)

# check GC content
# Check GC content summary
summary(combined_object@assays$ATAC@meta.features$GC.percent)

# calculate GC content, if missing.
combined_object <- RegionStats(
    object = combined_object,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    assay = "ATAC",
    verbose = T
)

# work on GC content
# Identify regions with missing GC content
missing_gc <- rownames(combined_object@assays$ATAC@meta.features)[
    is.na(combined_object@assays$ATAC@meta.features$GC.percent)
]

# Check seqlevels of regions with missing GC content
missing_seqnames <- sapply(strsplit(missing_gc, "-"), `[`, 1)
table(missing_seqnames)

# to fix. NA GC are in non-standard chr which I thought I cleaned.
# my object filtered out non-standard chromosomes. but  it is still showing up (it is the couase of NA value..)
# Extract seqnames from rownames in meta.features
all_seqnames <- sapply(strsplit(rownames(combined_object@assays$ATAC@meta.features), "-"), `[`, 1)

# Identify regions on standard chromosomes
valid_regions <- which(all_seqnames %in% seqlevels(Annotation(combined_object)))

# Subset the ATAC assay by valid regions to get ride of non-standard data (which have NA GC content value)
combined_object@assays$ATAC <- subset(
    combined_object@assays$ATAC,
    features = rownames(combined_object@assays$ATAC@meta.features)[valid_regions]
)

# now we should have no NA value.
summary(combined_object@assays$ATAC@meta.features$GC.percent)

pwm <- TFBSTools::getMatrixSet(
    x = JASPAR2020::JASPAR2020,
    opts = list(collection = "CORE", tax_group = "vertebrates", species = 10090) # 10090 is the NCBI Taxonomy ID for mouse
)

combined_object <- AddMotifs(
    object = combined_object,
    genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,  # Use the filtered genome
    pfm = pwm
)

# Check object
combined_object@assays$ATAC@motifs
summary(combined_object@assays$ATAC@meta.features$GC.percent)

# Get unique column names without ".1" duplicates
meta.feature <- combined_object@assays$ATAC@meta.features

# Remove columns that are duplicates (keeping the first occurrence)
meta.feature <- meta.feature[, !duplicated(colnames(meta.feature))]
meta.feature <- meta.feature[, !grepl("\\.1$", colnames(meta.feature))]

# Assign the cleaned metadata back
combined_object@assays$ATAC@meta.features <- meta.feature

# final check
summary(combined_object@assays$ATAC@meta.features$GC.percent)
head(combined_object[["ATAC"]]@motifs)
dim(combined_object[['ATAC']][[]])
head(combined_object[['ATAC']][[]])
table(seqnames(Annotation(combined_object)))
seqlevels(Annotation(combined_object))
seqinfo(Annotation(combined_object))

qs::qsave(combined_object, file = "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/qsave/filtered_50_cells_atac_cleaned.qs")