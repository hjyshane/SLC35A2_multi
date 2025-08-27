base_dir <- "~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
setwd(base_dir)

list <- list.files("./DEG/sig_csv/", full.names = T)
full <- data.frame()

for (file in list) {
    updata <- read_csv(file) %>%
        filter(Regulation == "Up") %>%
        arrange(p_val_adj, avg_log2FC)

    dndata <- read_csv(file) %>%
        filter(Regulation == "Down") %>%
        arrange(p_val_adj, avg_log2FC)

    data <- rbind(updata, dndata)
    full <- rbind(full, data)

    data <- data.frame()
}

write.csv(full, file = "./DEG/sig_csv/top10_up_down.csv")
write.csv(full, file = "./DEG/sig_csv/ful_deg.csv")
full_deg <- read_csv("./DEG/sig_csv/ful_deg.csv")
full_deg[1] <- NULL
colnames(full_deg)[1] <- "Gene"


list <- list.files("./GO/csv/", full.names = T)
full <- data.frame()

for (file in list) {

    file_name <- basename(file)
    parts <- strsplit(gsub(".csv", "", file_name), "_")[[1]]
    cell_type <- parts[1]
    comparison_name <- paste(parts[-1], collapse = "_")

    data <- read_csv(file) %>%
        mutate(
            Cell_type = cell_type,
            Comparison = comparison_name
        ) %>%
        arrange(p.adjust, Fold_enrichment) %>%
        slice(1:10)
    full <- rbind(full, data)
    }

write.csv(full, file = "./GO/csv/top10go.csv")


full_deg <- read_csv("./DEG/sig_csv/ful_deg.csv")
full_deg[1] <- NULL
colnames(full_deg)[1] <- "Gene"


cell_type_order <- c(
    "Astrocytes", "CGE", "ExcitatoryNeuronsCA1",
    "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsImmatureDG",
    "ExcitatoryNeuronsMatureDG", "L23ExcitatoryNeurons",
    "L45ExcitatoryNeurons", "L56ExcitatoryNeurons", "MGE",
    "Oligodendrocytes", "Microglia"
)

# Define desired comparison (adjust as needed)
target_comparison <- "PTZvsSAL_1hr"

# STEP 1: Filter DEG data for one comparison and p-adj threshold
ptz_long <- full_deg %>%
    filter(comparison == target_comparison, p_val_adj < 0.05) %>%
    distinct(gene, cell_type)

# STEP 2: Define and apply desired cell type order
cell_type_order <- c(
    "Astrocytes", "CGE", "ExcitatoryNeuronsCA1",
    "ExcitatoryNeuronsCA3", "ExcitatoryNeuronsImmatureDG",
    "ExcitatoryNeuronsMatureDG", "L23ExcitatoryNeurons",
    "L45ExcitatoryNeurons", "L56ExcitatoryNeurons",
    "MGE", "Oligodendrocytes", "Microglia"
)

ptz_long <- ptz_long %>%
    mutate(cell_type = factor(cell_type, levels = cell_type_order))

# STEP 3: Determine gene order by left-most (earliest) cell_type appearance
gene_order <- ptz_long %>%
    mutate(cell_rank = as.numeric(cell_type)) %>%
    group_by(gene) %>%
    summarise(first_cell = min(cell_rank)) %>%
    arrange(first_cell) %>%
    pull(gene)

# STEP 4: Apply gene factor level ordering
ptz_long <- ptz_long %>%
    mutate(gene = factor(gene, levels = gene_order))

ggplot(ptz_long, aes(x = cell_type, y = gene)) +
    geom_tile(fill = "black", width = 0.9, height = 0.2) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank()
    ) +
    labs(
        title = "DEG span across cell types (PTZvsSAL_1hr)",
        x = "Cell Type"
    ) +
    coord_fixed(ratio = 0.25)
