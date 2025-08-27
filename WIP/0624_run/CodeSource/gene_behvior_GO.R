# get file list
output_dir_base <-"~/PTZ_ATAC_scRNA_072024/WIP/0624_run/"
sig_gene_files <- list.files(path = file.path(output_dir_base, "DEG", "sig_csv"), pattern = ".csv$", full.names = TRUE)

# combine list
deg_list <- lapply(
    sig_gene_files,
    function(file) {
        read_csv(file) %>%
            mutate(condition = paste(cell_type, comparison, sep = "_"))
    }
)

deg_all <- bind_rows(deg_list)

deg_wide <- deg_all %>%
    # filter(comparison %in% c("PTZvsSAL_1hr", "PTZvsSAL_24hr", "24hrvs1hr_PTZ")) %>%
    dplyr::select(gene, cell_type, comparison, avg_log2FC) %>%
    pivot_wider(
        names_from = comparison,
        values_from = avg_log2FC
    )

deg_wide <- deg_wide %>%
    mutate(
        pattern = case_when(
            PTZvsSAL_1hr > 0 & PTZvsSAL_24hr > 0 ~ "Increase maintain",
            PTZvsSAL_1hr > 0 & PTZvsSAL_24hr < 0 ~ "Increase then Reversal",
            PTZvsSAL_1hr > 0 & is.na(PTZvsSAL_24hr) ~ "Increase then Baseline",
            PTZvsSAL_1hr < 0 & PTZvsSAL_24hr < 0 ~ "Sustained Decrease",
            PTZvsSAL_1hr < 0 & PTZvsSAL_24hr > 0 ~ "Decrease then Reversal",
            PTZvsSAL_1hr < 0 & is.na(PTZvsSAL_24hr) ~ "Decrease then Baseline",
            is.na(PTZvsSAL_1hr) & PTZvsSAL_24hr > 0 ~ "Late Upregulation",
            is.na(PTZvsSAL_1hr) & PTZvsSAL_24hr < 0 ~ "Late Downregulation",
            is.na(PTZvsSAL_1hr) & is.na(PTZvsSAL_24hr) & is.na(`24hrvs1hr_PTZ`) ~ "No change in PTZvsSAL",
            is.na(PTZvsSAL_1hr) & is.na(PTZvsSAL_24hr) & `24hrvs1hr_PTZ` < 0 ~ "Decrease in PTZ",
            is.na(PTZvsSAL_1hr) & is.na(PTZvsSAL_24hr) & `24hrvs1hr_PTZ` > 0 ~ "Increase in PTZ",
            TRUE ~ "Other")
    )
write.csv(deg_wide, "./DEG/gene_behavior.csv")


deg_fileterd <- deg_wide %>%
    filter(deg_wide$pattern %in% c("Late Downregulation", "Sustained Decrease" , "Decrease then Reversal",
                                   "Increase then Reversal", "Increase maintain" , "Late Upregulation"))

deg_filtered_gene <- deg_fileterd$gene
bg_gene <- full$gene

go_result <- clusterProfiler::enrichGO(
    gene = deg_filtered_gene,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    universe = bg_gene,
    pAdjustMethod = "bonferroni",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05)
