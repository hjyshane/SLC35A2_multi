    group_by(cell_type, Sample) %>%
    count() %>%
    as.data.frame()
colnames(cellcount) <- c("Cell_Type", "Group", "Cell_Count")
cell_countPlot <- ggplot(cellcount) +
    aes(x = Group, y = Cell_Count, fill = Cell_Type) +
    # stat identity will use fill = cell type as identity and position fill will make proportional bar. if fill stack is used, it will be just stacked.
    geom_bar(stat = "identity", position = "fill") +
    labs(x = "",
         y = "Proportion (%)",
         title = "Cell type proportions by treatment group",
         # legend title can be accessed by fill/color whatever you used to plot them.
         fill = "Cell Type") +
    theme(panel.background = element_blank()) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal()

ggsave(cell_countPlot, file = "./Plots/cellount.png", width = 8, height = 8, dpi = 300)
