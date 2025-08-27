# read this function with
# source('path/to/this/file/16_find_dg.R')

#' Run differential expression analysis for specific cell types and group comparisons
#'
#' This function performs DEG analysis using either Wilcox or MAST test,
#' generates volcano plots, and saves background and significant results per comparison.
#'
#' @param seurat_obj Seurat object containing the data.
#' @param cell_type_meta Column name in `meta.data` that holds cell type labels (string).
#' @param comparison_meta Column name in `meta.data` that defines grouping for comparison (string).
#' @param cell_type Cell type to subset and analyze.
#' @param comparisons A list of named comparisons (e.g., list(PTZvsSAL = list(name = "PTZvsSAL", group1 = "PTZ", group2 = "SAL"))).
#' @param minimum_cnt Minimum number of cells required per group (default = 3).
#' @param mode Test method: "wilcox", "MAST", or "LR"(for ATAC) (default = "wilcox").
#' @param p_value Adjusted p-value threshold for significance (default = 0.05).
#' @param fc_value Minimum absolute log2 fold change to call a gene significant (default = 0.2).
#' @param save Logical; whether to save results (default = TRUE).
#' @param bg_dir Directory to save full background DEG tables.
#' @param sig_dir Directory to save filtered significant DEG tables.
#' @param plot_dir Directory to save volcano plots.
#'
#' @return None. Results are optionally saved to CSVs and PNGs.
#' @export

run_dg <- function(
    seurat_obj,
    cell_type_meta,
    comparison_meta,
    cell_types,
    comparisons,
    minimum_cnt = 10,
    mode = "MAST", 
    p_value = 0.05,
    fc_value = 0.2,
    save = TRUE,
    bg_dir = NULL,
    sig_dir = NULL,
    plot_dir = NULL) {
  # Save check
  if (save) {if (is.null(bg_dir)) {stop("You must provide 'bg_dir' when save = TRUE.")}
    if (!dir.exists(bg_dir)) {dir.create(bg_dir, recursive = TRUE)}
  }
  if (save) {if (is.null(sig_dir)) {stop("You must provide 'sig_dir' when save = TRUE.")}
    if (!dir.exists(sig_dir)) {dir.create(sig_dir, recursive = TRUE)}
  }
  if (save) {if (is.null(plot_dir)) {stop("You must provide 'plot_dir' when save = TRUE.")}
    if (!dir.exists(plot_dir)) {dir.create(plot_dir, recursive = TRUE)}
  }
  for (ct in cell_types) {
    # For each cell type, loop through each comparison group
    for (comparison in comparisons) {
      # Subset for cell type
      
      # Set idents
      Seurat::Idents(seurat_obj) <- cell_type_meta
      
      # subsetting cell type in each condition
      cell_type_data <- subset(seurat_obj,
                               subset = cell_type == ct & Sample %in% c(comparison$group1, comparison$group2)) # cell_type is metadata column name for cell type, Sample is for group
      
      if (length(cell_type_meta >= minimum_cnt)) {
        # Set identities
        Seurat::DefaultAssay(cell_type_data) <- "RNA"
        Seurat::Idents(cell_type_data) <- "Sample"
        
        # 표현행렬(로그 정규화; MAST는 log expr 선호)과 메타
        emat <- GetAssayData(cell_type_data, slot = "data", assay = "RNA")   # genes x cells
        cdata <- cell_type_data@meta.data
        
        # CDR(검출 feature 수) 계산: nFeature_RNA 사용
        cdata$CDR <- cell_type_data$nFeature_RNA
        
        # 그룹을 이항 팩터로 정리
        cdata$group <- factor(cdata$Condition, levels = c(group2, group1))  # ref=group2
        
        # MAST SingleCellAssay 구성
        sca <- MAST::FromMatrix(exprsArray = as.matrix(emat),
                                cData      = cdata)
        
        # zlm 피팅: group + CDR (공변량 1개만)
        zf <- zlm(~ group + CDR, sca, ebayes = TRUE, method = "bayesglm")
        
        # 성분별 LRT 수행: group 계수에 대해
        # summary()는 component 열에 C(연속), D(이산), H(허들/합성)를 포함
        summary_z <- summary(zf, doLRT = "group" ~ 1)  # 전체 group 효과 검정
        res <- summary_z$datatable %>%
          as.data.frame()
        
        # component별 p값 추출
        # C: 연속, D: 이산, H: 허들(합성)
        p_C <- res %>% dplyr::filter(component == "C") %>%
          dplyr::select(primerid, `Pr(>Chisq)`) %>%
          rename(p_cont = `Pr(>Chisq)`)
        p_D <- res %>% dplyr::filter(component == "D") %>%
          dplyr::select(primerid, `Pr(>Chisq)`) %>%
          rename(p_disc = `Pr(>Chisq)`)
        p_H <- res %>% dplyr::filter(component == "H") %>%
          dplyr::select(primerid, `Pr(>Chisq)`) %>%
          rename(p_hurdle = `Pr(>Chisq)`)
        
        # 연속 성분의 계수(=효과크기 유사; log-scale) 가져오기
        # coef()에서 "groupgroup1" 항을 continuous(Beta)에서 추출
        coefs <- coef(zf)
        betaC <- coefs$C  # 연속 성분 계수 행렬 (genes x coefficients)
        if(!"groupgroup1" %in% colnames(betaC)){
          stop("Coefficient 'groupgroup1' not found; check group coding.")
        }
        logFC_C <- betaC[, "groupgroup1", drop=FALSE] %>%
          as.data.frame() %>%
          tibble::rownames_to_column("gene") %>%
          rename(logFC_cont = "groupgroup1")
        
        # 검출률 차이(이산 성분) 계수도 원하면 추출 가능(해석은 log-odds)
        betaD <- coefs$D
        logOR_D <- NULL
        if("groupgroup1" %in% colnames(betaD)){
          logOR_D <- betaD[, "groupgroup1", drop=FALSE] %>%
            as.data.frame() %>%
            tibble::rownames_to_column("gene") %>%
            rename(logOR_disc = "groupgroup1")
        }
        
        # 세 성분 p값 병합
        p_all <- p_C %>%
          left_join(p_D, by = c("primerid")) %>%
          left_join(p_H, by = c("primerid")) %>%
          rename(gene = primerid)
        
        # Seurat pct.1/2와 avg_log2FC도 참고용으로 붙임
        seurat_like <- FindMarkers(
          object          = cell_type_data,
          ident.1         = comparison$group1,
          ident.2         = comparison$group2,
          assay           = "RNA",
          test.use        = "MAST",
          logfc.threshold = 0,
          min.pct         = 0.25,
          min.cells.group = 10,
          latent.vars     = "nFeature_RNA"
        ) %>%
          tibble::rownames_to_column("gene") %>%
          dplyr::select(gene, avg_log2FC, pct.1, pct.2)
        
        bg_results <- p_all %>%
          left_join(logFC_C, by = "gene") %>%
          left_join(seurat_like, by = "gene")
        
        # Add metadata info
        bg_results$gene <- rownames(bg_results)
        bg_results$cell_type <- ct
        bg_results$comparison <- comparison$name
        
        # Add regulation info
        bg_results <- bg_results %>%
          mutate(Regulation = case_when(
            avg_log2FC > 0 ~ "Up",
            avg_log2FC < 0 ~ "Down",
            TRUE ~ "NoChange"))
        
        # Get significant genes - DEGs
        sig_results <- bg_results %>%
          filter(p_val_adj < p_value & abs(avg_log2FC) >= fc_value)
        
        # Create VolcanoPlot
        if (nrow(bg_results) == 0) {
          message("There are no genes/peaks to plot")
          next
        } else {

          # A. 연속 성분: p_cont vs logFC_cont
          dfC <- bg_results %>% dplyr::filter(!is.na(p_cont), !is.na(logFC_cont))
          plt_cont <- EnhancedVolcano::EnhancedVolcano(
            dfC,
            lab = dfC$gene,
            x = 'logFC_cont',
            y = 'p_cont',
            title = paste0('[C] Continuous: ', comparison$group1, ' vs ', comparison$group2),
            xlab = 'Continuous logFC (MAST beta, log-scale)',
            ylab = bquote(~-Log[10]~italic(P[C])),
            pCutoff = 0.05,
            FCcutoff = 0.25,        # 연속 성분용 효과크기 컷(상황에 맞게)
            pointSize = 1.0,
            labSize = 3.0,
            drawConnectors = TRUE,
            max.overlaps = 20
          )
          
          # B. 이산 성분: p_disc vs logOR_disc(로그 오즈; 검출률 차이의 효과크기)
          dfD <- bg_results %>% dplyr::filter(!is.na(p_disc), !is.na(logOR_disc))
          plt_disc <- EnhancedVolcano::EnhancedVolcano(
            dfD,
            lab = dfD$gene,
            x = 'logOR_disc',       # 검출률 차이의 효과크기(로그 오즈)
            y = 'p_disc',
            title = paste0('[D] Discrete: ', comparison$group1, ' vs ', comparison$group2),
            xlab = 'Discrete log-odds (detection probability)',
            ylab = bquote(~-Log[10]~italic(P[D])),
            pCutoff = 0.05,
            FCcutoff = 0.5,         # 예: |logOR| >= 0.5 (데이터에 맞게 조정)
            pointSize = 1.0,
            labSize = 3.0,
            drawConnectors = TRUE,
            max.overlaps = 20
          )
          
          # C. 합성 성분: p_hurdle vs (참고) Seurat avg_log2FC
          dfH <- bg_results %>% dplyr::filter(!is.na(p_hurdle), !is.na(avg_log2FC))
          plt_hurdle <- EnhancedVolcano::EnhancedVolcano(
            dfH,
            lab = dfH$gene,
            x = 'avg_log2FC',
            y = 'p_hurdle',
            title = paste0('[H] Hurdle (combined): ', comparison$group1, ' vs ', comparison$group2),
            subtitle = 'x=avg_log2FC(continuous), y=combined p',
            xlab = bquote(~Log[2]~ 'FC (Seurat summary)'),
            ylab = bquote(~-Log[10]~italic(P[H])),
            pCutoff = 0.05,
            FCcutoff = 0.25,
            pointSize = 1.0,
            labSize = 3.0,
            drawConnectors = TRUE,
            max.overlaps = 20
          )
          
        }
        # Save
        if (save) {
          if (nrow(bg_results) == 0) {
            message("There are no genes/peaks to found")
            next
          } else {
            write.csv(bg_results, file = file.path(bg_dir, paste0(ct, "_", comparison$name, ".csv")))
            ggsave(plot = plt_cont, file = file.path(plot_dir, paste0("Countinuous_", ct, "_", comparison$name, ".png")), width = 10, height = 12, dpi = 300)
            ggsave(plot = plt_disc, file = file.path(plot_dir, paste0("Discrete_", ct, "_", comparison$name, ".png")), width = 10, height = 12, dpi = 300)
            ggsave(plot = plt_hurdle, file = file.path(plot_dir, paste0("Hurdle_", ct, "_", comparison$name, ".png")), width = 10, height = 12, dpi = 300)
            
          }
          if (nrow(sig_results) == 0) {
            message("There are no deg genes/peaks to found")
            next
          } else {
            write.csv(sig_results, file = file.path(sig_dir, paste0(ct, "_", comparison$name, ".csv")))
          }
          
        }
      } else next
    }
  }
}



