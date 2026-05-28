# =============================================================================
# GO Enrichment Analysis for 4 Stoats DEG Files
# Method: clusterProfiler ORA + GSEA
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# --- Files to analyze --------------------------------------------------------
files <- c(
  "260510_Stoats_DEG_renamed (1).csv",
  "260513_Stoats_IAvIB_DEG_Fixed.csv",
  "260513_Stoats_NAvIA_DEG_Fixed.csv",
  "260513_Stoats_NBvIB_DEG_Fixed.csv"
)

names(files) <- c(
  "Original_DEG",
  "IAvIB",
  "NAvIA",
  "NBvIB"
)

# --- Thresholds --------------------------------------------------------------
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0

# --- Function to map genes ---------------------------------------------------
map_to_entrez <- function(gene_symbols, orgdb = org.Hs.eg.db) {
  bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = orgdb
  )
}

# --- Function to run GO ORA --------------------------------------------------
run_GO_ORA <- function(gene_ids, universe_ids, ontology = "BP") {
  enrichGO(
    gene          = gene_ids,
    universe      = universe_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

# --- Function to save plots safely ------------------------------------------
save_plot <- function(plot_obj, filename, width = 10, height = 8) {
  ggsave(filename, plot = plot_obj, width = width, height = height, dpi = 150)
  cat("Saved:", filename, "\n")
}

# --- Main analysis function --------------------------------------------------
run_GO_analysis <- function(file, comparison_name) {
  
  cat("\n============================================================\n")
  cat("Running GO analysis for:", comparison_name, "\n")
  cat("File:", file, "\n")
  cat("============================================================\n")
  
  # Create output folder
  outdir <- paste0("GO_", comparison_name)
  dir.create(outdir, showWarnings = FALSE)
  
  # Load DEG file
  deg <- read.csv(file, row.names = 1)
  
  cat("Total genes:", nrow(deg), "\n")
  cat("Columns:", paste(colnames(deg), collapse = ", "), "\n")
  
  # Significant DEGs
  sig <- deg[
    !is.na(deg$padj) &
      deg$padj < PADJ_CUTOFF &
      abs(deg$log2FoldChange) > LFC_CUTOFF,
  ]
  
  sig_up <- sig[sig$log2FoldChange > 0, ]
  sig_down <- sig[sig$log2FoldChange < 0, ]
  
  cat("Significant DEGs:", nrow(sig), "\n")
  cat("  Up-regulated:", nrow(sig_up), "\n")
  cat("  Down-regulated:", nrow(sig_down), "\n")
  
  # Map genes
  universe_map <- map_to_entrez(rownames(deg))
  universe_ids <- universe_map$ENTREZID
  
  sig_map <- map_to_entrez(rownames(sig))
  sig_up_map <- map_to_entrez(rownames(sig_up))
  sig_down_map <- map_to_entrez(rownames(sig_down))
  
  # GO ORA
  go_BP_all <- run_GO_ORA(sig_map$ENTREZID, universe_ids, "BP")
  go_MF_all <- run_GO_ORA(sig_map$ENTREZID, universe_ids, "MF")
  go_CC_all <- run_GO_ORA(sig_map$ENTREZID, universe_ids, "CC")
  
  go_BP_up <- run_GO_ORA(sig_up_map$ENTREZID, universe_ids, "BP")
  go_BP_down <- run_GO_ORA(sig_down_map$ENTREZID, universe_ids, "BP")
  
  # GSEA
  gene_list <- deg$log2FoldChange
  names(gene_list) <- rownames(deg)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_list_entrez <- gene_list[names(gene_list) %in% universe_map$SYMBOL]
  names(gene_list_entrez) <- universe_map$ENTREZID[
    match(names(gene_list_entrez), universe_map$SYMBOL)
  ]
  
  gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
  gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
  
  set.seed(42)
  
  gsea_BP <- gseGO(
    geneList      = gene_list_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  # Save CSV results
  write.csv(as.data.frame(go_BP_all),
            file.path(outdir, paste0(comparison_name, "_GO_ORA_BP_allDEGs.csv")),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_MF_all),
            file.path(outdir, paste0(comparison_name, "_GO_ORA_MF_allDEGs.csv")),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_CC_all),
            file.path(outdir, paste0(comparison_name, "_GO_ORA_CC_allDEGs.csv")),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_BP_up),
            file.path(outdir, paste0(comparison_name, "_GO_ORA_BP_upDEGs.csv")),
            row.names = FALSE)
  
  write.csv(as.data.frame(go_BP_down),
            file.path(outdir, paste0(comparison_name, "_GO_ORA_BP_downDEGs.csv")),
            row.names = FALSE)
  
  write.csv(as.data.frame(gsea_BP),
            file.path(outdir, paste0(comparison_name, "_GSEA_BP_results.csv")),
            row.names = FALSE)
  
  # Save plots
  if (nrow(as.data.frame(go_BP_all)) > 0) {
    p_dot <- dotplot(
      go_BP_all,
      showCategory = 20,
      title = paste0("GO Biological Process — ", comparison_name)
    ) +
      theme_bw(base_size = 11)
    
    save_plot(
      p_dot,
      file.path(outdir, paste0(comparison_name, "_GO_dotplot_BP_allDEGs.png")),
      width = 10,
      height = 9
    )
  }
  
  if (nrow(as.data.frame(go_BP_up)) > 0 &&
      nrow(as.data.frame(go_BP_down)) > 0) {
    
    df_up <- as.data.frame(go_BP_up)[
      1:min(10, nrow(as.data.frame(go_BP_up))),
      c("Description", "p.adjust", "Count")
    ]
    
    df_down <- as.data.frame(go_BP_down)[
      1:min(10, nrow(as.data.frame(go_BP_down))),
      c("Description", "p.adjust", "Count")
    ]
    
    df_up$direction <- "Up-regulated"
    df_down$direction <- "Down-regulated"
    
    df_combined <- rbind(df_up, df_down)
    df_combined$Description <- factor(
      df_combined$Description,
      levels = rev(unique(df_combined$Description))
    )
    
    p_bar <- ggplot(
      df_combined,
      aes(x = Count, y = Description, fill = direction)
    ) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c(
        "Up-regulated" = "#E63946",
        "Down-regulated" = "#457B9D"
      )) +
      labs(
        title = paste0("Top GO-BP Terms by Direction — ", comparison_name),
        x = "Gene Count",
        y = NULL,
        fill = ""
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "top")
    
    save_plot(
      p_bar,
      file.path(outdir, paste0(comparison_name, "_GO_barplot_UP_vs_DOWN.png")),
      width = 11,
      height = 7
    )
  }
  
  if (nrow(as.data.frame(go_BP_all)) > 1) {
    go_BP_sim <- pairwise_termsim(go_BP_all)
    
    p_emap <- emapplot(go_BP_sim, showCategory = 30) +
      ggtitle(paste0("GO-BP Enrichment Map — ", comparison_name))
    
    save_plot(
      p_emap,
      file.path(outdir, paste0(comparison_name, "_GO_emapplot_BP.png")),
      width = 12,
      height = 10
    )
  }
  
  if (nrow(as.data.frame(gsea_BP)) > 0) {
    p_ridge <- ridgeplot(gsea_BP, showCategory = 20) +
      labs(title = paste0("GSEA GO-BP Ridge Plot — ", comparison_name)) +
      theme_bw(base_size = 10)
    
    save_plot(
      p_ridge,
      file.path(outdir, paste0(comparison_name, "_GSEA_ridgeplot_BP.png")),
      width = 11,
      height = 9
    )
  }
  
  if (nrow(as.data.frame(gsea_BP)) > 0) {
    p_gsea_dot <- dotplot(
      gsea_BP,
      showCategory = 20,
      split = ".sign",
      title = paste0("GSEA GO-BP — ", comparison_name)
    ) +
      facet_grid(. ~ .sign) +
      theme_bw(base_size = 10)
    
    save_plot(
      p_gsea_dot,
      file.path(outdir, paste0(comparison_name, "_GSEA_dotplot_BP.png")),
      width = 13,
      height = 9
    )
  }
  
  cat("\nFinished:", comparison_name, "\n")
}

# --- Run GO analysis for all four files --------------------------------------
for (comparison_name in names(files)) {
  run_GO_analysis(
    file = files[comparison_name],
    comparison_name = comparison_name
  )
}

cat("\nDone! GO analysis completed for all four files.\n")