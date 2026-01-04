# ============================================================
# PATHWAY ENRICHMENT / GSEA desde DE (CSV) — usando msigdbr + fgsea
# - Input: summary_tables_final/DE_pseudobulk_{ct}_Healthy_vs_Cirrhosis.csv
# - Output: summary_tables_final/GSEA_* + figures_final/GSEA_*
# - Extra: summary_tables_final/GSEA_MSigDB_summary_topPathways_Level1.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

PROJECT_ROOT <- getwd()
SUM_DIR <- file.path(PROJECT_ROOT, "summary_tables_final")
FIG_DIR <- file.path(PROJECT_ROOT, "figures_final")
dir.create(SUM_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Level1 celltypes
celltypes <- c("T","Mono","NK","B","DC")

# ------------------------------------------------------------
# Helpers: msigdbr "compatible" (viejo/nuevo)
# ------------------------------------------------------------
get_msig <- function(collection, subcollection=NULL, species="Homo sapiens") {
  fml <- names(formals(msigdbr::msigdbr))
  
  # msigdbr >= 10 usa collection/subcollection
  if ("collection" %in% fml) {
    if (is.null(subcollection)) {
      msig <- msigdbr::msigdbr(species = species, collection = collection)
    } else {
      msig <- msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
    }
  } else {
    # msigdbr viejo usa category/subcategory
    if (is.null(subcollection)) {
      msig <- msigdbr::msigdbr(species = species, category = collection)
    } else {
      msig <- msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
    }
  }
  
  # Por si acaso: si el paquete devuelve todo y no filtró, filtramos por columnas reales
  # (nombres típicos: gs_collection / gs_subcollection)
  if (!is.null(collection)) {
    if ("gs_collection" %in% colnames(msig)) {
      msig <- msig %>% filter(.data$gs_collection == collection)
    }
  }
  if (!is.null(subcollection)) {
    if ("gs_subcollection" %in% colnames(msig)) {
      msig <- msig %>% filter(.data$gs_subcollection == subcollection)
    } else if ("gs_subcat" %in% colnames(msig)) {
      msig <- msig %>% filter(.data$gs_subcat == subcollection)
    } else if ("subcategory" %in% colnames(msig)) {
      msig <- msig %>% filter(.data$subcategory == subcollection)
    }
  }
  
  msig
}

# ------------------------------------------------------------
# GSEA runner
# ------------------------------------------------------------
run_gsea_one <- function(de_csv, ct, msig_collection="H", msig_subcollection=NULL, out_prefix,
                         minSize=10, maxSize=500, nproc=1, eps=1e-50, top_plot_n=15) {
  
  de <- read_csv(de_csv, show_col_types = FALSE)
  
  # Esperado: gene, log2FC_Healthy_vs_Cirrhosis, pval, FDR
  req <- c("gene","log2FC_Healthy_vs_Cirrhosis","pval","FDR")
  if (!all(req %in% colnames(de))) {
    stop("Faltan columnas en DE: ", paste(setdiff(req, colnames(de)), collapse=", "))
  }
  
  # Ranking: signo(log2FC) * -log10(pval)
  # (positivo => más alto en Healthy; negativo => más alto en Cirrhosis)
  rnk <- de %>%
    transmute(
      gene = as.character(gene),
      stat = sign(log2FC_Healthy_vs_Cirrhosis) * (-log10(pmax(pval, 1e-300)))
    ) %>%
    distinct(gene, .keep_all = TRUE) %>%
    filter(is.finite(stat)) %>%
    arrange(desc(stat), gene)
  
  stats <- rnk$stat
  names(stats) <- rnk$gene
  
  # Gene sets MSigDB
  msig <- get_msig(collection = msig_collection, subcollection = msig_subcollection)
  
  if (!all(c("gs_name","gene_symbol") %in% colnames(msig))) {
    stop("msigdbr no devolvió columnas esperadas gs_name/gene_symbol. Colnames: ",
         paste(colnames(msig), collapse=", "))
  }
  
  pathways <- split(msig$gene_symbol, msig$gs_name)
  
  # FGSEA (forzamos nproc=1 para reducir warnings tipo serialize(...) en Windows)
  fg <- fgsea(
    pathways = pathways,
    stats = stats,
    minSize = minSize,
    maxSize = maxSize,
    eps = eps,
    nproc = nproc
  ) %>% arrange(padj)
  
  out_csv <- file.path(SUM_DIR, paste0(out_prefix, "_", ct, ".csv"))
  write_csv(fg, out_csv)
  
  # Plot simple: top N por padj
  top <- fg %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(top_plot_n)
  if (nrow(top) == 0) {
    message("[WARN] ", ct, " (", out_prefix, "): no hay filas para plot (padj NA). CSV guardado igualmente: ", out_csv)
    return(invisible(fg))
  }
  
  top$pathway <- factor(top$pathway, levels = rev(top$pathway))
  
  p <- ggplot(top, aes(x = NES, y = pathway)) +
    geom_point(aes(size = -log10(padj))) +
    labs(title = paste0("GSEA (", out_prefix, ") — ", ct),
         x = "NES", y = NULL) +
    theme_bw(base_size = 11)
  
  out_png <- file.path(FIG_DIR, paste0(out_prefix, "_", ct, ".png"))
  ggsave(out_png, p, width = 9, height = 6, dpi = 300)
  
  message("[OK] ", ct, " -> ", out_csv, " | ", out_png)
  invisible(fg)
}

# ------------------------------------------------------------
# Run for Level1 celltypes
# ------------------------------------------------------------
all_res <- list()

for (ct in celltypes) {
  de_csv <- file.path(SUM_DIR, paste0("DE_pseudobulk_", ct, "_Healthy_vs_Cirrhosis.csv"))
  if (!file.exists(de_csv)) {
    message("[SKIP] No existe: ", de_csv)
    next
  }
  
  # 1) Hallmark
  fg_h <- run_gsea_one(
    de_csv, ct,
    msig_collection = "H",
    msig_subcollection = NULL,
    out_prefix = "GSEA_MSigDB_Hallmark"
  )
  all_res[[length(all_res) + 1]] <- fg_h %>%
    mutate(celltype = ct, msig = "Hallmark")
  
  # 2) Reactome (C2 CP:REACTOME)
  fg_r <- run_gsea_one(
    de_csv, ct,
    msig_collection = "C2",
    msig_subcollection = "CP:REACTOME",
    out_prefix = "GSEA_MSigDB_Reactome"
  )
  all_res[[length(all_res) + 1]] <- fg_r %>%
    mutate(celltype = ct, msig = "Reactome")
}

# ------------------------------------------------------------
# Summary CSV único (top pathways por celltype, Hallmark + Reactome)
# ------------------------------------------------------------
if (length(all_res) > 0) {
  res_all <- bind_rows(all_res) %>%
    filter(!is.na(padj)) %>%
    mutate(
      direction = if_else(NES >= 0, "Enriched_in_Healthy", "Enriched_in_Cirrhosis")
    )
  
  pick_top <- function(df, n_each_dir = 10, padj_thr = 0.05) {
    sig <- df %>% filter(padj <= padj_thr)
    if (nrow(sig) == 0) {
      # fallback: sin significativos, devolver top por padj
      return(df %>% arrange(padj) %>% head(2 * n_each_dir))
    }
    
    top_pos <- sig %>% filter(NES > 0) %>% arrange(padj) %>% head(n_each_dir)
    top_neg <- sig %>% filter(NES < 0) %>% arrange(padj) %>% head(n_each_dir)
    
    bind_rows(top_pos, top_neg) %>% arrange(padj)
  }
  
  summary_tbl <- res_all %>%
    group_by(celltype, msig) %>%
    group_modify(~ pick_top(.x, n_each_dir = 10, padj_thr = 0.05)) %>%
    ungroup() %>%
    select(celltype, msig, pathway, direction, NES, padj, pval, size) %>%
    arrange(celltype, msig, padj)
  
  out_summary <- file.path(SUM_DIR, "GSEA_MSigDB_summary_topPathways_Level1.csv")
  write_csv(summary_tbl, out_summary)
  
  message("[OK] Summary CSV -> ", out_summary)
} else {
  message("[WARN] No se generó summary: all_res vacío (¿faltan DE_pseudobulk_*?).")
}

