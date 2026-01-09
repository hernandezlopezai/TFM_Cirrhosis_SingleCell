# PATHWAY ENRICHMENT / GSEA desde DE (CSV) — usando msigdbr + fgsea
# - Input: results/summary_tables/DE_pseudobulk_{ct}_Healthy_vs_Cirrhosis.csv
# - Fallback: summary_tables_final/DE_pseudobulk_{ct}_Healthy_vs_Cirrhosis.csv
# - Output CSV: results/summary_tables/GSEA_* + results/summary_tables/GSEA_MSigDB_summary_topPathways_Level1.csv
# - Output figs: figures/GSEA/GSEA_*

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})


# Repo

find_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = FALSE)
  p <- start
  repeat {
    if (file.exists(file.path(p, ".git")) || file.exists(file.path(p, "README.md"))) {
      return(p)
    }
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

PROJECT_ROOT <- find_project_root(getwd())

RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
FIGURES_DIR <- file.path(PROJECT_ROOT, "figures")

SUM_DIR <- file.path(RESULTS_DIR, "summary_tables")
FIG_OUT_DIR <- file.path(FIGURES_DIR, "GSEA")

LEGACY_SUM_DIR <- file.path(PROJECT_ROOT, "summary_tables_final")
LEGACY_FIG_DIR <- file.path(PROJECT_ROOT, "figures_final") # solo lectura (si hiciera falta)

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUM_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("PROJECT_ROOT  : ", PROJECT_ROOT)
message("SUM_DIR       : ", SUM_DIR)
message("FIG_OUT_DIR   : ", FIG_OUT_DIR)
message("LEGACY_SUM_DIR: ", LEGACY_SUM_DIR)
message("LEGACY_FIG_DIR: ", LEGACY_FIG_DIR)

# Level1 celltypes
celltypes <- c("T", "Mono", "NK", "B", "DC")


# Helpers

pick_first_existing <- function(paths, label = "file") {
  for (p in paths) {
    if (!is.null(p) && file.exists(p)) return(p)
  }
  stop("No encuentro ", label, ". Probé:\n", paste0(" - ", paths, collapse = "\n"))
}

# msigdbr
get_msig <- function(collection, subcollection = NULL, species = "Homo sapiens") {
  fml <- names(formals(msigdbr::msigdbr))

  # msigdbr >= 10 usar collection/subcollection
  if ("collection" %in% fml) {
    if (is.null(subcollection)) {
      msig <- msigdbr::msigdbr(species = species, collection = collection)
    } else {
      msig <- msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
    }
  } else {
    # msigdbr viejo usar category/subcategory
    if (is.null(subcollection)) {
      msig <- msigdbr::msigdbr(species = species, category = collection)
    } else {
      msig <- msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
    }
  }

  # Si el paquete devolvió de más, filtramos por columnas reales típicas
  if (!is.null(collection) && "gs_collection" %in% colnames(msig)) {
    msig <- msig %>% filter(.data$gs_collection == collection)
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

# Normalizar columnas DE
standardize_de_cols <- function(de) {
  gene_col <- intersect(c("gene", "Gene", "gene_symbol", "geneSymbol"), colnames(de))
  if (length(gene_col) == 0) stop("DE: no encuentro columna de gene. Colnames: ", paste(colnames(de), collapse = ", "))
  gene_col <- gene_col[1]

  log2fc_col <- intersect(
    c("log2FC_Healthy_vs_Cirrhosis", "log2FC", "log2fc", "logFC", "avg_log2FC"),
    colnames(de)
  )
  if (length(log2fc_col) == 0) stop("DE: no encuentro columna de log2FC. Colnames: ", paste(colnames(de), collapse = ", "))
  log2fc_col <- log2fc_col[1]

  pval_col <- intersect(c("pval", "p_val", "pvalue", "p_value"), colnames(de))
  if (length(pval_col) == 0) stop("DE: no encuentro columna de pval. Colnames: ", paste(colnames(de), collapse = ", "))
  pval_col <- pval_col[1]

  fdr_col <- intersect(c("FDR", "padj", "p_adj", "p_adjusted"), colnames(de))
  if (length(fdr_col) == 0) stop("DE: no encuentro columna de FDR/padj. Colnames: ", paste(colnames(de), collapse = ", "))
  fdr_col <- fdr_col[1]

  de %>%
    transmute(
      gene = as.character(.data[[gene_col]]),
      log2FC_Healthy_vs_Cirrhosis = suppressWarnings(as.numeric(.data[[log2fc_col]])),
      pval = suppressWarnings(as.numeric(.data[[pval_col]])),
      FDR  = suppressWarnings(as.numeric(.data[[fdr_col]])),
      .keep = "unused"
    )
}

# GSEA runner

run_gsea_one <- function(de_csv, ct,
                         msig_collection = "H",
                         msig_subcollection = NULL,
                         out_prefix,
                         minSize = 10,
                         maxSize = 500,
                         nproc = 1,
                         eps = 1e-50,
                         top_plot_n = 15) {

  de_raw <- readr::read_csv(de_csv, show_col_types = FALSE)
  de <- standardize_de_cols(de_raw)

  # Ranking: signo(log2FC) * -log10(pval)
  # (positivo => más alto en Healthy; negativo => más alto en Cirrhosis)
  rnk <- de %>%
    mutate(
      pval = pmax(pval, 1e-300),
      stat = sign(log2FC_Healthy_vs_Cirrhosis) * (-log10(pval))
    ) %>%
    select(gene, stat) %>%
    distinct(gene, .keep_all = TRUE) %>%
    filter(is.finite(stat)) %>%
    arrange(desc(stat), gene)

  stats <- rnk$stat
  names(stats) <- rnk$gene

  # Gene sets MSigDB
  msig <- get_msig(collection = msig_collection, subcollection = msig_subcollection)
  if (!all(c("gs_name", "gene_symbol") %in% colnames(msig))) {
    stop("msigdbr no devolvió gs_name/gene_symbol. Colnames: ", paste(colnames(msig), collapse = ", "))
  }
  pathways <- split(msig$gene_symbol, msig$gs_name)

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

  # Plot: top N por padj
  top <- fg %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(top_plot_n)
  if (nrow(top) == 0) {
    message("[WARN] ", ct, " (", out_prefix, "): sin filas para plot (padj NA). CSV guardado: ", out_csv)
    return(invisible(fg))
  }

  top$pathway <- factor(top$pathway, levels = rev(top$pathway))

  p <- ggplot(top, aes(x = NES, y = pathway)) +
    geom_point(aes(size = -log10(padj))) +
    labs(title = paste0("GSEA (", out_prefix, ") — ", ct), x = "NES", y = NULL) +
    theme_bw(base_size = 11)

  out_png <- file.path(FIG_OUT_DIR, paste0(out_prefix, "_", ct, ".png"))
  ggsave(out_png, p, width = 9, height = 6, dpi = 300)

  message("[OK] ", ct, " -> ", out_csv, " | ", out_png)
  invisible(fg)
}


# Run for Level1 celltypes

all_res <- list()

for (ct in celltypes) {
  de_csv <- pick_first_existing(
    c(
      file.path(SUM_DIR, paste0("DE_pseudobulk_", ct, "_Healthy_vs_Cirrhosis.csv")),
      file.path(LEGACY_SUM_DIR, paste0("DE_pseudobulk_", ct, "_Healthy_vs_Cirrhosis.csv"))
    ),
    label = paste0("DE_pseudobulk_", ct, "_Healthy_vs_Cirrhosis.csv")
  )

  # Hallmark
  fg_h <- run_gsea_one(
    de_csv, ct,
    msig_collection = "H",
    msig_subcollection = NULL,
    out_prefix = "GSEA_MSigDB_Hallmark",
    nproc = 1
  )
  all_res[[length(all_res) + 1]] <- fg_h %>%
    mutate(celltype = ct, msig = "Hallmark")

  # Reactome
  fg_r <- run_gsea_one(
    de_csv, ct,
    msig_collection = "C2",
    msig_subcollection = "CP:REACTOME",
    out_prefix = "GSEA_MSigDB_Reactome",
    nproc = 1
  )
  all_res[[length(all_res) + 1]] <- fg_r %>%
    mutate(celltype = ct, msig = "Reactome")
}


# Summary CSV único

if (length(all_res) > 0) {
  res_all <- bind_rows(all_res) %>%
    filter(!is.na(padj)) %>%
    mutate(direction = if_else(NES >= 0, "Enriched_in_Healthy", "Enriched_in_Cirrhosis"))

  pick_top <- function(df, n_each_dir = 10, padj_thr = 0.05) {
    sig <- df %>% filter(padj <= padj_thr)
    if (nrow(sig) == 0) {
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