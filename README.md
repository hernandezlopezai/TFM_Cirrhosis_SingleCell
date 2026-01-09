# Exploring Circulating Immune Profile Remodeling in Cirrhosis: Reproducible scRNA-seq Pipeline

Repository containing a reproducible scRNA-seq analysis pipeline for an immune (PBMC) cohort comparing Healthy vs Cirrhosis.

Public repository: `github.com/hernandezlopezai/TFM_Cirrhosis_SingleCell`

---

## Overview

1) QC -> normalization -> HVG -> PCA  
2) Integration / batch correction (Harmony)  
3) Global clustering and Level1 annotation + lineage split  
4) Lineage-specific subclustering (Level2) for lineages with true Level2 resolution  
5) Global assembly (Level1 + Level2), summary tables/figures  
6) Final downstream steps:
   - RBC-out filtering (including robust exclusion of `RBC_and_HSC` when applicable)
   - Construction of `Level1_refined` from Level2 (T vs NK; Mono vs DC; DC3->DC hotfix)
   - `Level2_final` (minimal post-processing): `Conv_T_other -> CD4_Memory` (without converting NaN -> "nan")
   - Final UMAP (Harmony), default batch = `patientID`
   - QA (coherence, mixing, numeric dotplot matrix)
   - Composition (counts/props) + scCODA (`Level1_refined` and `Level2_final`)
   - Pseudobulk DE + volcano plots + tables
   - Pathway enrichment (GSEA / MSigDB Hallmark + Reactome) in R
   - Annex: lineage UMAPs + `Level2_final` dotplots (2 markers)

---

## Repository structure

- `notebooks/`  
  Pipeline notebooks (01->20), runnable in order.

- `src/`  
  Reusable repository code. Includes `src/paths.py` (paths) and `src/markers.py` (marker panels + helpers).

- `config/`  
  Configuration: `config/config.yaml` + JSON maps (e.g. `level2_map.json`).

- `R/`  
  R scripts (e.g. `pathway_enrichment.R` for GSEA).

- `data/`  
  Inputs (large data are not versioned).

- `results/`  
  Pipeline-generated artifacts/tables.

- `figures/`  
  Pipeline-generated figures.

---

## Repository conventions

Paths:
- Always resolved via `src/paths.py`.
- Typical pattern in notebooks:
  ```python
  from src.paths import project_paths
  from pathlib import Path

  P = project_paths(Path.cwd())
  PROJECT_ROOT = P.PROJECT_ROOT
  RESULTS_DIR = P.RESULTS_DIR
  FIGURES_DIR = P.FIGURES_DIR
  ```

Config:
- Parameters live in `config/config.yaml` (flat YAML).

Outputs:
- Tables/artifacts: `results/`
- Figures: `figures/`
- Lineages:
- `results/lineages/level1/`
- `results/lineages/level2/`
- `results/markers/level2/`
- `figures/level2/`

Fixed downstream decisions:
- `Level2_final` is used downstream (does not replace Level2): minimal mapping `Conv_T_other -> CD4_Memory`, stored as JSON under `results/summary_tables/...`.
- NB13 computes the “final” Harmony by default using batch = `patientID` (optional override in config: `umap_final_harmony_batch_key`).

---

## Installation (Conda)

If `environment.yml` is available:
```bash
conda env create -f environment.yml
conda activate tfm-cirrhosis
```

If you need to generate it from the current environment (from the repository root):
```bash
conda activate tfm-cirrhosis
conda env export --from-history --no-builds > environment.yml
```

Reproducibility note:
- `--from-history` exports packages explicitly installed via conda, but may not capture packages installed via `pip`.
- If you used `pip install ...`, add a `pip:` block to `environment.yml` or save a `requirements-pip.txt`.

---

## Requirements

Python:
- Scanpy stack: `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`, etc.

R (for GSEA):
- `dplyr`, `readr`, `msigdbr`, `fgsea`, `ggplot2`

---

## Data

This repository does not include raw data or large objects. Place the required inputs under `data/` as defined in `config/config.yaml`.

---

## Running the pipeline

From the repository root:

```bash
jupyter lab
```


Run the notebooks in order (01->20):

1. `01_...` overview / setup  
2. `02_...` QC  
3. `03_...` normalization + HVG  
4. `04_...` Harmony embedding  
5. `05_...` neighbors + UMAP + L1 clustering  
6. `06_...` Level1 annotation + split  
7. `07_...` Level2 embedding/clustering per lineage  
8. `08_...` markers + DE per lineage  
9. `09_...` global assembly + plots/tables  
10. `10_...` final cleanup + Level1_refined + RBC-out  
11. `11_...` Conv_T_other / `Level2_final_map.json`  
12. `12_...` final global dotplot (Level2_final) + numeric matrix  
13. `13_...` final UMAP (Harmony) by `patientID`  
14. `14_...` QA coherence / HSCs  
15. `15_...` composition + boxplots  
16. `16_...` scCODA Level1_refined  
17. `17_...` scCODA Level2_final  
18. `18_...` pseudobulk DE + volcano plots  
19. `19_...` DE supplementary closure  
20. `20_...` ANNEX (lineage UMAPs + Level2_final dotplots)

---

## Pathway enrichment

The script lives in `R/pathway_enrichment.R`:
- Input: `results/summary_tables/DE_pseudobulk_{ct}_Healthy_vs_Cirrhosis.csv`
- Output CSV: `results/summary_tables/GSEA_*` and `results/summary_tables/GSEA_MSigDB_summary_topPathways_Level1.csv`
- Output figures: `figures/GSEA/`

Run from the repository root:

```bash
Rscript R/pathway_enrichment.R
```

---

## Main outputs

- Summary tables: `results/summary_tables/`
- Figures: `figures/`
- Lineages: `results/lineages/` and `figures/level2/`
- Annex: `figures/annex_lineage_umap_dotplots/`
- If present: `final_deliverables/` for export-ready materials (tables/PNGs/PDFs).

---

## Reproducibility

- The pipeline writes outputs to `results/` and `figures/`.
- To reproduce results, run notebooks in order (01->20) and use the parameters in `config/config.yaml`.
- Input data are not distributed in this repository; they must be placed under `data/` according to the configuration.

---

## License

This repository is distributed under the MIT License. See `LICENSE` for the full text.
