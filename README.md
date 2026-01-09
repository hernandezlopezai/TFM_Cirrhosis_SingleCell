# Exploring Circulating Immune Profile Remodeling in Cirrhosis: Reproducible scRNA-seq Pipeline

Repositorio del pipeline de análisis scRNA-seq para un cohorte inmune (PBMC) comparando Healthy vs Cirrhosis.

Repositorio público: `github.com/hernandezlopezai/TFM_Cirrhosis_SingleCell`

---

## Qué hace el pipeline (resumen)

1) QC -> normalización -> HVG -> PCA  
2) Integración / batch correction (Harmony)  
3) Clustering global y anotación Level1 + split por linajes  
4) Subclustering por linaje (Level2) en linajes con Level2 real
5) Ensamblaje global (Level1 + Level2), tablas/figuras globales  
6) Downstream final:
   - Filtrado RBC-out (incluye exclusión robusta de `RBC_and_HSC` cuando aplica)
   - Construcción de `Level1_refined` desde Level2 (T vs NK; Mono vs DC; hotfix DC3->DC)
   - `Level2_final` (post-procesado mínimo): `Conv_T_other -> CD4_Memory` (sin convertir NaN->"nan")
   - UMAP final (Harmony) por `patientID` por defecto
   - QA (coherencia, mixing, dotplot numérico)
   - Composición (counts/props) + scCODA (Level1_refined y Level2_final)
   - DE pseudobulk + volcanos + tablas
   - Pathway enrichment (GSEA / MSigDB Hallmark + Reactome) en R
   - Annex: UMAPs por linaje + dotplots Level2_final (2 marcadores)

---

## Estructura del repo

- `notebooks/`  
  Notebooks del pipeline (01->20), ejecutables en orden.

- `src/`  
  Código reutilizable del repo. Incluye `src/paths.py` (rutas) y `src/markers.py` (paneles + helpers).

- `config/`  
  Configuración: `config/config.yaml` + mapas JSON (p.ej. `level2_map.json`).

- `R/`  
  Scripts en R (p.ej. `pathway_enrichment.R` para GSEA).

- `data/`  
  Inputs (no se versionan datos grandes).

- `results/`  
  Artefactos/tablas generados por el pipeline.

- `figures/`  
  Figuras generadas por el pipeline.

---

## Convenciones del repo

Rutas:
- Se resuelven SIEMPRE con `src/paths.py`.
- Patrón típico en notebooks:
  ```python
  from src.paths import project_paths
  from pathlib import Path

  P = project_paths(Path.cwd())
  PROJECT_ROOT = P.PROJECT_ROOT
  RESULTS_DIR = P.RESULTS_DIR
  FIGURES_DIR = P.FIGURES_DIR
  ```

Config:
- Parámetros en `config/config.yaml` (YAML plano).

Outputs:
- Tablas/artefactos: `results/`
- Figuras: `figures/`
- Linajes:
  - `results/lineages/level1/`
  - `results/lineages/level2/`
  - `results/markers/level2/`
  - `figures/level2/`

Decisiones downstream fijadas:
- `Level2_final` se usa downstream (no reemplaza Level2): mapping mínimo `Conv_T_other -> CD4_Memory`, almacenado como JSON en `results/summary_tables/...`.
- NB13 hace Harmony “final” por defecto con batch=`patientID` (override opcional en config: `umap_final_harmony_batch_key`).

---

## Instalación (Conda)

Si ya se dispone de `environment.yml`:
```bash
conda env create -f environment.yml
conda activate tfm-cirrhosis
```

Si se necesita generar desde el entorno actual (en la raíz del repositorio):
```bash
conda activate tfm-cirrhosis
conda env export --from-history --no-builds > environment.yml
```

Nota sobre reproducibilidad:
- `--from-history` exporta lo “instalado explícitamente” por conda, pero puede no capturar paquetes instalados por `pip`.
- Si se usa `pip install ...`, hayq ue añadir un bloque `pip:` en `environment.yml` o guardar un `requirements-pip.txt`.

---

## Requisitos

Python:
- Scanpy stack: `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`, etc.

R (para GSEA):
- `dplyr`, `readr`, `msigdbr`, `fgsea`, `ggplot2`

---

## Datos

Este repositorio no incluye datos crudos ni objetos pesados. Coloca los inputs en `data/` según lo definido en `config/config.yaml`.

---

## Ejecución del pipeline

Desde la raíz del repo:
```bash
jupyter lab
```

Ejecutar los notebooks en orden (01->20):

1. `01_...` overview / setup  
2. `02_...` QC  
3. `03_...` normalization + HVG  
4. `04_...` Harmony embedding  
5. `05_...` neighbors + UMAP + clustering L1  
6. `06_...` Level1 annotation + split  
7. `07_...` Level2 embedding/clustering por linaje  
8. `08_...` markers + DE por linaje  
9. `09_...` assembly global + plots/tablas  
10. `10_...` limpieza final + Level1_refined + RBC-out  
11. `11_...` Conv_T_other / Level2_final_map.json  
12. `12_...` dotplot global final (Level2_final) + matriz numérica  
13. `13_...` UMAP final (Harmony) por `patientID`  
14. `14_...` QA coherence / HSCs  
15. `15_...` composición + boxplots  
16. `16_...` scCODA Level1_refined  
17. `17_...` scCODA Level2_final  
18. `18_...` pseudobulk DE + volcanos  
19. `19_...` cierre suplementaria DE  
20. `20_...` ANNEX (UMAPs por linaje + dotplots Level2_final)

---

## Pathway enrichment

El script está en `R/pathway_enrichment.R`:
- Input: `results/summary_tables/DE_pseudobulk_{ct}_Healthy_vs_Cirrhosis.csv`
- Output CSV: `results/summary_tables/GSEA_*` y `results/summary_tables/GSEA_MSigDB_summary_topPathways_Level1.csv`
- Output figs: `figures/GSEA/`

Ejecutar desde la raíz del repo:
```bash
Rscript R/pathway_enrichment.R
```

---

## Outputs principales

- Tablas resumen: `results/summary_tables/`
- Figuras: `figures/`
- Linajes: `results/lineages/` y `figures/level2/`
- Annex: `figures/annex_lineage_umap_dotplots/`
- Si existe: `final_deliverables/` para materiales “pegables” (tablas/PNGs/PDFs).

---

## Reproducibilidad

- El pipeline escribe outputs en `results/` y `figures/`.
- Para reproducir resultados, ejecuta los notebooks en orden (01→20) y usa los parámetros de `config/config.yaml`.
- Los datos de entrada no se distribuyen en este repositorio; deben colocarse en `data/` según la configuración.

---

## Licencia

Este repositorio se distribuye bajo la licencia MIT. Ver el archivo `LICENSE` para el texto completo.



