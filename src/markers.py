"""
markers.py — Panel de marcadores canónicos basado en el Inflammation PBMCs Atlas.

Reglas:
- Usar símbolos HGNC en mayúsculas.
- Listas cortas (3–10 genes) por tipo para no sobrecargar dotplots.
- Solo marcadores positivos robustos; la información de “negativo” se deja en comentarios.
- Archivo autocontenido, sin dependencias externas.

Secciones:
- geneMarkers_level1: {Level1: [SYMBOLS...]}
- geneMarkers_level2: {Level1: {Level2: [SYMBOLS...]}}
- GENE_ALIASES: mapeo simple de alias → símbolo canónico
- symbols_to_varnames(adata, symbols): helper para mapear símbolos a var_names
"""

# -------------------------------------------------------------------
# Level 1 (linajes) — derivados de marker_genes_dict["lineages"]
# y del esquema de linajes de tu proyecto.
# -------------------------------------------------------------------

geneMarkers_level1 = {
    # B cells (sin plasmáticas)
    "B": [
        "MS4A1",  # CD20
        "CD79A",
        "CD19",
        "TCL1A",
        "FCER2",
    ],

    # Plasmablasts / plasma cells
    "Plasma": [
        "MZB1",
        "JCHAIN",
        "XBP1",
        "SDC1",   # CD138
        "PRDM1",
        "IRF4",
    ],

    # Monocytes + dendritic cells
    "Mono_and_DC": [
        "LYZ",
        "CST3",
        "CD14",
        "FCGR3A",
        "CLEC9A",
        "CD1C",
    ],

    # Plasmacytoid DC
    "pDC": [
        "IL3RA",
        "LILRA4",
        "IRF7",
        "IRF8",
        "GZMB",
        "JCHAIN",
    ],

    # T + NK + ILC (gran linaje linfoide)
    "T_and_NK": [
        "CD3D",
        "CD4",
        "CD8A",
        "NKG7",
        "GNLY",
        "KLRD1",
    ],

    # Poblaciones “uncertain” / poco conformes con un linaje clásico
    "UTC": [
        "TMSB10",
        "MALAT1",
        "RPL13A",
    ],

    # ILC como linaje separado (aunque biológicamente estén cerca de T/NK)
    "ILC": [
        "KIT",
        "NCR1",
        "KLRG1",
    ],

    # Plaquetas
    "Platelets": [
        "PPBP",
        "PF4",
        "NRGN",
    ],

    # Eritrocitos + progenitores hematopoyéticos
    "RBC_and_HSC": [
        "HBA1",
        "HBB",
        "CD34",
        "KIT",
    ],
}

# -------------------------------------------------------------------
# Level 2 (subtipos dentro de cada linaje)
# Basados directamente en marker_genes_dict del Inflammation Atlas.
# Se han simplificado nombres y listas para adaptarlos al TFM.
# -------------------------------------------------------------------

geneMarkers_level2 = {
    # -----------------------------------------------------------
    # B cells
    # -----------------------------------------------------------
    "B": {
        "B_Naive": [
            "TCL1A",
            "FCER2",
            "IGHD",
            "IGHM",
            "CCR7",
            "SELL",
        ],
        "B_Memory": [
            "CD27",
            "TNFRSF13B",
            "BANK1",
            "IGHG1",
            "IGHA1",
            "CD74",
        ],
        "B_Immature": [
            "CD19",
            "RAG1",
            "RAG2",
            "CD9",
            "SOX4",
        ],
        "B_Atypical": [
            "TBX21",
            "ITGAX",
        ],
        "B_Activated": [
            "CD69",
            "CD83",
            "NFKB1",
            "NFKB2",
        ],
        "B_ISG": [
            "ISG15",
            "IFI6",
            "IFITM1",
            "IFITM2",
        ],
    },

    # -----------------------------------------------------------
    # Plasma cells
    # -----------------------------------------------------------
    "Plasma": {
        "Plasma": [
            "MZB1",
            "SDC1",   # CD138
            "JCHAIN",
            "DERL3",
            "XBP1",
            "PRDM1",
            "IRF4",
        ],
    },

    # -----------------------------------------------------------
    # Monocytes + DC
    # -----------------------------------------------------------
    "Mono_and_DC": {
        "Classical_Mono": [
            "CD14",
            "S100A8",
            "S100A9",
            "LYZ",
            "VCAN",
            "FCN1",
        ],
        "NonClassical_Mono": [
            "FCGR3A",
            "CX3CR1",
            "HLA-DRB1",
            "HLA-DRA",
        ],
        "cDC1": [
            "CLEC9A",
            "XCR1",
            "IDO1",
            "CLNK",
            "ZNF366",
        ],
        "cDC2": [
            "CD1C",
            "FCER1A",
            "CLEC10A",
        ],
        "DC3": [
            "CD1C",
            "S100A8",
            "S100A9",
            "ANXA1",
        ],
        "DC4": [
            "ITGAX",
            "FCGR3A",
            "SERPINA1",
            "LILRB2",
            "SIGLEC10",
        ],
        "DC5": [
            "AXL",
            "SIGLEC6",
            "CD22",
            "DAB2",
        ],
        "aDC": [
            "CCL19",
            "CCR7",
            "IL7R",
            "AIRE",
        ],
    },

    # -----------------------------------------------------------
    # pDC
    # -----------------------------------------------------------
    "pDC": {
        "pDC": [
            "IL3RA",
            "IRF7",
            "LILRA4",
            "IRF8",
            "JCHAIN",
            "GZMB",
        ],
    },

    # -----------------------------------------------------------
    # T + NK linaje
    # Aquí reflejamos lo que pide tu tutora:
    # CD4/CD8 Naive, Effector, Treg, Th, MAIT, γδT, NK, etc.
    # -----------------------------------------------------------
    "T_and_NK": {
        # “Naive/CM” + CD4
        "CD4_Naive": [
            "CD4",
            "CCR7",
            "TCF7",
            "LEF1",
            "SELL",
            "KLF2",
        ],
        # “Naive/CM” + CD8
        "CD8_Naive": [
            "CD8A",
            "CD8B",
            "CCR7",
            "TCF7",
            "LEF1",
            "SELL",
        ],
        # T helper / Th (follicular, Th1/Th17 mix)
        "CD4_Effector_Th": [
            "CXCR3",
            "GATA3",
            "RORC",
            "RORA",
            "IL17A",
            "IL17F",
            "CXCR6",
            "IFNG",
            "CXCR5",
            "CXCL13",
        ],
        # Effector / memory CD8 (citotóxicos)
        "CD8_Effector_Cytotoxic": [
            "GZMK",
            "GZMH",
            "PRF1",
            "NKG7",
            "CCL5",
            "CCL4",
            "ITGAE",
            "KLRG1",
            "CTSW",
        ],
        # Exhausted CD8 / T exhausted
        "Exhausted_T": [
            "HAVCR2",
            "LAG3",
            "PDCD1",
            "TIGIT",
            "TOX",
            "TOX2",
            "LAYN",
            "CTLA4",
        ],
        # T reguladoras
        "Treg": [
            "FOXP3",
            "CTLA4",
            "IL2RA",
            "ICOS",
            "TIGIT",
            "IKZF2",
            "GATA3",
            "CCR7",
        ],
        # γδT
        "GammaDelta_T": [
            "TRGC1",
            "TRGC2",
            "TRDC",
        ],
        # MAIT
        "MAIT": [
            "KLRB1",
            "IL7R",
            "SLC4A10",
        ],
        # NK clásicas (CD16hi / CD56dim, etc.)
        "NK": [
            "NCAM1",
            "FCGR3A",
            "CX3CR1",
            "GNLY",
            "KLRC2",
            "KLRD1",
            "KLRK1",
            "NKG7",
        ],
        # Respuesta interferón (ISG-high)
        "IFN_response_T": [
            "IFI16",
            "IFI35",
            "IFI44",
            "IFI44L",
            "IFI6",
            "IFIH1",
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "IFIT5",
            "ISG15",
        ],
        # Proliferativas (ciclo celular)
        "Proliferative_T": [
            "MKI67",
            "TOP2A",
            "STMN1",
            "UBE2C",
            "PCLAF",
            "CENPF",
            "CDK1",
        ],
    },

    # -----------------------------------------------------------
    # ILC como linaje separado (si aparece como Level1)
    # -----------------------------------------------------------
    "ILC": {
        "ILC": [
            "KIT",
            "NCR1",
            "KLRG1",
        ],
    },

    # -----------------------------------------------------------
    # Platelets
    # -----------------------------------------------------------
    "Platelets": {
        "Platelets": [
            "PPBP",
            "PF4",
            "NRGN",
        ],
    },

    # -----------------------------------------------------------
    # RBC + HSC
    # -----------------------------------------------------------
    "RBC_and_HSC": {
        "RBC": [
            "HBA1",
            "HBB",
        ],
        "HSCs": [
            "CD34",
            "KIT",
        ],
    },

    # -----------------------------------------------------------
    # UTC: tu tutora te ha dicho explícitamente
    # "UTC = MAIT + γδT"
    # Si algunos clusters se anotan como UTC en Level1,
    # estos marcadores permiten identificar qué subtipo son.
    # -----------------------------------------------------------
    "UTC": {
        "MAIT": [
            "KLRB1",
            "IL7R",
            "SLC4A10",
        ],
        "GammaDelta_T": [
            "TRGC1",
            "TRGC2",
            "TRDC",
        ],
    },
}

# -------------------------------------------------------------------
# Aliases → símbolos canónicos (por si en var_names hay anotaciones mixtas)
# -------------------------------------------------------------------

GENE_ALIASES = {
    "SDC1/CD138": "SDC1",
    # puedes añadir más si ves alias raros en tu objeto
}

# -------------------------------------------------------------------
# Helper ultra-ligero para mapear símbolos -> var_names
# -------------------------------------------------------------------

def symbols_to_varnames(adata, symbols):
    """
    Devuelve una lista de var_names presentes para los símbolos pedidos.
    - Si adata.var tiene columna 'symbol', la usa como referencia.
    - Si no, asume que adata.var_names ya son símbolos.
    - Aplica GENE_ALIASES cuando corresponde.
    """
    syms = []
    for s in symbols:
        s = GENE_ALIASES.get(s, s)
        syms.append(s)

    if "symbol" in adata.var.columns:
        # Construimos índice símbolo -> var_name (index)
        m = (
            adata.var["symbol"].dropna()
            .reset_index()
            .drop_duplicates("symbol", keep="first")
            .set_index("symbol")["index"]
            .to_dict()
        )
        return [m[s] for s in syms if s in m]
    else:
        present = set(adata.var_names)
        return [s for s in syms if s in present]
