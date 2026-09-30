# Module Score Gene Sets Reference

This document provides ready-to-use gene sets for `sc.tl.score_genes()` module
scoring, organized by cell type. These are commonly used literature-derived
marker gene sets for major cell lineages encountered in scRNA-seq analyses.

When using these in the pipeline, only genes **present in the dataset** will be
used (genes not found are silently skipped, requiring ≥ 2 matches per module).

---

## Smooth Muscle Cells (SMCs / VSMCs)

```python
SMC_MODULES = {
    "Contractile": ["MYH11", "ACTA2", "TAGLN", "CNN1", "MYL9", "TPM2", "CALD1", "LMOD1", "SYNPO2"],
    "Synthetic_Fibromyocyte": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL14A1", "BGN", "FN1", "VCAN", "TCF21", "AEBP1"],
    "Macrophage_like": ["PTPRC", "CD68", "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "LGALS3", "CTSB", "AIF1"],
    "Foamy": ["APOE", "FABP5", "LPL", "ABCA1", "CD36", "MSR1", "TREM2", "NPC2"],
    "Osteogenic": ["RUNX2", "BMP2", "BMP4", "SPP1", "ALPL", "IBSP", "DLX5", "SOST"],
    "Chondrogenic": ["SOX9", "COL2A1", "ACAN", "COMP", "MGP", "CRTAC1"],
    "Inflammatory": ["CCL2", "CXCL8", "CXCL2", "IL6", "NFKBIA", "JUN", "FOS", "VCAM1", "MMP2"],
}
```

## Macrophages / Monocytes

```python
MACROPHAGE_MODULES = {
    "M1_Proinflammatory": ["TNF", "IL1B", "IL6", "NOS2", "CXCL10", "CCL2", "CD80", "IRF5", "STAT1"],
    "M2_AntiInflammatory": ["MRC1", "CD163", "TGFB1", "IL10", "CCL18", "ARG1", "PPARG", "STAT6"],
    "Foamy_Macrophage": ["APOE", "FABP5", "LPL", "ABCA1", "CD36", "MSR1", "TREM2", "NPC2", "PLIN2"],
    "TREM2_hi": ["TREM2", "GPNMB", "SPP1", "FABP5", "LGALS3", "CD9"],
    "Resident_Tissue": ["C1QA", "C1QB", "C1QC", "LYVE1", "FOLR2", "RNASE1", "F13A1"],
    "Inflammatory_Monocyte": ["S100A8", "S100A9", "S100A12", "FCN1", "VCAN", "CD14"],
}
```

## T Cells

```python
TCELL_MODULES = {
    "Cytotoxic_CD8": ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "IFNG"],
    "Helper_CD4": ["CD4", "IL7R", "LTB", "MAL", "LDHB", "CCR7"],
    "Regulatory": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "IL10"],
    "Exhausted": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"],
    "Memory": ["IL7R", "S100A4", "CCR7", "TCF7", "SELL", "LEF1"],
    "Proliferating": ["MKI67", "TOP2A", "STMN1", "PCNA", "TYMS", "CDK1"],
}
```

## Endothelial Cells

```python
ENDOTHELIAL_MODULES = {
    "Arterial": ["GJA5", "EFNB2", "DLL4", "HEY1", "SOX17", "CXCL12", "BMX"],
    "Venous": ["ACKR1", "PLVAP", "VWF", "SELE", "CLU", "NR2F2"],
    "Lymphatic": ["PROX1", "LYVE1", "PDPN", "FLT4", "CCL21"],
    "Angiogenic": ["CXCR4", "PGF", "VEGFA", "ESM1", "APLN", "ANGPT2"],
    "EndMT": ["ACTA2", "TAGLN", "SNAI1", "SNAI2", "CDH2", "VIM", "FN1"],
}
```

## Fibroblasts

```python
FIBROBLAST_MODULES = {
    "Myofibroblast": ["ACTA2", "TAGLN", "MYL9", "TPM2", "POSTN", "CTGF"],
    "Inflammatory": ["CCL2", "CXCL1", "CXCL2", "IL6", "IL33", "PDPN", "HAS1"],
    "Matrix_Producing": ["COL1A1", "COL1A2", "COL3A1", "FN1", "DCN", "LUM", "BGN"],
    "Lipofibroblast": ["PLIN2", "FABP4", "FABP5", "PPARG", "ADRP", "LPL"],
}
```

## Common Pathway Gene Sets

```python
PATHWAY_MODULES = {
    "Ferroptosis": ["GPX4", "SLC7A11", "ACSL4", "LPCAT3", "TFRC", "FTH1", "FTL",
                    "HMOX1", "NQO1", "SLC3A2", "GCLC", "GCLM", "GSS", "NFE2L2"],
    "NF_kB_Signaling": ["NFKBIA", "NFKB1", "RELA", "TNF", "IL1B", "IKBKB",
                        "TRAF2", "BCL2", "BIRC3", "CXCL8", "CCL2", "ICAM1"],
    "Oxidative_Phosphorylation": ["NDUFA1", "NDUFB1", "SDHA", "UQCRB", "COX7A2",
                                  "ATP5F1A", "ATP5F1B", "ATP5MC1", "CYCS"],
    "Apoptosis": ["BAX", "BAK1", "BCL2", "BCL2L1", "CASP3", "CASP9",
                  "CYCS", "APAF1", "BID", "TP53"],
    "EMT": ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1",
            "ZEB2", "MMP2", "MMP9", "ACTA2"],
    "Cell_Cycle": ["MKI67", "TOP2A", "PCNA", "CDK1", "CDK2", "CCNB1",
                   "CCND1", "CCNE1", "MCM2", "STMN1"],
    "Fatty_Acid_Metabolism": ["FASN", "SCD", "ACACA", "ACLY", "HADHA", "HADHB",
                              "CPT1A", "CPT2", "ACADM", "ACADVL", "ECHS1"],
    "Glycolysis": ["HK1", "HK2", "GPI", "PFKL", "ALDOA", "GAPDH",
                   "PGK1", "ENO1", "PKM", "LDHA"],
}
```

---

## How to Use

In your config.json `PATHWAY_GENE_SETS` field, specify gene sets as:

```json
{
  "PATHWAY_GENE_SETS": {
    "Ferroptosis": ["GPX4", "SLC7A11", "ACSL4", "LPCAT3", "TFRC", "FTH1", "FTL"],
    "NF_kB": ["NFKBIA", "NFKB1", "RELA", "TNF", "IL1B"]
  }
}
```

The pipeline will compute `sc.tl.score_genes()` for each set and include the
scores in all downstream analyses (pseudotime kinetics, subtype comparisons,
sensitivity DEG correlations).
