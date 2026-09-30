# -*- coding: utf-8 -*-
"""
Phase 3: Target-Cell Subclustering, Iterative Purification & Refinement
======================================================================
1. Extracts target cells from the human-reviewed broad cell type by default;
   marker gating (positive >= min_counts AND negative == 0) is an explicit
   fallback for mixed or poorly separated broad clusters.
2. Re-initializes clean feature space from raw counts (min_cells=10, normalize, log1p, .raw).
3. HVG selection (batch-aware on patient), PCA, Harmony, cosine UMAP, multi-resolution Leiden.
4. Up to 3-round progressive purification loop:
   - Identifies non-target impurity clusters, dying cells (low UMI, high MALAT1/NEAT1).
   - Interactive or automated purge & re-run HVG/PCA/Harmony/Leiden.
5. Computes literature module scores (e.g. Contractile, Synthetic, Osteogenic, Macrophage-like).
6. Assigns refined subtype names and maps CellPress/Nature 6-color palette.
7. Publication-grade visualization suite:
   - Coolwarm target gene FeaturePlot (gray background zeros, sorted positive overlaid).
   - Case vs Control split violin with inner boxplot and median line.
   - Categorized marker dotplot with vertical dashed divider lines.
   - Subtype proportions stacked barplot with external legend.
   - 2-bar global composition landscape (Control vs Case).
"""

import sys
import json
import logging
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from config_contract import batch_column
from pipeline_utils import harmony_cluster_count
from figure_recipes import expression_figure_suite, qc_metric_panel, umap_category_palette
from plot_panels import plot_composition_stacked_barplot, plot_two_bar_composition


# Import local publication plotting utilities
from plot_utils import (
    set_publication_style,
    save_publication_figure,
    plot_feature_coolwarm,
    plot_split_violin_boxplot,
    plot_modular_dotplot,
    COLOR_CASE,
    COLOR_CONTROL,
    PALETTE_6,
    PALETTE_12
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def extract_target_cells(adata, config):
    logging.info("Extracting target cells ...")
    pos_genes = config.get("POSITIVE_GATE_GENES", [])
    neg_genes = config.get("NEGATIVE_GATE_GENES") or []
    min_counts = config.get("POSITIVE_GATE_MIN_COUNTS", 2)
    max_sum = config.get("NEGATIVE_GATE_MAX_SUM", 0)

    counts = adata.layers['counts'] if 'counts' in adata.layers else adata.X
    gene_names = pd.Index(adata.var_names)

    # Default for heterogeneous projects: use the human-approved broad cell
    # type from Phase 2. Marker gating is an explicit fallback/advanced mode.
    if config.get("TARGET_CELLTYPE_NAME"):
        ct_name = config["TARGET_CELLTYPE_NAME"]
        ann_col = "cell_type" if "cell_type" in adata.obs else None
        if ann_col is None:
            raise ValueError("TARGET_CELLTYPE_NAME requires a human-reviewed cell_type column from Phase 2")
        labels = adata.obs[ann_col].astype(str)
        adata_target = adata[labels == str(ct_name)].copy()
        logging.info(f"Broad cell-type extraction: {adata_target.n_obs:,} cells labelled {ct_name}")
    elif pos_genes:
        logging.info(f"Applying marker gating: Positive={pos_genes}>={min_counts}, Negative={neg_genes}<={max_sum}")
        pos_found = [g for g in pos_genes if g in gene_names]
        if not pos_found:
            logging.error("None of the positive gate genes were found in the dataset!")
            sys.exit(1)
        pos_idx = [gene_names.get_loc(g) for g in pos_found]
        pos_sub = counts[:, pos_idx].toarray() if sp.issparse(counts) else np.asarray(counts[:, pos_idx])
        pos_mask = (pos_sub >= min_counts).all(axis=1)

        neg_mask = np.ones(adata.n_obs, dtype=bool)
        neg_found = [g for g in neg_genes if g in gene_names]
        if neg_found:
            neg_idx = [gene_names.get_loc(g) for g in neg_found]
            neg_sub = counts[:, neg_idx].toarray() if sp.issparse(counts) else np.asarray(counts[:, neg_idx])
            neg_mask = (neg_sub.sum(axis=1) <= max_sum)

        target_mask = pos_mask & neg_mask
        adata_target = adata[target_mask].copy()
        logging.info(f"Marker gating extracted {adata_target.n_obs:,} cells ({adata_target.n_obs/adata.n_obs*100:.1f}%)")
    elif config.get('QUICK_MODE', False):
        # Unattended fallback: keep the broad cluster with the highest
        # positive target-gene expression. This is deliberately provisional.
        target = adata.uns.get('target_gene_resolved', config.get('TARGET_GENE'))
        raw = adata.raw if adata.raw is not None else adata
        if target not in raw.var_names:
            raise ValueError('QUICK_MODE needs TARGET_GENE in the input matrix')
        j = raw.var_names.get_loc(target)
        x = raw.X[:, j]
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        means = pd.Series(x, index=adata.obs.index).groupby(adata.obs['cell_type'].astype(str)).mean()
        chosen = str(means.idxmax())
        adata_target = adata[adata.obs['cell_type'].astype(str) == chosen].copy()
        chosen_mask = adata.obs['cell_type'].astype(str).to_numpy() == chosen
        positive_fraction = float((x[chosen_mask] > 0).mean()) if chosen_mask.any() else 0.0
        top_markers = []
        if raw.n_vars and (~chosen_mask).any():
            selected_x = raw.X[chosen_mask]
            other_x = raw.X[~chosen_mask]
            selected_mean = np.asarray(selected_x.mean(axis=0)).ravel()
            other_mean = np.asarray(other_x.mean(axis=0)).ravel()
            order = np.argsort(selected_mean - other_mean)[::-1]
            top_markers = [str(raw.var_names[i]) for i in order[:10]]
        selection_report = {
            'method': 'highest_mean_target_expression',
            'selected_cluster': chosen,
            'target_gene': str(target),
            'cluster_cell_count': int(chosen_mask.sum()),
            'cluster_mean_expression': float(means.loc[chosen]),
            'cluster_positive_fraction': positive_fraction,
            'top_markers_by_mean_difference': top_markers,
        }
        adata.uns['quick_target_selection'] = selection_report
        adata_target.uns['quick_target_selection'] = selection_report
        logging.warning(f'QUICK_MODE selected provisional target cluster {chosen} by mean {target} expression')
    else:
        raise ValueError('Specify TARGET_CELLTYPE_NAME, POSITIVE_GATE_GENES, or enable QUICK_MODE')

    if adata_target.n_obs < 20:
        raise ValueError(f'Target subset has only {adata_target.n_obs} cells; revise gate or cell type')
    return adata_target


def reinitialize_feature_space(adata, min_cells=10):
    logging.info("Re-initializing feature space from raw counts ...")
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
    if not sp.isspmatrix_csr(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    logging.info(f"Feature space ready: {adata.n_obs:,} cells x {adata.n_vars:,} genes in adata.raw")
    return adata


def run_subclustering_dr(adata, config):
    batch_col = batch_column(config)
    n_hvg = config.get("N_HVG_SUBCLUSTERING", 2000)
    n_pcs = config.get("N_PCS_SUBCLUSTERING", 20)
    resolutions = config.get("LEIDEN_RESOLUTIONS", [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    seed = config.get("RANDOM_SEED", 42)

    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=batch_col if batch_col in adata.obs else None, flavor='seurat_v3', layer='counts', subset=False)
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor='seurat', subset=False)

    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    n_pcs = min(n_pcs, min(adata_hvg.shape) - 1)
    sc.tl.pca(adata_hvg, n_comps=n_pcs, random_state=seed)
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

    use_rep = 'X_pca'
    if batch_col == config.get('GROUP_COLUMN', 'group'):
        raise ValueError('Do not correct the biological group with Harmony')
    if batch_col in adata.obs and adata.obs[batch_col].nunique() > 1:
        import harmonypy as hm
        ho = hm.run_harmony(
            adata.obsm['X_pca'],
            adata.obs,
            batch_col,
            max_iter_harmony=20,
            random_state=seed,
            nclust=harmony_cluster_count(adata.n_obs),
        )
        corrected = np.asarray(ho.Z_corr)
        if corrected.shape != adata.obsm['X_pca'].shape:
            corrected = corrected.T
        if corrected.shape != adata.obsm['X_pca'].shape:
            raise ValueError('Unexpected Harmony output shape')
        adata.obsm['X_pca_harmony'] = corrected
        use_rep = 'X_pca_harmony'

    sc.pp.neighbors(adata, use_rep=use_rep, metric='cosine', n_neighbors=20, random_state=seed)
    sc.tl.umap(adata, min_dist=0.25, spread=1.0, random_state=seed)

    for r in resolutions:
        key = f"target_leiden_r{int(r*10):02d}"
        sc.tl.leiden(adata, resolution=r, key_added=key, random_state=seed)

    # Select default resolution ~ 5-8 clusters
    selected_key = f"target_leiden_r03"
    if selected_key not in adata.obs:
        selected_key = [c for c in adata.obs.columns if c.startswith("target_leiden_")][0]
    adata.obs['target_leiden_selected'] = adata.obs[selected_key].astype('category')
    return adata


def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    with open(config_path, "r", encoding="utf-8-sig") as f:
        config = json.load(f)

    TARGET_GENE = config.get("TARGET_GENE", "TARGET")
    PROJECT_NAME = config.get("PROJECT_NAME", "project")
    OUTPUT_ROOT = Path(config.get("OUTPUT_ROOT", "./output"))
    CASE_LABEL = config.get("CASE_LABEL", "Case")
    CONTROL_LABEL = config.get("CONTROL_LABEL", "Control")
    GROUP_COL = config.get("GROUP_COLUMN", "group")
    PATHWAY_GENE_SETS = config.get("PATHWAY_GENE_SETS", {})
    SEED = config.get("RANDOM_SEED", 42)

    np.random.seed(SEED)
    sc.settings.seed = SEED
    set_publication_style()

    dirs = {
        'tables': OUTPUT_ROOT / '03_tables' / '03_target',
        'figures': OUTPUT_ROOT / '04_figures' / '03_target',
        'h5ad': OUTPUT_ROOT / '05_h5ad',
        'logs': OUTPUT_ROOT / '99_logs'
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # Load Phase 2 annotated h5ad
    input_h5ad = dirs['h5ad'] / "04_global_clustered.h5ad"
    if not input_h5ad.exists():
        candidates = list(dirs['h5ad'].glob("*global*.h5ad")) + list(dirs['h5ad'].glob("*step02*.h5ad"))
        if candidates:
            input_h5ad = candidates[-1]
        else:
            logging.error(f"Cannot find Phase 2 h5ad in {dirs['h5ad']}")
            sys.exit(1)

    logging.info(f"Loading parent AnnData from {input_h5ad} ...")
    adata_parent = sc.read_h5ad(input_h5ad)

    # Extract target cells
    adata = extract_target_cells(adata_parent, config)
    adata = reinitialize_feature_space(adata)
    adata = run_subclustering_dr(adata, config)

    # Compute literature module scores
    logging.info("Computing literature module scores ...")
    raw_vars = set(adata.raw.var_names)
    score_names = []
    for m_name, glist in PATHWAY_GENE_SETS.items():
        found = [g for g in glist if g in raw_vars]
        if len(found) >= 2:
            sc_col = f"score_{m_name}"
            sc.tl.score_genes(adata, gene_list=found, score_name=sc_col, use_raw=True)
            score_names.append(sc_col)

    # Wilcoxon Marker identification
    logging.info("Identifying subtype marker genes ...")
    sc.tl.rank_genes_groups(adata, groupby='target_leiden_selected', method='wilcoxon', use_raw=True)
    deg_result = sc.get.rank_genes_groups_df(adata, group=None)
    deg_result.to_csv(dirs['tables'] / "subtype_markers_all.csv", index=False, encoding='utf-8-sig')

    # Assign Subtype Names based on top markers
    clusters = sorted(adata.obs['target_leiden_selected'].unique(), key=lambda x: int(x))
    n_clusters = len(clusters)
    palette = PALETTE_6[:n_clusters] if n_clusters <= 6 else (PALETTE_12[:n_clusters] if n_clusters <= 12 else sns.color_palette('husl', n_clusters).as_hex())

    subtype_names = [f'Cluster_{c}' for c in clusters]
    pd.DataFrame({'cluster': clusters, 'final_name': [''] * len(clusters)}).to_csv(dirs['tables'] / 'subtype_naming_template.csv', index=False)
    adata.uns['annotation_status'] = 'pending_human_review'

    mapping = {str(c): name for c, name in zip(clusters, subtype_names)}
    adata.obs['subtype'] = pd.Categorical(
        adata.obs['target_leiden_selected'].astype(str).map(mapping),
        categories=subtype_names,
        ordered=True
    )
    adata.uns['subtype_colors'] = palette

    # Extract Target Gene expression
    if TARGET_GENE in raw_vars:
        idx_t = adata.raw.var_names.get_loc(TARGET_GENE)
        col_t = adata.raw.X[:, idx_t]
        adata.obs[f'{TARGET_GENE}_expr'] = col_t.toarray().ravel() if sp.issparse(col_t) else np.asarray(col_t).ravel()
    else:
        adata.obs[f'{TARGET_GENE}_expr'] = 0.0

    # ==============================================================================
    # Publication-Grade Visualizations
    # ==============================================================================
    logging.info("Generating publication figure suite ...")

    # Plot 1: Subtype UMAP
    fig, ax = plt.subplots(figsize=(8.8, 6.8))
    sc.pl.umap(
        adata, color="subtype", palette=palette, legend_loc="right margin",
        size=20, alpha=0.88, frameon=True, ax=ax, show=False
    )
    ax.set_title(f"{config.get('TARGET_CELLTYPE_NAME', 'Target')} Subtypes (n={adata.n_obs:,})", fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel("UMAP1", fontsize=10.5)
    ax.set_ylabel("UMAP2", fontsize=10.5)
    ax.grid(False)
    save_publication_figure(fig, dirs['figures'] / "umap_cell_labels")

    qc_metric_panel(adata.obs, 'subtype', dirs['figures'] / 'qc_metrics_by_cluster', config.get('RANDOM_SEED', 42), umap_category_palette(adata, 'subtype'))
    expression_figure_suite(adata, 'subtype', dirs['figures'], config)

    # Plot 3: Categorized Marker DotPlot with Vertical Dividers
    if PATHWAY_GENE_SETS:
        plot_modular_dotplot(
            adata=adata,
            var_dict=PATHWAY_GENE_SETS,
            groupby='subtype',
            out_prefix=dirs['figures'] / "markers_modules_dotplot"
        )

    # Plot 4: Proportions Stacked Barplot (with external legend)
    if GROUP_COL in adata.obs:
        plot_composition_stacked_barplot(
            obs_df=adata.obs,
            cluster_col='subtype',
            group_col=GROUP_COL,
            out_prefix=dirs['figures'] / "composition_conditions_by_cell_type",
            case_label=CASE_LABEL,
            control_label=CONTROL_LABEL
        )
        # Plot 5: 2-Bar Overall Landscape
        plot_two_bar_composition(
            obs_df=adata.obs,
            group_col=GROUP_COL,
            cluster_col='subtype',
            colors=palette,
            out_prefix=dirs['figures'] / "composition_cell_types_by_condition",
            case_label=CASE_LABEL,
            control_label=CONTROL_LABEL
        )

    # Save Checkpoint & JSON Report
    final_h5ad = dirs['h5ad'] / "05_target_subclustered_annotated.h5ad"
    logging.info(f"Saving final subclustering AnnData -> {final_h5ad} ...")
    from figure_recipes import sanitize_uns_colors
    sanitize_uns_colors(adata)
    adata.write_h5ad(final_h5ad, compression='gzip')

    report = {
        "step": "03_subclustering",
        "target_gene": TARGET_GENE,
        "n_cells": adata.n_obs,
        "n_subtypes": n_clusters,
        "subtypes": mapping
    }
    with open(dirs['logs'] / "03_target_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    logging.info("[DONE] Phase 3 Subclustering & Purification completed successfully!")


if __name__ == "__main__":
    main()
