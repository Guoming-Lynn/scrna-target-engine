# -*- coding: utf-8 -*-
"""
Phase 4: Pseudotime Trajectory & Disease Progression Dynamics
============================================================
Publication-grade trajectory analysis:
1. PAGA connectivity matrix heatmap & graph topology.
2. Diffusion Map (15 DCs) & root cell identification in mature subtype.
3. DPT pseudotime computation [0.0, 1.0].
4. Trajectory UMAP (with root star) + Subtype violin with inner boxplot.
5. Case vs Control disease progression dynamics:
   - KDE cell density curves with shaded ribbons.
   - 20-bin binned stacked percentage barplot with external legend.
6. Functional module kinetics:
   - Mean +/- SEM error ribbons along pseudotime.
   - Dual y-axis (twinx) for Target Gene vs Pathway Defense score.
7. Staged phenotypic marker cascade heatmap (50 quantile windows, Z-scores [-2.5, 2.5]).
8. Genome-wide Spearman correlation with pseudotime.
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
from scipy.stats import spearmanr
from scipy.sparse.csgraph import connected_components

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

# Import local publication plotting utilities
from plot_panels import plot_pseudotime_progression_panel, plot_staged_cascade_heatmap
from plot_utils import (
    set_publication_style,
    save_publication_figure,
    COLOR_CASE,
    COLOR_CONTROL,
    PALETTE_6,
    PALETTE_12
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    with open(config_path, "r", encoding="utf-8-sig") as f:
        config = json.load(f)

    TARGET_GENE = config.get("TARGET_GENE", "TARGET")
    PROJECT_NAME = config.get("PROJECT_NAME", "project")
    OUTPUT_ROOT = Path(config.get("OUTPUT_ROOT", "./output"))
    CASE_LABEL = config.get("CASE_LABEL", "Case")
    CONTROL_LABEL = config.get("CONTROL_LABEL", "Control")
    RANDOM_SEED = config.get("RANDOM_SEED", 42)
    PATHWAY_GENE_SETS = config.get("PATHWAY_GENE_SETS", {})

    np.random.seed(RANDOM_SEED)
    sc.settings.seed = RANDOM_SEED
    set_publication_style()

    # Directories
    dirs = {
        'tables': OUTPUT_ROOT / '03_tables' / 'optional_pseudotime',
        'figures': OUTPUT_ROOT / '04_figures' / 'optional_pseudotime',
        'h5ad': OUTPUT_ROOT / '05_h5ad',
        'logs': OUTPUT_ROOT / '99_logs'
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # 1. Load annotated AnnData from Phase 3
    input_h5ad = dirs['h5ad'] / "05_target_subclustered_annotated.h5ad"
    if not input_h5ad.exists():
        # Fallback to general annotated file
        candidate_files = list(dirs['h5ad'].glob("*annotated*.h5ad"))
        if candidate_files:
            input_h5ad = candidate_files[-1]
        else:
            logging.error(f"Cannot find input annotated AnnData in {dirs['h5ad']}")
            sys.exit(1)

    logging.info(f"Loading AnnData from {input_h5ad} ...")
    adata = sc.read_h5ad(input_h5ad)
    logging.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    cluster_key = "subtype" if "subtype" in adata.obs else "smc_celltype_manual"
    if cluster_key not in adata.obs:
        cats = [c for c in adata.obs.columns if c.startswith("leiden")]
        cluster_key = cats[0] if cats else "leiden"

    subtypes = list(adata.obs[cluster_key].cat.categories) if hasattr(adata.obs[cluster_key], 'cat') else list(adata.obs[cluster_key].unique())
    n_subtypes = len(subtypes)
    colors = PALETTE_6 if n_subtypes <= 6 else PALETTE_12[:n_subtypes]
    color_map = {st: c for st, c in zip(subtypes, colors)}

    raw_vars = set(adata.raw.var_names) if adata.raw is not None else set(adata.var_names)
    X_src = adata.raw.X if adata.raw is not None else adata.X
    src_vars = adata.raw.var_names if adata.raw is not None else adata.var_names

    # 2. Neighbors, Diffusion Map & PAGA
    logging.info("Computing neighbors, Diffusion Map (15 DCs) and PAGA topology ...")
    rep = 'X_pca_harmony' if 'X_pca_harmony' in adata.obsm else 'X_pca'
    sc.pp.neighbors(adata, n_neighbors=25, use_rep=rep, metric='cosine', random_state=RANDOM_SEED)

    # Connected components guard
    n_comps, comp_labels = connected_components(adata.obsp['connectivities'], directed=False)
    if n_comps > 1:
        logging.warning(f"Graph has {n_comps} components. Restricting to largest component.")
        largest_comp = pd.Series(comp_labels).value_counts().index[0]
        adata = adata[comp_labels == largest_comp].copy()
        sc.pp.neighbors(adata, n_neighbors=25, use_rep=rep, metric='cosine', random_state=RANDOM_SEED)

    sc.tl.diffmap(adata, n_comps=15, random_state=RANDOM_SEED)
    sc.tl.paga(adata, groups=cluster_key)

    # Plot 1: PAGA Heatmap & Graph
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    paga_mat = adata.uns['paga']['connectivities'].toarray()
    short_labels = [str(st)[:22] for st in subtypes]
    sns.heatmap(paga_mat, cmap='YlOrRd', annot=True, fmt='.2f', xticklabels=short_labels, yticklabels=short_labels, ax=ax)
    ax.set_title("PAGA Connectivity Matrix", fontsize=12, fontweight='bold', pad=10)
    plt.xticks(rotation=35, ha='right')
    save_publication_figure(fig, dirs['figures'] / "01a_PAGA_connectivity_matrix_heatmap")

    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    sc.pl.paga(adata, threshold=0.05, show=False, frameon=True, title="PAGA Graph Topology", ax=ax)
    save_publication_figure(fig, dirs['figures'] / "01b_PAGA_graph_topology")

    # 3. Root Cell Selection & DPT Calculation
    logging.info("Selecting trajectory Root Cell ...")
    # Determine origin subtype: either subtype with highest mature/contractile score, or subtype 0
    score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
    if score_cols:
        root_metric_col = score_cols[0]
        best_local_idx = adata.obs[root_metric_col].argmax()
    else:
        best_local_idx = 0

    root_barcode = adata.obs_names[best_local_idx]
    adata.uns['iroot'] = best_local_idx
    adata.obs['is_DPT_root'] = (adata.obs_names == root_barcode)

    logging.info(f"Root cell selected: {root_barcode}")
    sc.tl.dpt(adata, n_dcs=15)
    adata.obs['dpt_pseudotime'] = np.clip(adata.obs['dpt_pseudotime'].values, 0.0, 1.0)

    # Pseudotime summary table
    pt_stat = adata.obs.groupby(cluster_key, observed=True)['dpt_pseudotime'].agg(
        n_cells='size', mean_pt='mean', median_pt='median', std_pt='std'
    ).reset_index().sort_values('median_pt')
    pt_stat.to_csv(dirs['tables'] / "03_DPT_pseudotime_summary_by_subtype.csv", index=False, encoding='utf-8-sig')

    # Plot 4: DPT UMAP (with root star) + Subtype violin plot
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2))
    scat = axes[0].scatter(
        adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
        c=adata.obs['dpt_pseudotime'], cmap='viridis', s=18, alpha=0.88, rasterized=True
    )
    rx, ry = adata.obsm['X_umap'][best_local_idx]
    axes[0].scatter(rx, ry, color='#d62728', edgecolors='white', s=140, marker='*', zorder=10, label='Root cell')
    cbar = plt.colorbar(scat, ax=axes[0], fraction=0.046, pad=0.03)
    cbar.set_label('DPT Pseudotime', fontsize=10)
    axes[0].set_title("Trajectory: Diffusion Pseudotime (DPT)", fontsize=12.5, fontweight='bold', pad=10)
    axes[0].set_xlabel("UMAP1", fontsize=10.5)
    axes[0].set_ylabel("UMAP2", fontsize=10.5)
    axes[0].legend(loc='lower left', frameon=True, fontsize=9.5)
    axes[0].grid(alpha=0.18, linestyle='--')

    from figure_recipes import ranked_violin
    ranked_violin(adata.obs, cluster_key, 'dpt_pseudotime', axes[1],
                  palette=dict(zip(map(str, subtypes), colors)),
                  title='Pseudotime distribution across subtypes')
    save_publication_figure(fig, dirs['figures'] / "04_UMAP_DPT_pseudotime_and_subtype_violin")

    # 4. Disease Progression Dynamics (KDE + 20-bin Stacked)
    group_col = 'group' if 'group' in adata.obs else 'condition'
    if group_col in adata.obs:
        logging.info("Plotting disease progression dynamics along pseudotime ...")
        plot_pseudotime_progression_panel(
            obs_df=adata.obs,
            pt_col='dpt_pseudotime',
            group_col=group_col,
            case_label=CASE_LABEL,
            control_label=CONTROL_LABEL,
            out_prefix=dirs['figures'] / "05_AC_vs_PA_progression_density_and_binned_ratio",
            n_bins=20
        )

    # 5. Functional Module & Target Gene Kinetics
    logging.info("Computing pseudotime kinetics trends ...")
    bins = np.linspace(0.0, 1.0, 21)
    bin_centers = (bins[:-1] + bins[1:]) / 2.0
    adata.obs['pt_bin'] = pd.cut(adata.obs['dpt_pseudotime'], bins=bins, labels=range(20), include_lowest=True)

    # Extract Target Gene expression
    if TARGET_GENE in src_vars:
        idx_t = src_vars.get_loc(TARGET_GENE)
        col_t = X_src[:, idx_t]
        adata.obs[f'{TARGET_GENE}_expr'] = col_t.toarray().ravel() if sp.issparse(col_t) else np.asarray(col_t).ravel()
    else:
        adata.obs[f'{TARGET_GENE}_expr'] = 0.0

    fig, axes = plt.subplots(1, 2, figsize=(16, 5.8))
    # Left: Module kinetics curves with SEM
    if score_cols:
        for idx_sc, sc_name in enumerate(score_cols[:5]):
            binned = adata.obs.groupby('pt_bin', observed=False)[sc_name].agg(['mean', 'sem']).reindex(range(20)).ffill().bfill()
            mv = binned['mean'].values
            sv = binned['sem'].values
            cl = PALETTE_6[idx_sc % len(PALETTE_6)]
            axes[0].plot(bin_centers, mv, color=cl, linewidth=2.2, label=sc_name.replace('score_', ''))
            axes[0].fill_between(bin_centers, mv - sv, mv + sv, color=cl, alpha=0.15)
        axes[0].set_title("Module Score Dynamics along Trajectory", fontsize=12.5, fontweight='bold', pad=10)
        axes[0].set_xlabel("DPT Pseudotime (Early -> Late)", fontsize=11)
        axes[0].set_ylabel("Standardized Module Score", fontsize=11)
        axes[0].legend(loc='upper right', frameon=True, fontsize=9)
        axes[0].grid(alpha=0.25, linestyle='--')
        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)

    # Right: Target Gene kinetics in Case vs Control + twinx pathway
    if group_col in adata.obs and TARGET_GENE in src_vars:
        case_trend = adata.obs[adata.obs[group_col] == CASE_LABEL].groupby('pt_bin', observed=False)[f'{TARGET_GENE}_expr'].mean().reindex(range(20))
        ctrl_trend = adata.obs[adata.obs[group_col] == CONTROL_LABEL].groupby('pt_bin', observed=False)[f'{TARGET_GENE}_expr'].mean().reindex(range(20))

        v_case = case_trend.notna()
        v_ctrl = ctrl_trend.notna()
        axes[1].plot(bin_centers[v_case], case_trend[v_case], color=COLOR_CASE, linewidth=2.5, marker='o', markersize=4.5, label=f'{TARGET_GENE} in {CASE_LABEL} (Case)')
        axes[1].plot(bin_centers[v_ctrl], ctrl_trend[v_ctrl], color=COLOR_CONTROL, linewidth=2.5, linestyle='--', marker='s', markersize=4.5, label=f'{TARGET_GENE} in {CONTROL_LABEL} (Control)')

        if len(score_cols) > 1:
            ax2 = axes[1].twinx()
            pathway_trend = adata.obs.groupby('pt_bin', observed=False)[score_cols[1]].mean().reindex(range(20)).ffill().bfill()
            ax2.plot(bin_centers, pathway_trend, color='#7570b3', linewidth=2.0, linestyle=':', label=score_cols[1].replace('score_', ''))
            ax2.set_ylabel(score_cols[1].replace('score_', ''), fontsize=10.5, color='#7570b3')
            ax2.tick_params(axis='y', labelcolor='#7570b3')
            ax2.spines['top'].set_visible(False)

        axes[1].set_title(f"{TARGET_GENE} Expression Kinetics along Pseudotime", fontsize=12.5, fontweight='bold', pad=10)
        axes[1].set_xlabel("DPT Pseudotime", fontsize=11)
        axes[1].set_ylabel(f"{TARGET_GENE} Expression", fontsize=11)
        axes[1].legend(loc='upper left', frameon=True, fontsize=9)
        axes[1].grid(alpha=0.25, linestyle='--')
        axes[1].spines['top'].set_visible(False)

    save_publication_figure(fig, dirs['figures'] / "06_pseudotime_module_and_gene_trends")

    # 6. Staged Phenotypic Marker Cascade Heatmap
    logging.info("Generating staged marker cascade heatmap ...")
    if PATHWAY_GENE_SETS:
        plot_staged_cascade_heatmap(
            adata=adata,
            staged_genes_dict=PATHWAY_GENE_SETS,
            pt_col='dpt_pseudotime',
            out_prefix=dirs['figures'] / "07_pseudotime_marker_cascade_heatmap",
            n_windows=50
        )

    # 7. Genome-wide Spearman Correlation with Pseudotime
    logging.info("Computing genome-wide Spearman correlation with pseudotime ...")
    pt_vec = adata.obs['dpt_pseudotime'].values
    corr_results = []
    sample_genes = list(src_vars)[:min(2000, len(src_vars))]
    for g in sample_genes:
        idx_g = src_vars.get_loc(g)
        col_g = X_src[:, idx_g]
        vals = col_g.toarray().ravel() if sp.issparse(col_g) else np.asarray(col_g).ravel()
        if vals.std() > 1e-4:
            rho, pval = spearmanr(pt_vec, vals)
            corr_results.append({'gene': g, 'spearman_rho': rho, 'p_value': pval})

    corr_df = pd.DataFrame(corr_results).dropna().sort_values('spearman_rho', ascending=False)
    corr_df.to_csv(dirs['tables'] / "08_pseudotime_genome_wide_spearman_correlation.csv", index=False, encoding='utf-8-sig')

    # 8. Save AnnData Checkpoint
    final_h5ad = dirs['h5ad'] / "06_target_DPT_pseudotime.h5ad"
    logging.info(f"Saving pseudotime AnnData checkpoint -> {final_h5ad} ...")
    adata.write_h5ad(final_h5ad, compression='gzip')

    # 9. Save JSON Run Report
    report = {
        "step": "04_pseudotime",
        "target_gene": TARGET_GENE,
        "n_cells": adata.n_obs,
        "root_cell": root_barcode,
        "n_subtypes": len(subtypes),
        "pt_min": float(adata.obs['dpt_pseudotime'].min()),
        "pt_max": float(adata.obs['dpt_pseudotime'].max()),
        "top_positive_driver": corr_df.iloc[0]['gene'] if len(corr_df) > 0 else None,
        "top_negative_driver": corr_df.iloc[-1]['gene'] if len(corr_df) > 0 else None
    }
    with open(dirs['logs'] / "04_step04_pseudotime_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    logging.info("[DONE] Phase 4 Trajectory Analysis completed successfully!")


if __name__ == "__main__":
    main()
