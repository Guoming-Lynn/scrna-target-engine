"""Composition and trajectory panels moved out of plot_utils."""

from pathlib import Path
from typing import Dict, List, Tuple, Union
import logging

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from plot_utils import (
    CMAP_DIVERGING,
    COLOR_CASE,
    COLOR_CONTROL,
    save_publication_figure,
    set_publication_style,
)


# ==============================================================================
# 4. Proportions Stacked Barplots (Legend Positioned Outside)
# ==============================================================================

def plot_composition_stacked_barplot(
    obs_df: pd.DataFrame,
    cluster_col: str,
    group_col: str,
    out_prefix: Union[Path, str],
    case_label: str = "AC",
    control_label: str = "PA",
    title: str = "Subtype Composition by Condition"
):
    """
    Stacked percentage barplot with legend placed outside to prevent bar overlap.
    """
    set_publication_style()
    out_prefix = Path(out_prefix)
    crosstab = pd.crosstab(obs_df[cluster_col], obs_df[group_col], normalize='index') * 100

    cols = [c for c in [case_label, control_label] if c in crosstab.columns]
    if len(cols) < 2:
        cols = list(crosstab.columns)

    palette = [COLOR_CASE, COLOR_CONTROL] if cols == [case_label, control_label] else None

    fig, ax = plt.subplots(figsize=(8.8, 5.5))
    crosstab[cols].plot(
        kind='bar',
        stacked=True,
        color=palette,
        edgecolor='black',
        linewidth=0.6,
        ax=ax
    )
    ax.set_title(title, fontsize=12.5, fontweight='bold', pad=10)
    ax.set_xlabel("Subtype / Cluster", fontsize=11)
    ax.set_ylabel("Percentage (%)", fontsize=11)
    ax.set_xticklabels(crosstab.index, rotation=30, ha='right', fontsize=9.5)
    ax.legend(
        [f"{c} (Case)" if c == case_label else f"{c} (Control)" for c in cols],
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        frameon=True,
        fontsize=9.5
    )
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_publication_figure(fig, out_prefix)


def plot_two_bar_composition(
    obs_df: pd.DataFrame,
    group_col: str,
    cluster_col: str,
    colors: List[str],
    out_prefix: Union[Path, str],
    case_label: str = "AC",
    control_label: str = "PA",
    title: str = "Overall Landscape: Control vs Case"
):
    """
    2-bar stacked composition comparing Control vs Case using exact UMAP colors.
    """
    set_publication_style()
    out_prefix = Path(out_prefix)
    comp = pd.crosstab(obs_df[group_col], obs_df[cluster_col], normalize='index') * 100

    order = [c for c in [control_label, case_label] if c in comp.index]
    if len(order) >= 2:
        comp = comp.reindex(order)

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    comp.plot(
        kind='bar',
        stacked=True,
        color=colors[:comp.shape[1]],
        edgecolor='black',
        linewidth=0.7,
        width=0.45,
        ax=ax
    )
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel("Condition", fontsize=11.5, fontweight='bold')
    ax.set_ylabel("Proportion (%)", fontsize=11.5, fontweight='bold')
    ax.set_xticklabels([f"{c} (Control)" if c == control_label else f"{c} (Case)" for c in comp.index], rotation=0, fontsize=10.5, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.legend(title="Subtypes", bbox_to_anchor=(1.03, 1), loc='upper left', frameon=True, fontsize=9.5, title_fontsize=10.5)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_publication_figure(fig, out_prefix)


# ==============================================================================
# 5. Pseudotime Disease Progression (KDE Density + 20-Bin Stacked Ratio)
# ==============================================================================

def plot_pseudotime_progression_panel(
    obs_df: pd.DataFrame,
    pt_col: str,
    group_col: str,
    case_label: str,
    control_label: str,
    out_prefix: Union[Path, str],
    n_bins: int = 20
):
    """
    Two-panel disease progression along pseudotime:
    - Left: Kernel Density Estimation (KDE) distribution curves with shaded area.
    - Right: n_bins stacked barplot of Case vs Control proportion.
    """
    set_publication_style()
    out_prefix = Path(out_prefix)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

    # Left: KDE Density
    ctrl_data = obs_df[obs_df[group_col] == control_label]
    case_data = obs_df[obs_df[group_col] == case_label]

    sns.kdeplot(
        data=ctrl_data,
        x=pt_col,
        color=COLOR_CONTROL,
        label=f'{control_label} (Control, n={len(ctrl_data):,})',
        fill=True,
        alpha=0.35,
        linewidth=2.2,
        ax=axes[0]
    )
    sns.kdeplot(
        data=case_data,
        x=pt_col,
        color=COLOR_CASE,
        label=f'{case_label} (Case, n={len(case_data):,})',
        fill=True,
        alpha=0.35,
        linewidth=2.2,
        ax=axes[0]
    )
    axes[0].set_title("Progression Density along Trajectory", fontsize=12.5, fontweight='bold', pad=10)
    axes[0].set_xlabel("Pseudotime", fontsize=11)
    axes[0].set_ylabel("Cell Density", fontsize=11)
    axes[0].legend(loc='upper right', frameon=True, fontsize=10)
    axes[0].grid(alpha=0.25, linestyle='--')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # Right: Binned Proportions
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2.0
    pt_binned = pd.cut(obs_df[pt_col], bins=bins, labels=range(n_bins), include_lowest=True)
    counts = pd.crosstab(pt_binned, obs_df[group_col]).reindex(range(n_bins), fill_value=0)
    totals = counts.sum(axis=1).replace(0, np.nan)
    pct_case = (counts[case_label] / totals * 100).fillna(0).values if case_label in counts else np.zeros(n_bins)
    pct_ctrl = (counts[control_label] / totals * 100).fillna(0).values if control_label in counts else np.zeros(n_bins)

    bin_df = pd.DataFrame({'pct_case': pct_case, 'pct_ctrl': pct_ctrl}, index=bin_centers.round(2).astype(str))
    bin_df.plot(
        kind='bar',
        stacked=True,
        color=[COLOR_CASE, COLOR_CONTROL],
        edgecolor='black',
        linewidth=0.5,
        width=0.85,
        ax=axes[1]
    )
    axes[1].set_title(f"{case_label} vs {control_label} across {n_bins} Pseudotime Bins", fontsize=12.5, fontweight='bold', pad=10)
    axes[1].set_xlabel("Pseudotime Bin (Early -> Late)", fontsize=11)
    axes[1].set_ylabel("Percentage (%)", fontsize=11)
    axes[1].legend([f"{case_label} (Case)", f"{control_label} (Control)"], bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, fontsize=10)
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha='right', fontsize=8.5)
    axes[1].grid(axis='y', linestyle='--', alpha=0.3)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    save_publication_figure(fig, out_prefix)


# ==============================================================================
# 6. Staged Phenotypic Marker Cascade Heatmap
# ==============================================================================

def plot_staged_cascade_heatmap(
    adata,
    staged_genes_dict: Dict[str, List[str]],
    pt_col: str,
    out_prefix: Union[Path, str],
    n_windows: int = 50,
    zscore_clip: Tuple[float, float] = (-2.5, 2.5)
):
    """
    Phenotypic transition cascade heatmap:
    - 50 quantile windows along pseudotime.
    - Z-score normalized expression clipped to [-2.5, 2.5].
    - Left-side color bar marking biological stages.
    - 'coolwarm' colormap.
    """
    set_publication_style()
    out_prefix = Path(out_prefix)
    raw_vars = set(adata.raw.var_names) if adata.raw is not None else set(adata.var_names)

    all_genes = []
    gene_to_category = {}
    for cat_name, glist in staged_genes_dict.items():
        valid = [g for g in glist if g in raw_vars]
        for g in valid:
            all_genes.append(g)
            gene_to_category[g] = cat_name

    if not all_genes:
        logging.warning("No genes found for cascade heatmap")
        return

    pt_series = adata.obs[pt_col]
    window_labels = pd.qcut(pt_series, q=n_windows, labels=False, duplicates='drop')

    # Extract expression matrix
    X_src = adata.raw.X if adata.raw is not None else adata.X
    src_vars = adata.raw.var_names if adata.raw is not None else adata.var_names
    idx_genes = [src_vars.get_loc(g) for g in all_genes]

    sub_X = X_src[:, idx_genes]
    dense_X = sub_X.toarray() if sp.issparse(sub_X) else np.asarray(sub_X)
    df = pd.DataFrame(dense_X, columns=all_genes, index=adata.obs_names)
    df['pt_window'] = window_labels.values

    window_means = df.groupby('pt_window').mean()
    zscore_mat = window_means.apply(lambda col: (col - col.mean()) / (col.std() if col.std() > 0 else 1.0), axis=0).T
    zscore_mat = zscore_mat.clip(zscore_clip[0], zscore_clip[1])

    # Plot
    fig = plt.figure(figsize=(14.5, max(8.0, len(all_genes) * 0.28)))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.03, 0.97], wspace=0.04)

    cat_names = list(staged_genes_dict.keys())
    cat_palette = sns.color_palette("Set2", n_colors=len(cat_names)).as_hex()
    cat_map = {c: col for c, col in zip(cat_names, cat_palette)}

    # Category Bar
    ax_cat = fig.add_subplot(gs[0, 0])
    cat_indices = np.array([cat_palette.index(cat_map[gene_to_category[g]]) for g in all_genes])[:, None]
    ax_cat.imshow(cat_indices, aspect='auto', cmap=mcolors.ListedColormap(cat_palette))
    ax_cat.set_xticks([])
    ax_cat.set_yticks([])
    for sp_name in ax_cat.spines.values():
        sp_name.set_visible(False)

    # Main Heatmap
    ax_hm = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        zscore_mat,
        cmap=CMAP_DIVERGING,
        cbar_kws={'label': 'Z-scored Expression', 'fraction': 0.03, 'pad': 0.02},
        ax=ax_hm,
        rasterized=True
    )
    ax_hm.set_title("Staged Phenotypic Marker Cascade along Trajectory", fontsize=13, fontweight='bold', pad=12)
    ax_hm.set_xlabel(f"Pseudotime Window (Early -> Late, {n_windows} Quantiles)", fontsize=11)
    ax_hm.set_ylabel("Phenotypic Marker Genes", fontsize=11)
    ax_hm.set_yticklabels(all_genes, rotation=0, fontsize=8.5)
    ax_hm.set_xticks(np.linspace(0, n_windows - 1, 6))
    ax_hm.set_xticklabels([f"{v:.1f}" for v in np.linspace(0.0, 1.0, 6)], rotation=0)

    # Add category legend
    legend_handles = [plt.Line2D([0], [0], color=cat_map[c], lw=5, label=c) for c in cat_names]
    ax_hm.legend(handles=legend_handles, bbox_to_anchor=(1.12, 1), loc='upper left', frameon=True, fontsize=9.5, title="Phenotype Stage")

    save_publication_figure(fig, out_prefix)


