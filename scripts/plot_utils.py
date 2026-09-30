# -*- coding: utf-8 -*-
"""
scRNA-seq Target-Gene Pipeline: Publication-Grade Aesthetic Plotting Utilities
=============================================================================
Encapsulates all visual styling, color palettes, and specialized plotting functions
for consistent, reviewable figures across pipeline runs.

Styles & Standards:
- Font: DejaVu Sans / Arial / Helvetica, 300 DPI PNG + vector PDF.
- Spines: Top/Right spines removed or subtle #555555.
- Colors: Case = #E69F00 (orange), Control = #0072B2 (blue).
- CellPress / Nature 6-color & 12-color palettes.
- FeaturePlot: Zero-expression cells in background gray (#E6E6E6), positive cells
  sorted and plotted on top with 'coolwarm', 99th percentile cutoff.
- Split Violin: Case vs Control split with inner boxplot and median line.
- Modular DotPlot: 'coolwarm' colormap, black dot edges, vertical dashed separators.
- Cascade Heatmap: 50 quantile windows, Z-score [-2.5, 2.5], side stage color bar.
- Trajectory: KDE density + 20-bin stacked proportions + dual-axis kinetics ribbons.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import logging

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from figure_output import figure_path, finalize_metadata, write_figure_index


# ==============================================================================
# 0. Global Palettes & Aesthetic Standards
# ==============================================================================

# Condition contrast
COLOR_CASE = "#E69F00"       # Okabe-Ito orange (Case)
COLOR_CONTROL = "#0072B2"    # Okabe-Ito blue (Control)
PALETTE_CONDITION = {"Case": COLOR_CASE, "Control": COLOR_CONTROL}

# Standard 6-cluster palette (CellPress / Nature classic)
PALETTE_6 = [
    "#4374B3",  # Steel Blue
    "#E25A53",  # Coral Red
    "#73A973",  # Sage Green
    "#E8963D",  # Ochre Orange
    "#9467BD",  # Muted Purple
    "#D95F02"   # Burnt Orange
]

# Standard 12-cluster palette
PALETTE_12 = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
    "#2ca02c", "#98df8a", "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5", "#8c564b", "#e377c2"
]

# Diverging & Sequential Colormaps
CMAP_DIVERGING = "coolwarm"
CMAP_SEQUENTIAL = "viridis"
CMAP_MODULE = "YlOrRd"


def set_publication_style():
    """Apply global matplotlib aesthetic settings for Cell/Nature-grade output."""
    from matplotlib import font_manager
    installed = {f.name for f in font_manager.fontManager.ttflist}
    fallbacks = [name for name in ['Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC'] if name in installed]
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica'] + fallbacks
    plt.rcParams['axes.edgecolor'] = '#333333'
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['axes.labelsize'] = 10.5
    plt.rcParams['xtick.labelsize'] = 9.5
    plt.rcParams['ytick.labelsize'] = 9.5
    plt.rcParams['legend.fontsize'] = 9.5
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['figure.titlesize'] = 13.5
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'
    mpl.rcParams.update({'pdf.fonttype': 42, 'svg.fonttype': 'none'})


def save_publication_figure(fig: plt.Figure, out_prefix: Union[Path, str], dpi: int = 300, svg: bool = False):
    """
    Save figure simultaneously as PNG (300 DPI) and vector PDF with tight layout.
    Cleans up figure memory after saving.
    """
    out_prefix = figure_path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    finalize_metadata(fig, out_prefix)
    fig.savefig(out_prefix.with_suffix('.png'), dpi=dpi, bbox_inches='tight')
    fig.savefig(out_prefix.with_suffix('.pdf'), bbox_inches='tight')
    if svg:
        fig.savefig(out_prefix.with_suffix('.svg'), bbox_inches='tight')
    audits = [ax._plot_audit for ax in fig.axes if hasattr(ax, '_plot_audit')]
    if audits:
        import json
        out_prefix.with_suffix('.plot.json').write_text(json.dumps(audits, indent=2, ensure_ascii=False), encoding='utf-8')
    plt.close(fig)
    write_figure_index(out_prefix.parent)


# ==============================================================================
# 1. Advanced FeaturePlot (Coolwarm with Background-Zero & Overlaid Highs)
# ==============================================================================

def plot_feature_coolwarm(
    umap_coords: np.ndarray,
    values: np.ndarray,
    feature_name: str,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    point_size: float = 16,
    alpha: float = 0.88,
    vmax_percentile: float = 99.0,
    vmin_percentile: float = 1.0,
    color_limits=None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Publication-grade single-cell feature plot:
    - Zero-expression cells plotted first in light neutral gray (#E6E6E6).
    - Expressing cells sorted by value and plotted on top with 'coolwarm'.
    - Color scale clamped to 99th percentile of positive cells to prevent outlier squashing.
    """
    set_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 5.5))
        created_fig = True
    else:
        fig = ax.figure

    values = np.asarray(values, dtype=float)
    zero_mask = (values <= 0) | ~np.isfinite(values)
    pos_mask = ~zero_mask

    # 1. Draw zero/dropout background cells
    ax.scatter(
        umap_coords[zero_mask, 0],
        umap_coords[zero_mask, 1],
        c='#E6E6E6',
        s=point_size * 0.8,
        alpha=0.55,
        rasterized=True,
        label='Zero'
    )

    # 2. Draw expressing cells on top, sorted ascending so highest values appear on the very top
    if pos_mask.sum() > 0:
        pos_idx = np.where(pos_mask)[0]
        sorted_order = pos_idx[np.argsort(values[pos_idx])]
        vmin_val, vmax_val = color_limits if color_limits is not None else np.percentile(values[pos_mask], [vmin_percentile, vmax_percentile])
        if vmax_val <= vmin_val:
            vmin_val, vmax_val = 0., max(float(vmax_val), 1e-3)

        scat = ax.scatter(
            umap_coords[sorted_order, 0],
            umap_coords[sorted_order, 1],
            c=values[sorted_order],
            cmap=CMAP_DIVERGING,
            vmin=vmin_val,
            vmax=vmax_val,
            s=point_size,
            alpha=alpha,
            rasterized=True
        )
        cbar = fig.colorbar(scat, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label(f'{feature_name} log1p expression\npositive-cell p{vmin_percentile:g}-p{vmax_percentile:g}', fontsize=9.5)
        cbar.ax.tick_params(labelsize=8.5)
        ax._plot_audit = {'kind': 'feature', 'n_cells': len(values),
                          'positive_cells': int(pos_mask.sum()), 'nonfinite': int((~np.isfinite(values)).sum()),
                          'vmin': float(vmin_val), 'vmax': float(vmax_val),
                          'below_scale': int((values[pos_mask] < vmin_val).sum()),
                          'above_scale': int((values[pos_mask] > vmax_val).sum()),
                          'scale_only': True, 'shared_scale': color_limits is not None}

    ax.set_title(title if title else f"{feature_name} Expression", fontsize=12, fontweight='bold', pad=8)
    ax.set_xlabel("UMAP1", fontsize=10.5)
    ax.set_ylabel("UMAP2", fontsize=10.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp_name in ax.spines.values():
        sp_name.set_color('#666666')

    return fig, ax


# ==============================================================================
# 2. Split Violin with Inner Boxplot (Case vs Control across Subtypes)
# ==============================================================================

def plot_split_violin_boxplot(
    obs_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    group_col: str,
    case_label: str,
    control_label: str,
    title: Optional[str] = None,
    y_label: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    palette: Optional[Dict[str, str]] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Split violin plot comparing Case vs Control with inner boxplot and median line.
    """
    set_publication_style()
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(9.5, 5.5))
        created_fig = True
    else:
        fig = ax.figure

    if palette is None:
        palette = {case_label: COLOR_CASE, control_label: COLOR_CONTROL}

    # Match the HADHA reference figure: stable subtype order, explicit group
    # order, clipped KDE tails, and a compact box/median inside each half violin.
    # Seaborn otherwise infers order from the first appearance in the dataframe,
    # which makes otherwise identical runs change their x-axis arrangement.
    plot_df = obs_df[[x_col, y_col, group_col]].copy()
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors='coerce')
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col, group_col])
    panel_counts = plot_df.groupby([x_col, group_col], observed=True)[y_col].transform('size')
    plot_df = plot_df[panel_counts > 100].copy()
    if plot_df.empty:
        ax.text(.5, .5, 'No categories with >100 cells', ha='center', transform=ax.transAxes)
        return fig, ax
    order = list(pd.unique(plot_df[x_col]))
    if hasattr(obs_df[x_col].dtype, 'categories'):
        order = [x for x in obs_df[x_col].dtype.categories if x in set(order)]

    sns.violinplot(
        data=plot_df,
        x=x_col,
        y=y_col,
        hue=group_col,
        order=order,
        hue_order=[case_label, control_label],
        palette=palette,
        split=True,
        inner='box',
        cut=0,
        linewidth=1.2,
        saturation=0.95,
        ax=ax
    )

    ax.set_title(title if title else f"{y_col}: {case_label} vs {control_label}", fontsize=12, fontweight='bold', pad=8)
    ax.set_xlabel("Cell Subtype / Cluster", fontsize=10.5)
    ax.set_ylabel(y_label if y_label else y_col, fontsize=10.5)
    ax.legend(title="Group", loc='upper right', frameon=True, fontsize=9.5)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return fig, ax


# ==============================================================================
# 3. Categorized Marker DotPlot with Vertical Divider Lines
# ==============================================================================

def plot_modular_dotplot(
    adata,
    var_dict: Dict[str, List[str]],
    groupby: str,
    out_prefix: Union[Path, str],
    cmap: str = 'coolwarm',
    standard_scale: str = 'var',
    dot_max: float = 0.85,
    dot_min: float = 0.05
):
    """
    Publication dotplot with categorized gene modules, 'coolwarm' colormap,
    clean background (no grid), and vertical dashed divider lines separating modules.
    """
    import scanpy as sc
    set_publication_style()
    out_prefix = figure_path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    raw_vars = set(adata.raw.var_names) if adata.raw is not None else set(adata.var_names)

    # Filter dictionary to keep only existing genes
    filtered_dict = {}
    dividers = []
    current_gene_count = 0

    for cat_name, glist in var_dict.items():
        valid_genes = [g for g in glist if g in raw_vars]
        if valid_genes:
            filtered_dict[cat_name] = valid_genes
            current_gene_count += len(valid_genes)
            dividers.append(current_gene_count - 0.5)

    if not filtered_dict:
        print("[Warning] No genes found for dotplot!")
        return

    # Pop the last divider (not needed at the end of the plot)
    if dividers:
        dividers.pop()

    dp = sc.pl.dotplot(
        adata,
        var_names=filtered_dict,
        groupby=groupby,
        use_raw=True if adata.raw is not None else False,
        standard_scale=standard_scale,
        cmap=cmap,
        dot_max=dot_max,
        dot_min=dot_min,
        return_fig=True,
        show=False
    )
    dp.style(grid=False, cmap=cmap, dot_edge_lw=0.4, dot_edge_color='#333333')
    dp.make_figure()

    ax_main = dp.ax_dict['mainplot_ax']
    ax_main.grid(False)

    # Add vertical dashed lines separating gene modules
    for div_x in dividers:
        ax_main.axvline(div_x, color='#444444', linestyle='--', linewidth=0.9, alpha=0.75)

    dp.savefig(out_prefix.with_suffix('.png'), dpi=300, bbox_inches='tight')
    dp.savefig(out_prefix.with_suffix('.pdf'), bbox_inches='tight')
    finalize_metadata(dp.fig, out_prefix)
    write_figure_index(out_prefix.parent)
    plt.close(dp.fig)


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
        print("[Warning] No genes found for cascade heatmap!")
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


# ==============================================================================
# 7. Volcano Plot (with Labeled Top Genes)
# ==============================================================================

from scientific_plots import plot_volcano
