"""Reusable QC and expression panels. Plotting only; never changes cell selection."""
from pathlib import Path
import re
import textwrap
from figure_output import figure_path, write_figure_index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from plot_utils import (COLOR_CASE, COLOR_CONTROL, PALETTE_12,
                        plot_feature_coolwarm, save_publication_figure,
                        set_publication_style)


def umap_category_palette(adata, key):
    """Bind colors to categorical labels once, reuse before/after rank changes."""
    from plot_utils import PALETTE_6
    import seaborn as sns
    if not isinstance(adata.obs[key].dtype, pd.CategoricalDtype):
        adata.obs[key] = adata.obs[key].astype('category')
    categories = list(adata.obs[key].cat.categories)
    colors = adata.uns.get(f'{key}_colors')
    if colors is None:
        n = len(categories)
        colors = PALETTE_6[:n] if n <= 6 else (PALETTE_12[:n] if n <= 12 else sns.color_palette('husl', n).as_hex())
        adata.uns[f'{key}_colors'] = colors
    from matplotlib.colors import to_hex
    colors = [to_hex(c) for c in colors]
    adata.uns[f'{key}_colors'] = colors
    if len(colors) != len(categories):
        raise ValueError(f'{key}_colors does not match categorical labels; repair the UMAP palette explicitly')
    return dict(zip(map(str, categories), colors))


def sanitize_uns_colors(adata):
    """Normalize display palettes immediately before writing a checkpoint."""
    from matplotlib.colors import to_hex
    for key, colors in list(adata.uns.items()):
        if key.endswith('_colors'):
            adata.uns[key] = [to_hex(c) for c in colors]


def ranked_violin(df, category, value, ax, palette=None, color=None, seed=42,
                  points=100, alpha=.35, title=None):
    """Full-data KDE/IQR; only jitter points are sampled (<=100/category).

    Descending median, descending mean for ties, then category name. Explicit
    zero-based coordinates keep violin, white box, dots and ticks aligned.
    """
    from scipy.stats import gaussian_kde
    set_publication_style()
    data = df[[category, value]].copy()
    data[value] = pd.to_numeric(data[value], errors='coerce')
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    data[category] = data[category].astype(str)
    n_invalid = len(df) - len(data)
    if palette is None and color is None:
        palette = {name: PALETTE_12[i % len(PALETTE_12)] for i, name in enumerate(sorted(data[category].unique()))}
    summary = data.groupby(category, observed=True)[value].agg(['median', 'mean', 'size'])
    omitted = [{'category': str(name), 'n': int(row['size']), 'reason': 'n <= 100'}
               for name, row in summary[summary['size'] <= 100].iterrows()]
    summary = summary[summary['size'] > 100]
    data = data[data[category].isin(summary.index)]
    summary = summary.sort_index().sort_values(['median', 'mean'], ascending=False, kind='stable')
    order = summary.index.tolist()
    rng = np.random.default_rng(seed)
    audit = []
    for i, name in enumerate(order):
        y = data.loc[data[category] == name, value].to_numpy(dtype=float)
        fill = color or (palette or {}).get(name, PALETTE_12[i % len(PALETTE_12)])
        # Constant/all-zero groups remain visible as box/median and points.
        if len(y) >= 10 and np.ptp(y) > 1e-12:
            grid = np.linspace(y.min(), y.max(), 160)
            density = gaussian_kde(y)(grid)
            width = .39 * density / density.max()
            ax.fill_betweenx(grid, i-width, i+width, facecolor=fill,
                             edgecolor='#444444', linewidth=.7, alpha=.9, zorder=1)
        sampled = rng.choice(y, size=min(points, len(y)), replace=False)
        ax.scatter(i+rng.uniform(-.20, .20, len(sampled)), sampled, s=7,
                   c='#333333', alpha=alpha, linewidths=0, rasterized=True, zorder=2)
        ax.boxplot([y], positions=[i], widths=.16, patch_artist=True,
                   showfliers=False, manage_ticks=False,
                   boxprops={'facecolor': 'white', 'edgecolor': '#222222', 'linewidth': .8, 'zorder': 3},
                   medianprops={'color': '#D62728', 'linewidth': 1.7, 'zorder': 4},
                   whiskerprops={'color': '#222222', 'linewidth': .8, 'zorder': 3},
                   capprops={'color': '#222222', 'linewidth': .8, 'zorder': 3})
        audit.append({'category': name, 'color': fill, 'n': len(y), 'points': len(sampled),
                      'density_drawn': bool(len(y) >= 10 and np.ptp(y) > 1e-12),
                      'median': float(np.median(y)), 'mean': float(np.mean(y))})
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(['\n'.join(textwrap.wrap(x, 16, break_long_words=False)) + f'\nn={int(summary.loc[x, "size"]):,}' for x in order],
                       rotation=30 if len(order) > 12 else 0,
                       ha='right' if len(order) > 12 else 'center')
    ax.set_xlim(-.6, max(len(order)-.4, .6))
    ax.set_xlabel(category.replace('_', ' '))
    pretty = {'n_genes_by_counts': 'Detected genes', 'n_genes': 'Detected genes', 'total_counts': 'Total UMI counts', 'pct_counts_mt': 'Mitochondrial counts (%)', 'pct_mt': 'Mitochondrial counts (%)', 'pct_hb': 'Hemoglobin counts (%)', 'pct_ribo': 'Ribosomal counts (%)', 'log10_total_counts': 'log10(UMI counts + 1)', 'doublet_score': 'Doublet score'}
    ax.set_ylabel(pretty.get(value, value.replace('_', ' ')))
    ax.set_title(pretty.get(title or value, title or value))
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    if not len(data):
        ax.text(.5, .5, 'No categories with >100 cells', ha='center', transform=ax.transAxes)
    ax._plot_audit = {'kind': 'ranked_violin', 'seed': seed, 'scatter_alpha': alpha,
                      'n_input': len(df), 'n_nonfinite_or_missing': n_invalid,
                      'minimum_cells_exclusive': 100, 'omitted_categories': omitted,
                      'n_omitted_by_size': sum(row['n'] for row in omitted),
                      'density_and_box': 'all finite cells; box=IQR, whiskers=1.5 IQR',
                      'ordering': 'median descending; mean descending for ties', 'groups': audit}
    return order


def qc_metric_panel(obs, key, destination, seed=42, palette=None):
    metrics = [m for m in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score'] if m in obs]
    if not metrics:
        return
    fig, axes = plt.subplots(len(metrics), 1, figsize=(max(10, 1.4*obs[key].nunique()), 3.7*len(metrics)), squeeze=False)
    for ax, metric in zip(axes.ravel(), metrics):
        ranked_violin(obs, key, metric, ax, palette=palette, seed=seed,
                      title=metric)
    save_publication_figure(fig, destination)


def expression_figure_suite(adata, key, directory, config):
    """All-cell feature/violin, case-over-control violin, separate condition pairs."""
    directory = Path(directory)
    gene = config['TARGET_GENE']
    group = config.get('GROUP_COLUMN', 'group')
    seed = config.get('RANDOM_SEED', 42)
    source = adata.raw if adata.raw is not None else adata
    if gene not in source.var_names:
        return
    x = source[:, gene].X
    values = x.toarray().ravel() if hasattr(x, 'toarray') else np.asarray(x).ravel()
    obs = adata.obs.copy()
    value = f'{gene} log1p expression'
    obs[value] = values
    categories = sorted(obs[key].astype(str).unique())
    palette = umap_category_palette(adata, key)
    positive = values[np.isfinite(values) & (values > 0)]
    pct = config.get('FEATURE_PERCENTILES', [1, 99])
    if not (len(pct) == 2 and 0 <= pct[0] < pct[1] <= 100):
        raise ValueError('FEATURE_PERCENTILES must be [lower, upper] within 0..100')
    limits = tuple(np.percentile(positive, pct)) if positive.size else (0., 1.)
    width = max(10, 1.4*len(categories))
    fig, ax = plot_feature_coolwarm(adata.obsm['X_umap'], values, gene,
                                  color_limits=limits, vmin_percentile=pct[0], vmax_percentile=pct[1])
    save_publication_figure(fig, directory / f'target_{gene}_umap')
    fig, ax = plt.subplots(figsize=(width, 5.2))
    ranked_violin(obs, key, value, ax, palette=palette, seed=seed, title=f'{gene} | all cells')
    save_publication_figure(fig, directory / f'target_{gene}_violin')
    if group not in obs:
        return
    conditions = [(config.get('CASE_LABEL', 'Case'), COLOR_CASE, 'case'),
                  (config.get('CONTROL_LABEL', 'Control'), COLOR_CONTROL, 'control')]
    fig, axes = plt.subplots(2, 1, figsize=(width, 10.5), sharey=True)
    for ax, (condition, color, role) in zip(axes, conditions):
        frame = obs[obs[group].astype(str) == condition]
        ranked_violin(frame, key, value, ax, palette=palette, seed=seed,
                      title=f'{condition} ({role}, n={len(frame):,} cells)')
    save_publication_figure(fig, directory / f'target_{gene}_by_group_violin')
    # Retain shared embedding coordinates, axis bounds, expression scale and y limits.
    xy = np.asarray(adata.obsm['X_umap'])
    lower, upper = np.nanmin(xy, axis=0), np.nanmax(xy, axis=0)
    pad = np.maximum((upper-lower)*.04, .01)
    ymax = max(float(np.nanmax(values))*1.06, .1)
    for condition, color, role in conditions:
        mask = obs[group].astype(str).eq(condition).to_numpy()
        fig, axes = plt.subplots(1, 2, figsize=(width+6, 5.7), gridspec_kw={'width_ratios': [6, width]})
        plot_feature_coolwarm(xy[mask], values[mask], gene, ax=axes[0],
                             title=f'{condition} | {gene} (n={mask.sum():,})',
                             color_limits=limits, vmin_percentile=pct[0], vmax_percentile=pct[1])
        axes[0].set_xlim(lower[0]-pad[0], upper[0]+pad[0])
        axes[0].set_ylim(lower[1]-pad[1], upper[1]+pad[1])
        ranked_violin(obs.loc[mask], key, value, axes[1], palette=palette, seed=seed,
                      title=f'{condition} | {gene} by {key.replace("_", " ")}')
        axes[1].set_ylim(-.03*ymax, ymax)
        safe = re.sub(r'[^\w.-]+', '_', condition)
        save_publication_figure(fig, directory / f'target_{gene}_{role}_{safe}_feature_violin')


def cell_count_bars(obs, destination, case, control, title='Clean cell counts per sample', sample_order=None):
    """One sample, one centered bar. No hue dodge or unobserved zero combinations."""
    counts = obs.groupby('sample_id', observed=True).size()
    order = list(sample_order) if sample_order is not None else sorted(counts.index.astype(str))
    counts = counts.reindex(order, fill_value=0)
    colors = {case: COLOR_CASE, control: COLOR_CONTROL}
    sample_group = obs.groupby('sample_id', observed=True)['group'].agg(lambda x: str(x.iloc[0]))
    if (obs.groupby('sample_id', observed=True)['group'].nunique() > 1).any():
        raise ValueError('Each sample_id must have exactly one condition for cell counts')
    fig, ax = plt.subplots(figsize=(max(8, len(order)*.95), 5))
    bars = ax.bar(np.arange(len(order)), counts, width=.65,
                  color=[colors.get(sample_group.get(s), '#999999') for s in order])
    rotation = 0 if len(order) <= 8 else 30
    ax.set_xticks(np.arange(len(order)), order, rotation=rotation,
                  ha='center' if rotation == 0 else 'right', rotation_mode='anchor')
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.margins(y=.15)
    ax.set_ylabel('Cells')
    ax.set_title(title)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[Patch(color=c, label=g) for g, c in colors.items() if g in set(sample_group)], frameon=False)
    save_publication_figure(fig, destination)
    counts.rename('n_cells').to_csv(figure_path(destination).with_suffix('.csv'))
    write_figure_index(Path(destination).parent)


def qc_scatter(obs, destination, title, color_col='pct_mt', color_limits=None, xy_limits=None):
    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    genes = 'n_genes' if 'n_genes' in obs else 'n_genes_by_counts'
    vmin, vmax = color_limits if color_limits is not None else (None, None)
    scatter = ax.scatter(obs['total_counts'], obs[genes], c=obs[color_col],
                         cmap='viridis', vmin=vmin, vmax=vmax, s=3, alpha=.55, linewidths=0, rasterized=True)
    fig.colorbar(scatter, ax=ax, label=color_col)
    ax.set(xlabel='Total counts', ylabel='Detected genes', title=f'{title} (n={len(obs):,})')
    if (obs['total_counts'] > 0).all():
        ax.set_xscale('log')
        ax.set_xlabel('Total counts (log scale)')
    if xy_limits is not None:
        left = xy_limits[2] if len(xy_limits) > 2 else max(float(obs['total_counts'].min())*.9, .1)
        ax.set_xlim(left if ax.get_xscale() == 'log' else 0, xy_limits[0])
        ax.set_ylim(0, xy_limits[1])
    ax.spines[['top', 'right']].set_visible(False)
    save_publication_figure(fig, destination)


def scrublet_histograms(distributions, stats, destination):
    """Observed/simulated densities and both thresholds, in a compact grid."""
    n = len(stats)
    if not n:
        return
    columns = min(3, n)
    rows = int(np.ceil(n/columns))
    fig, axes = plt.subplots(rows, columns, figsize=(columns*4.5, rows*3.6), squeeze=False)
    for ax in axes.ravel():
        ax.set_visible(False)
    for ax, row in zip(axes.ravel(), stats):
        ax.set_visible(True)
        observed, simulated = distributions[row['sample_id']]
        used, auto = row['effective_threshold'], row['auto_threshold']
        bins = np.linspace(0, max(np.max(observed), np.max(simulated), used)*1.05, 41)
        ax.hist(observed, bins=bins, density=True, color=COLOR_CONTROL, alpha=.7, label='Observed')
        ax.hist(simulated, bins=bins, density=True, color=COLOR_CASE, alpha=.4, label='Simulated doublets')
        if auto is not None and np.isfinite(auto):
            ax.axvline(auto, color='#777777', linestyle=':', label=f'Auto {auto:.3f}')
        ax.axvline(used, color='#222222', linestyle='--', label=f'Used {used:.3f}')
        ax.set(title=f"{row['sample_id']}: {row['doublets']:,} doublets", xlabel='Doublet score', ylabel='Density')
        ax.legend(fontsize=8, frameon=False)
        ax.spines[['top', 'right']].set_visible(False)
    save_publication_figure(fig, destination)


def scrublet_comparison(before, after, destination):
    """Original two-panel post-QC vs final comparison colored by doublet score."""
    limits = np.nanpercentile(before['doublet_score'], [1, 99])
    if limits[0] == limits[1]:
        limits = (0, max(float(limits[1]), .01))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True, sharey=True)
    for ax, (obs, label) in zip(axes, [(before, 'After QC'), (after, 'Final clean')]):
        genes = 'n_genes' if 'n_genes' in obs else 'n_genes_by_counts'
        sca = ax.scatter(obs.total_counts, obs[genes], c=obs.doublet_score,
                         cmap='viridis', vmin=limits[0], vmax=limits[1], s=3,
                         alpha=.6, linewidths=0, rasterized=True)
        ax.set(title=f'{label} (n={len(obs):,})', xlabel='Total counts', ylabel='Detected genes')
        ax.spines[['top', 'right']].set_visible(False)
        fig.colorbar(sca, ax=ax, label='Doublet score (shared p1-p99)')
    save_publication_figure(fig, destination)
