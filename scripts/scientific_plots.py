"""Plot-only scientific display helpers; preserve the statistical input tables."""
import textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from figure_output import figure_path


def plot_volcano(deg_df, gene_col, log2fc_col, pval_col, out_prefix,
                 title='Differential expression', fc_cutoff=.25, p_cutoff=.05,
                 top_n_labels=10, exclude_genes=None):
    from plot_utils import set_publication_style, save_publication_figure
    set_publication_style()
    prefix = figure_path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    data = deg_df.copy()
    excluded = [exclude_genes] if isinstance(exclude_genes, str) else list(exclude_genes or [])
    data = data[~data[gene_col].astype(str).str.upper().isin({x.upper() for x in excluded})].copy()
    x = pd.to_numeric(data[log2fc_col], errors='coerce')
    p = pd.to_numeric(data[pval_col], errors='coerce')
    valid = np.isfinite(x) & np.isfinite(p) & p.between(0,1)
    data = data.loc[valid].copy()
    x, p = x.loc[valid].to_numpy(), p.loc[valid].to_numpy()
    if not len(data):
        raise ValueError('No valid genes remain for the volcano plot')
    y = -np.log10(np.clip(p, 1e-300, 1))
    up = (x >= fc_cutoff) & (p < p_cutoff)
    down = (x <= -fc_cutoff) & (p < p_cutoff)
    # Bound the display only. Boundary triangles explicitly mark saturated points.
    informative = up | down
    # Nonsignificant near-zero-expression genes can have enormous unstable FC.
    # Use the significant cloud for display bounds when it exists, retaining
    # every other gene as a regular point or an explicit saturated triangle.
    basis = informative if informative.any() else np.ones(len(x), dtype=bool)
    xlim = max(fc_cutoff*1.5, min(6., max(1.5, float(np.quantile(np.abs(x[basis]), .99)))))
    ylim = max(-np.log10(p_cutoff)*1.2, min(80., max(5., float(np.quantile(y[basis], .99)))))
    px, py = np.clip(x,-xlim,xlim), np.minimum(y,ylim)
    clipped = (px != x) | (py != y)
    fig, ax = plt.subplots(figsize=(10,6.5))
    for mask, color, label in [(~(up|down),'#bcbcbc','Not significant'),(up,'#d62728','Up'),(down,'#2166ac','Down')]:
        ax.scatter(px[mask & ~clipped], py[mask & ~clipped], s=12, c=color,
                   alpha=.65, linewidths=0, rasterized=True, label=f'{label} ({mask.sum():,})')
        ax.scatter(px[mask & clipped], py[mask & clipped], s=19, c=color,
                   alpha=.8, marker='^', linewidths=.25, edgecolors='#444444', rasterized=True)
    ax.axvline(fc_cutoff, color='#555555', lw=.8, ls='--')
    ax.axvline(-fc_cutoff, color='#555555', lw=.8, ls='--')
    ax.axhline(-np.log10(p_cutoff), color='#555555', lw=.8, ls='--')
    for sign, mask in [(-1,down),(1,up)]:
        ids = np.flatnonzero(mask)
        ids = ids[np.argsort(-y[ids], kind='stable')[:max(0,top_n_labels//2)]]
        # Tied saturated points fan out in x order, avoiding crossed leaders.
        ids = ids[np.lexsort((sign*px[ids], py[ids]))]
        center = float(np.median(py[ids])) if len(ids) else 0.
        span = min(ylim*.38, max(0, len(ids)-1)*ylim*.085)
        low = np.clip(center-span/2, ylim*.08, ylim*.97-span)
        for idx, label_y in zip(ids, np.linspace(low, low+span, len(ids))):
            ax.annotate(str(data.iloc[idx][gene_col]), (px[idx],py[idx]),
                        xytext=(sign*xlim*1.12,label_y), fontsize=9,
                        ha='left' if sign > 0 else 'right', va='center',
                        bbox={'facecolor':'white', 'edgecolor':'none', 'pad':1.5},
                        arrowprops={'arrowstyle':'-', 'color':'#777777','lw':.6})
    ax.set(xlim=(-xlim*1.65,xlim*1.65), ylim=(-ylim*.025,ylim*1.06),
           xlabel='log2 fold change', ylabel='−log10 adjusted P value', title=title)
    ax.legend(loc='upper center', bbox_to_anchor=(.5,1.18), ncol=3, frameon=False)
    ax.text(.5,-.16, f'Triangles: {clipped.sum():,} genes at display limits; full values retained in tables.',
            ha='center', transform=ax.transAxes, fontsize=8.5)
    ax.spines[['top','right']].set_visible(False)
    data['display_x'],data['display_y'],data['display_clipped'] = px,py,clipped
    data.to_csv(prefix.with_suffix('.plot_data.csv'),index=False)
    ax._plot_audit = {'kind':'volcano','n_input':len(deg_df),'n_plotted':len(data),
                      'excluded_target_genes':excluded,'invalid_values':int((~valid).sum()),
                      'x_display_limit':xlim,'y_display_limit':ylim,'n_clipped':int(clipped.sum()),
                      'display_limit_basis':'significant-gene p99, bounded' if informative.any() else 'all-gene p99, bounded',
                      'fdr_cutoff':p_cutoff,'abs_log2fc_cutoff':fc_cutoff,
                      'n_up':int(up.sum()),'n_down':int(down.sum()), 'tables_modified':False}
    save_publication_figure(fig,prefix)


def plot_ora_results(results, destination, top_n=12):
    from plot_utils import save_publication_figure
    data=results.copy()
    data['Adjusted P-value']=pd.to_numeric(data['Adjusted P-value'],errors='coerce')
    data=data[np.isfinite(data['Adjusted P-value']) & data['Adjusted P-value'].between(0,1)]
    data=data.sort_values('Adjusted P-value',kind='stable').head(top_n).iloc[::-1]
    if data.empty:
        return
    labels=['\n'.join(textwrap.wrap(str(t),48,break_long_words=False)) for t in data['Term']]
    fig,ax=plt.subplots(figsize=(10,max(5,len(data)*.55)))
    ax.barh(np.arange(len(data)), -np.log10(np.clip(data['Adjusted P-value'],1e-300,1)), color='#0072B2')
    ax.set_yticks(np.arange(len(data)),labels,fontsize=9)
    ax.set(xlabel='−log10 adjusted P value',title='Over-representation analysis')
    ax.spines[['top','right']].set_visible(False)
    save_publication_figure(fig,destination)


def plot_gsea_results(results, destination):
    from plot_utils import save_publication_figure
    data=results.copy()
    for col in ['NES','FDR q-val']:
        data[col]=pd.to_numeric(data[col],errors='coerce')
    data=data[np.isfinite(data.NES)&np.isfinite(data['FDR q-val'])].sort_values('FDR q-val')
    data=pd.concat([data[data.NES>0].head(8),data[data.NES<0].head(8)]).sort_values('NES')
    if data.empty:
        return
    fig,ax=plt.subplots(figsize=(11,max(5,len(data)*.6)))
    bars=ax.barh(np.arange(len(data)),data.NES,color=['#d62728' if x>0 else '#2166ac' for x in data.NES])
    for bar,q in zip(bars,data['FDR q-val']):
        if q>=.05:
            bar.set_alpha(.3)
            bar.set_hatch('//')
    labels=['\n'.join(textwrap.wrap(str(t),45,break_long_words=False)) for t in data.Term]
    ax.set_yticks(np.arange(len(data)),labels,fontsize=9)
    ax.axvline(0,color='#555555',lw=.8)
    ax.set(xlabel='Normalized enrichment score (NES)',title='GSEA pathway scores')
    ax.text(.5,-.12,'Hatched bars: FDR ≥ 0.05. Positive/negative NES is enrichment direction.',ha='center',transform=ax.transAxes,fontsize=8.5)
    ax.spines[['top','right']].set_visible(False)
    save_publication_figure(fig,destination)
