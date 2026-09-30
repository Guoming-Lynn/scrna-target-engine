"""Stable figure names, lightweight layout checks and a browsable output index."""
from pathlib import Path
import csv
import html
import json
import re

STAGES = {'qc': '01_qc', 'global': '02_global', 'target': '03_target', 'deg': '04_deg'}

NAMES = {
    '25_UMAP_legend': 'umap_cell_labels',
    '25_UMAP_on_data': 'umap_cell_labels_on_map',
    '25_target_celltype_annotated_UMAP': 'umap_cell_labels',
    '18_cluster_QC': 'qc_metrics_by_cluster',
    'qc_metrics_violin': 'qc_metrics_by_cluster',
    '27_group_composition': 'composition_cell_types_by_condition',
    '27_tissue_group_composition_stacked_barplot': 'composition_cell_types_by_condition',
    '26_subtype_composition_case_vs_control': 'composition_conditions_by_cell_type',
    '16_marker_dotplot': 'markers_top_dotplot',
    '16_canonical_state_marker_dotplot': 'markers_modules_dotplot',
    'canonical_markers_dotplot': 'markers_canonical_dotplot',
    'umap_overview': 'umap_metadata_and_qc',
    'umap_leiden_multires': 'umap_resolution_grid',
    'final_clean_cell_counts': 'cell_counts_clean',
    'raw_cell_counts': 'cell_counts_raw',
    'pre_QC_scatter': 'qc_pre_counts_vs_genes',
    'post_QC_scatter': 'qc_post_counts_vs_genes',
    'scrublet_histograms_grid': 'scrublet_score_distributions',
    'QC_vs_final_scatter_genes_vs_UMI_colored_scrublet_score': 'scrublet_before_after_counts_vs_genes',
    '04_target_high_low_subtype_composition_stacked_barplot': 'composition_target_expression_strata',
    '08_Enrichment_Case_UP_barplot': 'ora_case_up_pathways',
    '08b_Enrichment_Case_DOWN_barplot': 'ora_case_down_pathways',
    '10_GSEA_top_pathways_NES_barplot': 'gsea_case_pathway_nes',
}


def figure_path(prefix):
    prefix = Path(prefix)
    name = prefix.name
    if name in NAMES:
        name = NAMES[name]
    elif name.startswith('pre_QC_') or name.startswith('post_QC_'):
        name = name.replace('pre_QC_', 'qc_pre_').replace('post_QC_', 'qc_post_')
    elif name.startswith('14_'):
        name = 'target_' + name[3:].replace('_feature_and_violin', '_feature_violin')
    elif name.startswith('violin_'):
        name = 'target_' + name[7:] + '_violin'
    elif re.match(r'^umap_[A-Z0-9]', name):
        name = 'target_' + name[5:] + '_umap'
    elif name.startswith('05_Volcano_'):
        name = 'volcano_' + name[len('05_Volcano_'):]
    else:
        name = re.sub(r'^\d+[a-z]?_', '', name)
    name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', name)
    return prefix.with_name(name)


def audit_layout(fig):
    """Check visible neighboring tick labels after rendering; no inference QA."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    problems = []
    for i, ax in enumerate(fig.axes):
        if not ax.get_visible():
            continue
        for direction, ticks in [('x', ax.get_xticklabels()), ('y', ax.get_yticklabels())]:
            boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in ticks if t.get_visible() and t.get_text()]
            for (first, a), (second, b) in zip(boxes, boxes[1:]):
                if a.overlaps(b):
                    problems.append({'axis': i, 'direction': direction, 'labels': [first, second]})
    return problems


def finalize_metadata(fig, prefix):
    prefix = Path(prefix)
    problems = audit_layout(fig)
    audits = [ax._plot_audit for ax in fig.axes if hasattr(ax, '_plot_audit')]
    captions = []
    for panel in audits:
        if panel.get('kind') == 'ranked_violin':
            captions.append('n denotes cells. Violin density and box summaries use all finite cells; '
                            'white box: 25th–75th percentiles; red line: median; whiskers: 1.5 IQR. '
                            'Only categories with >100 finite cells in this panel are displayed; data and analysis are unchanged. Only jitter points are sampled, at most 100/category. Constant groups omit density. No inferential test is shown.')
        elif panel.get('kind') == 'feature':
            captions.append('Color displays normalized log1p expression; zero-expression cells are gray. '
                            'Positive-cell percentile limits saturate color only; all cells and expression values are retained.')
        elif panel.get('kind') == 'volcano':
            captions.append(f"Red/blue indicate positive/negative fold change with adjusted P < {panel['fdr_cutoff']} and |log2FC| ≥ {panel['abs_log2fc_cutoff']}. "
                            f"The configured target is excluded from the display. {panel['n_clipped']} boundary triangles mark saturated display positions; full values are in plot_data.csv. "
                            'This is an exploratory cell-level contrast, not a test of donor-level disease effects.')
    if captions:
        prefix.with_suffix('.caption.txt').write_text('\n'.join(dict.fromkeys(captions)), encoding='utf-8')
    prefix.with_suffix('.layout.json').write_text(json.dumps({'tick_overlap_warnings': problems}, ensure_ascii=False, indent=2), encoding='utf-8')


def write_figure_index(directory):
    """Index actual outputs only, so optional/skipped panels are never promised."""
    directory = Path(directory)
    rows = []
    for png in sorted(directory.rglob('*.png')):
        def relative(suffix):
            candidate = png.with_suffix(suffix)
            return candidate.relative_to(directory).as_posix() if candidate.exists() else ''
        rows.append({'figure': png.relative_to(directory).with_suffix('').as_posix(),
                     'png': relative('.png'), 'pdf': relative('.pdf'),
                     'caption': relative('.caption.txt'), 'audit': relative('.plot.json'),
                     'data': relative('.plot_data.csv') or relative('.csv')})
    with (directory/'figure_index.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=['figure','png','pdf','caption','audit','data'])
        writer.writeheader()
        writer.writerows(rows)
    cards = []
    for row in rows:
        links = ' · '.join(f'<a href="{html.escape(row[k], quote=True)}">{k.upper()}</a>' for k in ['png','pdf','caption','audit','data'] if row[k])
        cards.append(f'<article><h2>{html.escape(row["figure"])}</h2>{links}<a href="{html.escape(row["png"], quote=True)}"><img loading="lazy" src="{html.escape(row["png"], quote=True)}"></a></article>')
    folders = '<nav>' + ' · '.join(f'<a href="{html.escape(p.name,quote=True)}/index.html">{html.escape(p.name)}</a>' for p in sorted(directory.iterdir()) if p.is_dir() and (p/'index.html').exists()) + '</nav>'
    (directory/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Figure index</title>'
        '<style>body{font:16px sans-serif;max-width:1200px;margin:32px auto;background:#fafafa}article{background:white;padding:24px;margin-bottom:24px}h2{font-size:18px}img{display:block;max-width:100%;max-height:700px;margin:20px auto}a{color:#0072b2}</style>'
        '<h1>Figure index</h1>'+folders+''.join(cards), encoding='utf-8')
    if directory.name in {'01_qc','02_global','03_target','04_deg','scrublet'}:
        write_figure_index(directory.parent)
