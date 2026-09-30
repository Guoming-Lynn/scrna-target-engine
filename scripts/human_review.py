"""Checkpoint-bound human decisions and count-based reclustering rounds.

Inspired by HADHA_NEW run_step02_round2.py, run_step03_round3.py and
update_all_requested_figures.py. Never import these original execution scripts.
"""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def validate(decision, source_hash, clusters):
    if decision.get('source_sha256') != source_hash:
        raise ValueError('Review belongs to a different checkpoint; use the current round template')
    if decision.get('approved') is not True:
        return False
    rows = decision.get('clusters', [])
    if len(rows) != len(clusters) or {str(r['cluster']) for r in rows} != set(clusters):
        raise ValueError('Review must contain exactly the current clusters')
    if any(r.get('action') not in ('keep', 'drop') for r in rows):
        raise ValueError('Human must select keep or drop for every cluster')
    kept = [r for r in rows if r['action'] == 'keep']
    if not kept:
        raise ValueError('Cannot drop all clusters')
    if not any(r['action'] == 'drop' for r in rows):
        if any(not str(r.get('final_name', '')).strip() for r in kept):
            raise ValueError('Final names are required when accepting a round')
    return True


def draw_review(adata, key, directory, config):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import scanpy as sc
    from plot_utils import (save_publication_figure, plot_modular_dotplot,
                            plot_two_bar_composition, plot_feature_coolwarm,
                            plot_split_violin_boxplot, PALETTE_6, PALETTE_12)
    import seaborn as sns
    from figure_recipes import qc_metric_panel, expression_figure_suite, umap_category_palette
    directory.mkdir(parents=True, exist_ok=True)
    count = adata.obs[key].nunique()
    colors = list(umap_category_palette(adata, key).values())
    for location, suffix in [('right margin', 'legend'), ('on data', 'on_data')]:
        fig, ax = plt.subplots(figsize=(8.8, 6.8))
        sc.pl.umap(adata, color=key, palette=colors, legend_loc=location, size=20, alpha=.88, frameon=True, ax=ax, show=False)
        ax.set_title(f'{key} (n={adata.n_obs:,})', fontsize=13.5, fontweight='bold')
        ax.grid(False)
        save_publication_figure(fig, directory / ('umap_cell_labels' if suffix == 'legend' else 'umap_cell_labels_on_map'))
    metrics = [c for c in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score', 'scrublet_score'] if c in adata.obs]
    if metrics:
        qc_metric_panel(adata.obs, key, directory / 'qc_metrics_by_cluster', config.get('RANDOM_SEED', 42), dict(zip(adata.obs[key].cat.categories.astype(str), colors)))
        adata.obs.groupby(key, observed=True)[metrics].median().to_csv(directory / 'cluster_QC_medians.csv')
    group = config.get('GROUP_COLUMN', 'group')
    if group in adata.obs:
        pd.crosstab(adata.obs[key], adata.obs[group]).to_csv(directory / 'cluster_group_counts.csv')
        plot_two_bar_composition(adata.obs, group, key, colors, directory / 'composition_cell_types_by_condition', config['CASE_LABEL'], config['CONTROL_LABEL'])

    expression_figure_suite(adata, key, directory, config)
    if count > 1 and adata.obs[key].value_counts().min() > 1:
        sc.tl.rank_genes_groups(adata, groupby=key, method='wilcoxon', use_raw=True)
        markers = sc.get.rank_genes_groups_df(adata, group=None)
        markers.to_csv(directory / 'markers_all.csv', index=False)
        top = markers[markers.logfoldchanges > 0].groupby('group', observed=True).head(20)
        top.to_csv(directory / 'markers_top20.csv', index=False)
        modules = {str(c): g.names.head(2).tolist() for c, g in top.groupby('group', observed=True)}
        if modules:
            plot_modular_dotplot(adata, modules, key, directory / 'markers_top_dotplot')


def review(checkpoint, phase, config):
    import shutil
    import pandas as pd
    import scanpy as sc
    checkpoint = Path(checkpoint)
    # QUICK_MODE is intentionally permissive: publish provisional cluster
    # labels and continue without waiting for a decision file.  This mode is
    # for rapid exploration; the regular checkpoint-bound review remains the
    # default for audited runs.
    if config.get('QUICK_MODE', False):
        import pandas as pd
        import scanpy as sc
        a = sc.read_h5ad(checkpoint)
        key = 'leiden' if phase == 'phase02' else 'target_leiden_selected'
        label = 'cell_type' if phase == 'phase02' else 'subtype'
        vals = a.obs[key].astype(str)
        names = {c: f'Cluster_{c}' for c in sorted(vals.unique())}
        a.obs[label] = pd.Categorical(vals.map(names), categories=list(names.values()), ordered=True)
        a.uns['annotation_status'] = 'quick_provisional'
        a.uns['quick_mode'] = True
        from figure_recipes import sanitize_uns_colors
        sanitize_uns_colors(a)
        stage = '02_global' if phase == 'phase02' else '03_target'
        stage_dir = Path(config['OUTPUT_ROOT']) / '04_figures' / stage
        figure_status = 'generated'
        figure_error = None
        try:
            # Quick mode publishes the same label and marker panels as a
            # reviewed run, while retaining provisional names.
            draw_review(a, label, stage_dir, config)
            from figure_output import write_figure_index
            write_figure_index(stage_dir)
        except Exception as exc:
            figure_status = 'failed_optional'
            figure_error = repr(exc)
        a.write_h5ad(checkpoint)
        table_dir = Path(config['OUTPUT_ROOT']) / '03_tables' / ('02_global' if phase == 'phase02' else '03_target')
        table_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({'cluster': list(names), 'final_name': list(names.values()), 'status': 'provisional'}).to_csv(
            table_dir / ('cell_type_naming_quick.csv' if phase == 'phase02' else 'subtype_naming_quick.csv'), index=False)
        report = {
            'phase': phase,
            'annotation_status': 'quick_provisional',
            'analysis_unit': 'cell',
            'donor_level_inference': False,
            'n_cells': int(a.n_obs),
            'n_clusters': int(len(names)),
            'mapping': names,
            'figure_status': figure_status,
        }
        if phase == 'phase03' and isinstance(a.uns.get('quick_target_selection'), dict):
            report['quick_target_selection'] = a.uns['quick_target_selection']
        if figure_error:
            report['figure_error'] = figure_error
        report_file = '02_global_review_report.json' if phase == 'phase02' else '03_target_report.json'
        def _jsonable(obj):
            if hasattr(obj, 'tolist') and not isinstance(obj, (str, bytes)):
                return obj.tolist()
            if hasattr(obj, 'item') and not isinstance(obj, (str, bytes)):
                return obj.item()
            raise TypeError(f'Object of type {type(obj).__name__} is not JSON serializable')
        (Path(config['OUTPUT_ROOT']) / '99_logs' / report_file).write_text(
            json.dumps(report, indent=2, default=_jsonable), encoding='utf-8')
        return 0
    root = Path(config['OUTPUT_ROOT']) / '00_metadata' / 'human_review' / phase
    root.mkdir(parents=True, exist_ok=True)
    source_hash = digest(checkpoint)
    current = root / 'current.json'
    state = json.loads(current.read_text()) if current.exists() else {}
    if state.get('accepted_sha256') == source_hash:
        return 0
    key = 'leiden' if phase == 'phase02' else 'target_leiden_selected'
    label = 'cell_type' if phase == 'phase02' else 'subtype'
    a = sc.read_h5ad(checkpoint)
    clusters = sorted(a.obs[key].astype(str).unique())
    if state.get('source_sha256') != source_hash:
        number = max([int(p.name.split('_')[1]) for p in root.glob('round_*') if p.is_dir()] or [0]) + 1
        rd = root / f'round_{number:02d}'
        rd.mkdir()
        shutil.copy2(checkpoint, rd / 'source.h5ad')
        template = {'source_sha256': source_hash, 'approved': False, 'clusters': [dict(cluster=c, action='', final_name='', reason='') for c in clusters]}
        (rd / 'decision.json').write_text(json.dumps(template, indent=2, ensure_ascii=False), encoding='utf-8')
        state = {'source_sha256': source_hash, 'round_dir': str(rd)}
        current.write_text(json.dumps(state, indent=2), encoding='utf-8')
        draw_review(a, key, rd / 'figures_and_tables', config)
        print(f'Needs human review: {rd / "decision.json"}', flush=True)
        print('Use final_name for broad cell types (TARGET_CELLTYPE_NAME).' if phase == 'phase02' else 'Select keep/drop; supply final subtype names when accepting.', flush=True)
        return 2
    rd = Path(state['round_dir'])
    decision = json.loads((rd / 'decision.json').read_text(encoding='utf-8-sig'))
    if not validate(decision, source_hash, clusters):
        print(f'Waiting for human review: {rd / "decision.json"}', flush=True)
        return 2
    drops = [str(r['cluster']) for r in decision['clusters'] if r['action'] == 'drop']
    if drops:
        mask = a.obs[key].astype(str).isin(drops)
        if (~mask).sum() < 20:
            raise ValueError('Fewer than 20 cells would remain')
        a.obs.loc[mask].to_csv(rd / 'removed_cells.csv')
        a = a[~mask].copy()
        if 'counts' not in a.layers:
            raise ValueError('Reclustering requires full-gene raw counts')
        a.obsm.clear()
        a.obsp.clear()
        a.uns.clear()
        for col in list(a.obs):
            if 'leiden' in col or col in ('cell_type', 'subtype') or col.startswith('score_'):
                del a.obs[col]
        from phase03_subclustering import reinitialize_feature_space, run_subclustering_dr
        settings = dict(config)
        if phase == 'phase02':
            settings['N_HVG_SUBCLUSTERING'] = config.get('N_HVG_GLOBAL', 3500)
            settings['N_PCS_SUBCLUSTERING'] = config.get('N_PCS_GLOBAL', 50)
        a = run_subclustering_dr(reinitialize_feature_space(a), settings)
        if phase == 'phase02':
            a.obs['leiden'] = a.obs['target_leiden_selected'].copy()
        a.obs[label] = pd.Categorical('Cluster_' + a.obs[key].astype(str))
        from figure_recipes import sanitize_uns_colors
        sanitize_uns_colors(a)
        a.write_h5ad(checkpoint)
        # New cluster IDs always require a fresh human decision.
        return review(checkpoint, phase, config)
    mapping = {str(r['cluster']): str(r['final_name']).strip() for r in decision['clusters']}
    a.obs[label] = pd.Categorical(a.obs[key].astype(str).map(mapping))
    # Renaming/merging categories invalidates the old positional color array.
    # Register a palette for final names once; later figure-only runs preserve it.
    a.uns.pop(f'{label}_colors', None)
    a.uns['human_review'] = {'source_sha256': source_hash, 'decision': str(rd / 'decision.json'), 'status': 'accepted'}
    draw_review(a, label, rd / 'final_figures', config)
    # Publish final reviewed panels to the normal stage folder. Preserve all
    # preliminary panels with the review round, so contaminants cannot be
    # mistaken for the final deliverables.
    stage = '02_global' if phase == 'phase02' else '03_target'
    stage_dir = Path(config['OUTPUT_ROOT']) / '04_figures' / stage
    archive = rd / 'pre_review_stage_figures'
    if stage_dir.exists() and not archive.exists():
        shutil.move(str(stage_dir), str(archive))
    shutil.copytree(rd / 'final_figures', stage_dir, dirs_exist_ok=True)
    table_dir = Path(config['OUTPUT_ROOT']) / '03_tables' / stage
    table_dir.mkdir(parents=True, exist_ok=True)
    for csv in (rd / 'final_figures').glob('*.csv'):
        if csv.name != 'figure_index.csv':
            shutil.copy2(csv, table_dir / csv.name)
    from figure_output import write_figure_index
    write_figure_index(stage_dir)
    pd.DataFrame({'cluster': list(mapping), 'final_name': list(mapping.values())}).to_csv(
        table_dir / ('cell_type_naming_final.csv' if phase == 'phase02' else 'subtype_naming_template.csv'), index=False)
    report = {'phase': phase, 'annotation_status': 'human_accepted', 'n_cells': a.n_obs,
              'n_subtypes' if phase == 'phase03' else 'n_cell_types': a.obs[label].nunique(), 'mapping': mapping}
    report_file = '03_target_report.json' if phase == 'phase03' else '02_global_review_report.json'
    (Path(config['OUTPUT_ROOT']) / '99_logs' / report_file).write_text(json.dumps(report, indent=2), encoding='utf-8')
    a.uns['annotation_status'] = 'human_accepted'
    from figure_recipes import sanitize_uns_colors
    sanitize_uns_colors(a)
    a.write_h5ad(checkpoint)
    state['accepted_sha256'] = digest(checkpoint)
    current.write_text(json.dumps(state, indent=2), encoding='utf-8')
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkpoint', required=True)
    parser.add_argument('--phase', choices=['phase02', 'phase03'], required=True)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    raise SystemExit(review(args.checkpoint, args.phase, json.loads(Path(args.config).read_text(encoding='utf-8-sig'))))
