#!/usr/bin/env python3
"""
Phase 2: Global Clustering and Annotation
scRNA-seq Target-Gene Analysis Pipeline
"""

import sys
import json
import logging
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from config_contract import batch_column
from plot_utils import save_publication_figure, COLOR_CASE, COLOR_CONTROL
from figure_recipes import qc_metric_panel, expression_figure_suite, umap_category_palette
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else 'config.json'
    with open(config_path, 'r', encoding='utf-8-sig') as f:
        config = json.load(f)

    # Config keys used
    TARGET_GENE = config.get('TARGET_GENE', 'TARGET')
    PROJECT_NAME = config.get('PROJECT_NAME', 'scRNA_Project')
    OUTPUT_ROOT = Path(config.get('OUTPUT_ROOT', './output'))
    BATCH_COLUMN = batch_column(config)
    GROUP_COLUMN = config.get('GROUP_COLUMN', 'group')
    CASE_LABEL = config.get('CASE_LABEL', 'Case')
    CONTROL_LABEL = config.get('CONTROL_LABEL', 'Control')
    RANDOM_SEED = config.get('RANDOM_SEED', 42)
    N_HVG_GLOBAL = config.get('N_HVG_GLOBAL', 2000)
    N_PCS_GLOBAL = config.get('N_PCS_GLOBAL', 30)
    LEIDEN_RESOLUTIONS = config.get('LEIDEN_RESOLUTIONS', [0.3, 0.5, 0.8, 1.2])
    
    # Standard colors

    
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Define outputs
    h5ad_dir = OUTPUT_ROOT / '05_h5ad'
    in_h5ad = h5ad_dir / '03_clean.h5ad'
    out_h5ad = h5ad_dir / '04_global_clustered.h5ad'
    
    table_dir = OUTPUT_ROOT / '03_tables' / '02_global'
    table_dir.mkdir(parents=True, exist_ok=True)
    
    fig_dir = OUTPUT_ROOT / '04_figures' / '02_global'
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    if not in_h5ad.exists():
        logging.error(f"Input file not found: {in_h5ad}")
        sys.exit(1)
        
    logging.info(f"Loading data from {in_h5ad}")
    adata = sc.read_h5ad(in_h5ad)
    
    sc.settings.seed = RANDOM_SEED
    sc.settings.verbosity = 3
    
    # 5. Normalize target_sum=1e4, log1p. Store full-gene matrix in adata.raw.
    logging.info("Normalizing and log-transforming data")
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
        
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    # 6. Select HVGs batch-aware on BATCH_COLUMN
    logging.info(f"Selecting {N_HVG_GLOBAL} highly variable genes, batch_key='{BATCH_COLUMN}'")
    try:
        sc.pp.highly_variable_genes(
            adata, 
            n_top_genes=N_HVG_GLOBAL, 
            flavor='seurat_v3', 
            batch_key=BATCH_COLUMN,
            layer='counts',
            subset=False
        )
    except Exception as e:
        logging.warning(f"seurat_v3 flavor failed: {e}. Falling back to seurat.")
        sc.pp.highly_variable_genes(
            adata, 
            n_top_genes=N_HVG_GLOBAL, 
            flavor='seurat', 
            batch_key=BATCH_COLUMN,
            subset=False
        )
        
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    adata_hvg.raw = None
    
    # 7. Scale, PCA
    logging.info(f"Scaling and computing PCA with {N_PCS_GLOBAL} components")
    sc.pp.scale(adata_hvg, max_value=10)
    N_PCS_GLOBAL = min(N_PCS_GLOBAL, min(adata_hvg.shape) - 1)
    sc.tl.pca(adata_hvg, n_comps=N_PCS_GLOBAL, svd_solver='arpack', random_state=RANDOM_SEED)
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca'].copy()
    del adata_hvg
    
    # 8. Harmony batch correction
    logging.info(f"Running Harmony batch correction on '{BATCH_COLUMN}'")
    import harmonypy as hm
    if BATCH_COLUMN in adata.obs and adata.obs[BATCH_COLUMN].nunique() > 1:
        ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, BATCH_COLUMN, max_iter_harmony=20, random_state=RANDOM_SEED)
        corrected = np.asarray(ho.Z_corr)
        if corrected.shape[0] != adata.n_obs and corrected.shape[1] == adata.n_obs:
            corrected = corrected.T
        if corrected.shape[0] != adata.n_obs:
            raise ValueError(f'Harmony output shape {corrected.shape} does not match cells {adata.n_obs}')
        adata.obsm['X_pca_harmony'] = corrected
    else:
        adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()
    
    # 9. Build neighbors and UMAP
    logging.info("Computing neighbors and UMAP on Harmony embeddings")
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=N_PCS_GLOBAL, use_rep='X_pca_harmony', metric='cosine', random_state=RANDOM_SEED)
    sc.tl.umap(adata, min_dist=0.35, random_state=RANDOM_SEED)
    
    # 10. Multi-resolution Leiden
    logging.info(f"Running Leiden clustering at resolutions: {LEIDEN_RESOLUTIONS}")
    for res in LEIDEN_RESOLUTIONS:
        res_key = f'leiden_r{str(res).replace(".", "_")}'
        sc.tl.leiden(adata, resolution=res, key_added=res_key, random_state=RANDOM_SEED)
        
    # 11. Select optimal resolution (first giving 10-15 clusters, or default to median/first)
    selected_res = LEIDEN_RESOLUTIONS[0]
    selected_key = f'leiden_r{str(selected_res).replace(".", "_")}'
    for res in LEIDEN_RESOLUTIONS:
        res_key = f'leiden_r{str(res).replace(".", "_")}'
        n_clusters = len(adata.obs[res_key].unique())
        logging.info(f"Resolution {res}: {n_clusters} clusters")
        if 10 <= n_clusters <= 15:
            selected_res = res
            selected_key = res_key
            break
            
    logging.info(f"Selected resolution {selected_res} (key: {selected_key}) for downstream analysis")
    adata.obs['leiden'] = adata.obs[selected_key]

    # Lightweight automatic annotation for unattended runs. Scores are
    # descriptive marker averages; uncertain clusters remain Unknown.
    canonical_markers = {
        'T cells': ['CD3D','CD3E','CD8A','CD4'],
        'B cells': ['CD79A','MS4A1','CD19'],
        'NK cells': ['GNLY','NKG7','FGFBP2'],
        'Monocytes/Macrophages': ['CD14','CD68','AIF1','CSF1R'],
        'Endothelial': ['VWF','PECAM1','CDH5'],
        'Fibroblasts': ['DCN','COL1A1','FBLN1'],
        'SMCs': ['ACTA2','TAGLN','MYH11'],
        'Epithelial': ['EPCAM','KRT18','KRT19'],
        'Mast cells': ['TPSAB1','CPA3'],
        'Plasma cells': ['MZB1','SDC1','IGHG1']
    }
    raw_vars = pd.Index(adata.raw.var_names)
    score_matrix = {}
    for label, genes in canonical_markers.items():
        found = [g for g in genes if g in raw_vars]
        if found:
            idx = [raw_vars.get_loc(g) for g in found]
            vals = adata.raw.X[:, idx]
            vals = vals.toarray() if hasattr(vals, 'toarray') else np.asarray(vals)
            score_matrix[label] = vals.mean(axis=1)
    if score_matrix:
        score_df = pd.DataFrame(score_matrix, index=adata.obs_names)
        score_df.groupby(adata.obs['leiden'], observed=True).mean().to_csv(table_dir / 'candidate_marker_scores_by_cluster.csv')
    adata.obs['cell_type'] = pd.Categorical('Cluster_' + adata.obs['leiden'].astype(str))
    adata.uns['annotation_status'] = 'pending_human_review'

    # 12. Compute Wilcoxon DEG for selected resolution
    logging.info("Computing DEGs (Wilcoxon) for selected clusters")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=True)
    
    deg_df = sc.get.rank_genes_groups_df(adata, group=None)
    deg_df.to_csv(table_dir / 'all_markers_wilcoxon.csv', index=False)
    
    # Top 20 wide format
    top20_dfs = []
    for clst in adata.obs['leiden'].cat.categories:
        clst_df = sc.get.rank_genes_groups_df(adata, group=clst)
        clst_df = clst_df[clst_df['pvals_adj'] < 0.05].head(20)[['names', 'scores']]
        clst_df.columns = [f'Cluster_{clst}_Gene', f'Cluster_{clst}_Score']
        clst_df = clst_df.reset_index(drop=True)
        top20_dfs.append(clst_df)
    
    top20_wide = pd.concat(top20_dfs, axis=1)
    top20_wide.to_csv(table_dir / 'top20_markers_wide.csv', index=False)
    
    # 13. Compute cluster audit table
    logging.info("Computing cluster audit metrics")
    audit_rows = []
    for clst in adata.obs['leiden'].cat.categories:
        clst_adata = adata[adata.obs['leiden'] == clst]
        n_cells = clst_adata.shape[0]
        med_genes = clst_adata.obs['n_genes_by_counts'].median() if 'n_genes_by_counts' in clst_adata.obs else clst_adata.obs.get('n_genes', pd.Series([0], index=clst_adata.obs.index)).median()
        med_umi = clst_adata.obs['total_counts'].median() if 'total_counts' in clst_adata.obs else 0
        med_mt = clst_adata.obs['pct_counts_mt'].median() if 'pct_counts_mt' in clst_adata.obs else clst_adata.obs.get('pct_mt', pd.Series([0], index=clst_adata.obs.index)).median()
        med_doublet = clst_adata.obs['doublet_score'].median() if 'doublet_score' in clst_adata.obs else 0
        
        row = {
            'Cluster': clst,
            'N_Cells': n_cells,
            'Median_Genes': med_genes,
            'Median_UMI': med_umi,
            'Median_MT_Pct': med_mt,
            'Median_DoubletScore': med_doublet
        }
        
        if BATCH_COLUMN in clst_adata.obs:
            sample_counts = clst_adata.obs[BATCH_COLUMN].value_counts()
            for s, c in sample_counts.items():
                row[f'Sample_{s}_Count'] = c
                
        if GROUP_COLUMN in clst_adata.obs:
            group_counts = clst_adata.obs[GROUP_COLUMN].value_counts()
            for g, c in group_counts.items():
                row[f'Group_{g}_Count'] = c
                
        audit_rows.append(row)
        
    audit_df = pd.DataFrame(audit_rows).fillna(0)
    audit_df.to_csv(table_dir / 'cluster_audit.csv', index=False)
    
    # 14. Scan TARGET_GENE expression
    logging.info(f"Scanning target gene: {TARGET_GENE}")
    target_stats = []
    
    if TARGET_GENE in adata.raw.var_names:
        target_idx = np.where(adata.raw.var_names == TARGET_GENE)[0][0]
        expr_array = adata.raw.X[:, target_idx].toarray().flatten() if hasattr(adata.raw.X, 'toarray') else adata.raw.X[:, target_idx]
        
        for clst in adata.obs['leiden'].cat.categories:
            mask = (adata.obs['leiden'] == clst).values
            expr_clst = expr_array[mask]
            
            pos_rate = (expr_clst > 0).mean() if len(expr_clst) > 0 else 0
            mean_expr = expr_clst.mean() if len(expr_clst) > 0 else 0
            
            row = {'Cluster': clst, 'Positive_Rate': pos_rate, 'Mean_Expr': mean_expr}
            
            if GROUP_COLUMN in adata.obs:
                case_mask = mask & (adata.obs[GROUP_COLUMN] == CASE_LABEL).values
                ctrl_mask = mask & (adata.obs[GROUP_COLUMN] == CONTROL_LABEL).values
                
                if case_mask.sum() > 0:
                    row['Case_Mean_Expr'] = expr_array[case_mask].mean()
                if ctrl_mask.sum() > 0:
                    row['Control_Mean_Expr'] = expr_array[ctrl_mask].mean()
                    
            target_stats.append(row)
            
        pd.DataFrame(target_stats).to_csv(table_dir / f'{TARGET_GENE}_cluster_expression.csv', index=False)
    else:
        logging.warning(f"Target gene {TARGET_GENE} not found in the dataset.")
        
    # 15. Generate publication plots
    logging.info("Generating publication plots")
    
    metrics = [m for m in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'doublet_score'] if m in adata.obs]
    qc_metric_panel(adata.obs, 'leiden', fig_dir / 'qc_metrics_by_cluster', RANDOM_SEED, umap_category_palette(adata, 'leiden'))
    if GROUP_COLUMN in adata.obs:
        adata.obs[GROUP_COLUMN] = adata.obs[GROUP_COLUMN].astype('category')
        adata.uns[f'{GROUP_COLUMN}_colors'] = [COLOR_CASE if str(g) == CASE_LABEL else COLOR_CONTROL if str(g) == CONTROL_LABEL else '#999999' for g in adata.obs[GROUP_COLUMN].cat.categories]

    plot_keys = ['leiden', BATCH_COLUMN]
    if GROUP_COLUMN in adata.obs:
        plot_keys.append(GROUP_COLUMN)
        
    for qc in metrics:
        plot_keys.append(qc)
        
    sc.pl.umap(adata, color=plot_keys, wspace=0.4, show=False)
    save_publication_figure(plt.gcf(), fig_dir / 'umap_metadata_and_qc')
    
    res_keys = [f'leiden_r{str(r).replace(".", "_")}' for r in LEIDEN_RESOLUTIONS]
    sc.pl.umap(adata, color=res_keys, wspace=0.4, show=False, ncols=2)
    save_publication_figure(plt.gcf(), fig_dir / 'umap_resolution_grid')
    
    canonical_markers = {
        'T cells': ['CD3D','CD3E','CD8A','CD4'],
        'B cells': ['CD79A','MS4A1','CD19'],
        'NK cells': ['GNLY','NKG7','FGFBP2'],
        'Monocytes/Macrophages': ['CD14','CD68','AIF1','CSF1R'],
        'Endothelial': ['VWF','PECAM1','CDH5'],
        'Fibroblasts': ['DCN','COL1A1','FBLN1'],
        'SMCs': ['ACTA2','TAGLN','MYH11'],
        'Epithelial': ['EPCAM','KRT18','KRT19'],
        'Mast cells': ['TPSAB1','CPA3'],
        'Plasma cells': ['MZB1','SDC1','IGHG1']
    }
    
    plot_markers = {}
    for ct, genes in canonical_markers.items():
        valid_genes = [g for g in genes if g in adata.raw.var_names]
        if valid_genes:
            plot_markers[ct] = valid_genes
            
    sc.pl.dotplot(adata, plot_markers, groupby='leiden', show=False, standard_scale='var')
    save_publication_figure(plt.gcf(), fig_dir / 'markers_canonical_dotplot')
    
    expression_figure_suite(adata, 'leiden', fig_dir, config)

    # 16. Save annotated h5ad checkpoint
    logging.info(f"Saving clustered h5ad to {out_h5ad}")
    from figure_recipes import sanitize_uns_colors
    sanitize_uns_colors(adata)
    adata.write_h5ad(out_h5ad)
    
    # 17. Save JSON report
    report = {
        'PROJECT_NAME': PROJECT_NAME,
        'TARGET_GENE': TARGET_GENE,
        'N_Cells': adata.shape[0],
        'N_HVG_GLOBAL': N_HVG_GLOBAL,
        'N_PCS_GLOBAL': N_PCS_GLOBAL,
        'LEIDEN_RESOLUTIONS': LEIDEN_RESOLUTIONS,
        'Selected_Resolution': selected_res,
        'N_Clusters': len(adata.obs['leiden'].unique()),
        'resolution_selection': 'heuristic_first_10_to_15_clusters',
        'resolution_note': '分辨率由启发式自动选择',
    }
    
    with open(OUTPUT_ROOT / '99_logs' / '02_global_report.json', 'w') as f:
        json.dump(report, f, indent=4)
        
    logging.info("Phase 2 complete!")

if __name__ == '__main__':
    main()
