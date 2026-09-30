# -*- coding: utf-8 -*-
"""
Stage 04: Target Gene Sensitivity DEG & Pathway Enrichment Analysis
=================================================================
1. Stratifies cells by Target Gene expression:
   - Binary: Negative (counts = 0) vs Positive (counts >= 1).
   - High vs Low within positives: Low (< Q75 log1p) vs High (>= Q75 log1p).
2. Subtype-wise Target Gene strata composition stacked barplot.
3. Case-internal DEG: {TARGET_GENE}_High vs {TARGET_GENE}_Low within Case group.
4. Control-internal DEG: {TARGET_GENE}_High vs {TARGET_GENE}_Low within Control group.
5. Identifies Case-specific DEGs (significant in Case but non-significant/discordant in Control).
6. Publication-grade Volcano Plots with labeled top DEGs.
7. ORA against a local GMT file. Gene lists are not sent to a network service.
8. Pre-ranked scores against the same local GMT file.
9. Exports full tables, 300 DPI PNG + PDF plots, and JSON run report.
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

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from local_enrichment import GseaNotRun, gmt_sha256, local_gsea_table, local_ora_table
from scientific_plots import plot_ora_results, plot_gsea_results

# Import local publication plotting utilities
from plot_utils import (
    set_publication_style,
    save_publication_figure,
    plot_volcano,
    COLOR_CASE,
    COLOR_CONTROL,
    PALETTE_6
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def run_local_ora(gene_list, table_prefix, figure_prefix, gmt_path, universe):
    if gmt_path is None or not Path(gmt_path).is_file():
        logging.warning("Local GMT is missing. Enrichment was not sent anywhere.")
        return None
    if not gene_list or len(gene_list) < 5:
        logging.warning("Gene list too small for ORA enrichment (< 5 genes).")
        return None
    res = local_ora_table(list(gene_list), Path(gmt_path), list(universe))
    if res is None or res.empty:
        return None
    res.to_csv(table_prefix.with_suffix('.csv'), index=False, encoding='utf-8-sig')
    try:
        plot_ora_results(res, figure_prefix)
    except Exception as exc:
        logging.warning("ORA plot failed: %s", exc)
    return res


def run_local_gsea(rnk_df, table_prefix, figure_prefix, gmt_path, seed):
    if gmt_path is None or not Path(gmt_path).is_file():
        logging.warning("Local GMT is missing. Enrichment was not sent anywhere.")
        return None
    logging.info("Running local pre-ranked scores against %s", gmt_path)
    res = local_gsea_table(rnk_df, Path(gmt_path), permutations=1000, seed=int(seed))
    if res is None or res.empty:
        return None
    res.to_csv(table_prefix.with_suffix('.csv'), index=False, encoding='utf-8-sig')
    try:
        plot_gsea_results(res, figure_prefix)
    except Exception as exc:
        logging.warning("GSEA plot failed: %s", exc)
    return res


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
    FDR_CUTOFF = config.get("DEG_FDR_CUTOFF", 0.05)
    FC_CUTOFF = config.get("DEG_LOG2FC_CUTOFF", 0.25)
    GMT_PATH = config.get("LOCAL_GMT")
    RANDOM_SEED = int(config.get("RANDOM_SEED", 42))

    set_publication_style()

    dirs = {
        'tables': OUTPUT_ROOT / '03_tables' / '04_deg',
        'figures': OUTPUT_ROOT / '04_figures' / '04_deg',
        'h5ad': OUTPUT_ROOT / '05_h5ad',
        'logs': OUTPUT_ROOT / '99_logs'
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    def write_deg_status(label, status, reason, n_high=0, n_low=0):
        payload = {
            'status': status, 'reason': reason, 'label': label,
            'target_gene': TARGET_GENE, 'n_high': int(n_high), 'n_low': int(n_low)
        }
        path = dirs['tables'] / f"deg_{TARGET_GENE}_{label}_status.json"
        path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding='utf-8')
        # A skipped volcano is still audited so dataset-level quality gating can reject it.
        audit = {'kind': 'volcano', 'n_plotted': 0, 'n_up': 0, 'n_down': 0,
                 'skip_reason': reason, 'status': status}
        (dirs['figures'] / f"volcano_{label}_{TARGET_GENE}_high_vs_low.plot.json").write_text(
            json.dumps(audit, indent=2, ensure_ascii=False), encoding='utf-8')
        return payload

    # 1. Load the final target-cell checkpoint. Pseudotime is not part of this pipeline.
    input_h5ad = dirs['h5ad'] / "05_target_subclustered_annotated.h5ad"
    if not input_h5ad.exists():
        candidates = list(dirs['h5ad'].glob("*.h5ad"))
        if candidates:
            input_h5ad = candidates[-1]
        else:
            logging.error(f"Cannot find input h5ad in {dirs['h5ad']}")
            sys.exit(1)

    logging.info(f"Loading AnnData from {input_h5ad} ...")
    adata = sc.read_h5ad(input_h5ad)
    logging.info(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    TARGET_GENE = str(adata.uns.get('target_gene_resolved', TARGET_GENE))

    # Check raw vs X
    src_vars = adata.raw.var_names if adata.raw is not None else adata.var_names
    X_src = adata.raw.X if adata.raw is not None else adata.X

    if TARGET_GENE not in src_vars:
        logging.error(f"Target gene {TARGET_GENE} not found in AnnData!")
        sys.exit(1)

    # 2. Stratification: Extract Expression & Strata
    idx_t = src_vars.get_loc(TARGET_GENE)
    col_t = X_src[:, idx_t]
    t_vals = col_t.toarray().ravel() if sp.issparse(col_t) else np.asarray(col_t).ravel()
    adata.obs[f'{TARGET_GENE}_expr'] = t_vals

    is_positive = t_vals > 0
    adata.obs['target_binary'] = np.where(is_positive, 'Positive', 'Negative')

    pos_vals = t_vals[is_positive]
    q75 = np.percentile(pos_vals, 75) if len(pos_vals) > 0 else 0.0

    conditions = [
        t_vals <= 0,
        (t_vals > 0) & (t_vals < q75),
        t_vals >= q75
    ]
    choices = ['Negative', 'Low', 'High']
    adata.obs['target_strata'] = np.select(conditions, choices, default='Negative')

    logging.info(f"{TARGET_GENE} Strata counts:")
    logging.info(f"{adata.obs['target_strata'].value_counts().to_dict()}")

    # 3. Subtype Strata Composition Barplot
    cluster_key = "subtype" if "subtype" in adata.obs else [c for c in adata.obs.columns if "celltype" in c or "leiden" in c][0]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    strata_comp = pd.crosstab(adata.obs[cluster_key], adata.obs['target_strata'], normalize='index').reindex(columns=['Negative', 'Low', 'High'], fill_value=0) * 100
    strata_comp.plot(
        kind='bar',
        stacked=True,
        color=['#d9d9d9', '#fdae6b', '#d94801'],
        edgecolor='black',
        linewidth=0.5,
        ax=ax
    )
    ax.set_title(f"{TARGET_GENE} Expression Strata across Subtypes", fontsize=12.5, fontweight='bold', pad=10)
    ax.set_xlabel("Subtype", fontsize=10.5)
    ax.set_ylabel("Proportion (%)", fontsize=10.5)
    ax.legend(title=f"{TARGET_GENE} Level", bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, fontsize=9.5)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xticks(rotation=30, ha='right', fontsize=9)
    save_publication_figure(fig, dirs['figures'] / "composition_target_expression_strata")

    # 4. Differential Expression Function
    def run_strata_deg(sub_adata, label):
        logging.info(f"Running Wilcoxon DEG: {TARGET_GENE}_High vs {TARGET_GENE}_Low in {label} ...")
        # Subset to Low and High only
        valid_cells = sub_adata.obs['target_strata'].isin(['High', 'Low'])
        adata_deg = sub_adata[valid_cells].copy()

        n_high = (adata_deg.obs['target_strata'] == 'High').sum()
        n_low = (adata_deg.obs['target_strata'] == 'Low').sum()
        if n_high < 5 or n_low < 5:
            logging.warning(f"Not enough High/Low cells in {label} (High: {n_high}, Low: {n_low}). Skipping DEG.")
            return pd.DataFrame(columns=['names','scores','logfoldchanges','pvals','pvals_adj','log2FC','p_val_adj']), write_deg_status(label, 'skipped', 'insufficient_high_low_cells', n_high, n_low)

        try:
            sc.tl.rank_genes_groups(adata_deg, groupby='target_strata', groups=['High'], reference='Low', method='wilcoxon', use_raw=True)
            res = sc.get.rank_genes_groups_df(adata_deg, group='High')
        except Exception as exc:
            logging.warning("DEG failed for %s: %s", label, exc)
            return pd.DataFrame(columns=['names','scores','logfoldchanges','pvals','pvals_adj','log2FC','p_val_adj']), write_deg_status(label, 'skipped', f'rank_genes_groups_failed: {exc}', n_high, n_low)
        if res is None or res.empty or 'logfoldchanges' not in res.columns or 'pvals_adj' not in res.columns:
            return pd.DataFrame(columns=['names','scores','logfoldchanges','pvals','pvals_adj','log2FC','p_val_adj']), write_deg_status(label, 'skipped', 'empty_or_malformed_deg_result', n_high, n_low)
        res['log2FC'] = res['logfoldchanges']
        res['p_val_adj'] = res['pvals_adj']
        res.to_csv(dirs['tables'] / f"deg_{TARGET_GENE}_{label}_high_vs_low_all.csv", index=False, encoding='utf-8-sig')
        res_interpretive = res[res['names'].map(lambda value: str(value).casefold() != TARGET_GENE.casefold())].copy()
        res_interpretive.to_csv(dirs['tables'] / f"deg_{TARGET_GENE}_{label}_high_vs_low_interpretive.csv", index=False, encoding='utf-8-sig')

        # Volcano plot
        plot_volcano(
            deg_df=res_interpretive,
            gene_col='names',
            log2fc_col='log2FC',
            pval_col='p_val_adj',
            out_prefix=dirs['figures'] / f"volcano_{label}_{TARGET_GENE}_high_vs_low",
            title=f"{label}: {TARGET_GENE}-High vs -Low Volcano Plot",
            fc_cutoff=FC_CUTOFF,
            p_cutoff=FDR_CUTOFF,
            exclude_genes=TARGET_GENE
        )
        return res, {'status': 'completed', 'label': label, 'n_high': int(n_high), 'n_low': int(n_low)}

    # 5. Overall, case-internal and control-internal exploratory DEGs
    res_all, audit_all = run_strata_deg(adata, "All_cells")
    res_case = pd.DataFrame(); audit_case = {'status': 'not_run'}
    res_ctrl = pd.DataFrame(); audit_ctrl = {'status': 'not_run'}
    if GROUP_COL in adata.obs:
        adata_case = adata[adata.obs[GROUP_COL] == CASE_LABEL].copy()
        adata_ctrl = adata[adata.obs[GROUP_COL] == CONTROL_LABEL].copy()
        res_case, audit_case = run_strata_deg(adata_case, f"{CASE_LABEL}_Case")
        res_ctrl, audit_ctrl = run_strata_deg(adata_ctrl, f"{CONTROL_LABEL}_Control")

    # 6. Pathway Enrichment (ORA & GSEA) on Case DEGs
    enrichment_report = {
        'backend': 'local_gmt',
        'status': 'not_run',
        'target_excluded': True,
        'gmt_sha256': gmt_sha256(Path(GMT_PATH)) if GMT_PATH and Path(GMT_PATH).is_file() else None,
    }
    target_rows = pd.DataFrame()
    sig_case = pd.DataFrame()
    if not res_case.empty and {'p_val_adj','log2FC','names'}.issubset(res_case.columns):
        sig_case_all = res_case[(res_case['p_val_adj'] < FDR_CUTOFF) & (res_case['log2FC'].abs() >= FC_CUTOFF)]
        # An empty boolean index drops the columns on this pandas build.
        if sig_case_all.empty:
            sig_case = sig_case_all.copy()
            target_rows = sig_case_all.copy()
        else:
            target_mask = sig_case_all['names'].map(
                lambda value: str(value).casefold() == TARGET_GENE.casefold()
            )
            sig_case = sig_case_all.loc[~target_mask].copy()
            target_rows = sig_case_all.loc[target_mask].copy()
        sig_case.to_csv(dirs['tables'] / f"deg_{CASE_LABEL}_significant.csv", index=False, encoding='utf-8-sig')

        up_genes = sig_case[sig_case['log2FC'] >= FC_CUTOFF]['names'].tolist()
        down_genes = sig_case[sig_case['log2FC'] <= -FC_CUTOFF]['names'].tolist()

        logging.info(f"Case DEGs: {len(up_genes)} Up-regulated, {len(down_genes)} Down-regulated.")

        universe = res_case["names"].astype(str).tolist()
        ora_up = run_local_ora(up_genes, dirs['tables'] / "ora_case_up", dirs['figures'] / "ora_case_up_pathways", GMT_PATH, universe)
        ora_down = run_local_ora(down_genes, dirs['tables'] / "ora_case_down", dirs['figures'] / "ora_case_down_pathways", GMT_PATH, universe)

        # GSEA Pre-ranked
        rnk = res_case[res_case['names'].map(lambda value: str(value).casefold() != TARGET_GENE.casefold())][['names', 'log2FC', 'p_val_adj']].copy()
        rnk['score'] = np.sign(rnk['log2FC']) * -np.log10(np.clip(rnk['p_val_adj'], 1e-300, 1.0))
        rnk_df = rnk[['names', 'score']].dropna().sort_values('score', ascending=False)
        rnk_df.to_csv(dirs['tables'] / "gsea_rank_input.rnk", sep='\t', index=False, header=False)
        try:
            gsea_res = run_local_gsea(rnk_df, dirs['tables'] / "gsea_results", dirs['figures'] / "gsea_case_pathway_nes", GMT_PATH, RANDOM_SEED)
            gsea_status = "completed" if gsea_res is not None else "no_significant_output"
        except GseaNotRun:
            logging.warning("GSEA 未运行")
            gsea_res = None
            gsea_status = "未运行"
        enrichment_report.update({
            'status': 'completed' if (ora_up is not None or ora_down is not None or gsea_res is not None) else 'no_significant_output',
            'n_case_sig_deg_including_target': int(len(sig_case_all)),
            'n_case_sig_deg_excluding_target': int(len(sig_case)),
            'target_defining_sig_rows': int(len(target_rows)),
            'ora_up_rows': int(len(ora_up)) if ora_up is not None else 0,
            'ora_up_significant_fdr_rows': int((ora_up['Adjusted P-value'] < FDR_CUTOFF).sum()) if ora_up is not None and 'Adjusted P-value' in ora_up else 0,
            'ora_down_rows': int(len(ora_down)) if ora_down is not None else 0,
            'ora_down_significant_fdr_rows': int((ora_down['Adjusted P-value'] < FDR_CUTOFF).sum()) if ora_down is not None and 'Adjusted P-value' in ora_down else 0,
            'gsea_rows': int(len(gsea_res)) if gsea_res is not None else 0,
            'gsea_significant_fdr_rows': int((gsea_res['FDR q-val'] < FDR_CUTOFF).sum()) if gsea_res is not None and 'FDR q-val' in gsea_res else 0,
            'gsea_status': gsea_status,
        })

    # 7. Save Checkpoint with Strata Annotations
    final_h5ad = dirs['h5ad'] / "07_target_sensitivity_DE_ready.h5ad"
    logging.info(f"Saving final sensitivity AnnData -> {final_h5ad} ...")
    from figure_recipes import sanitize_uns_colors
    sanitize_uns_colors(adata)
    adata.write_h5ad(final_h5ad, compression='gzip')

    # 8. Save JSON Run Report
    strata_by_sample = {}
    if 'sample_id' in adata.obs:
        strata_by_sample = {
            str(sample): {str(level): int(count) for level, count in row.items()}
            for sample, row in pd.crosstab(adata.obs['sample_id'], adata.obs['target_strata']).iterrows()
        }
    report = {
        "step": "05_sensitivity_deg",
        "target_gene": TARGET_GENE,
        "n_cells": adata.n_obs,
        "n_target_positive": int(is_positive.sum()),
        "pct_target_positive": float(is_positive.mean() * 100),
        "target_strata_counts": {str(k): int(v) for k, v in adata.obs['target_strata'].value_counts().items()},
        "target_strata_by_sample": strata_by_sample,
        "analysis_unit": "cell",
        "donor_level_inference": False,
        "target_gene_excluded_from_interpretive_deg_and_enrichment": True,
        "n_target_defining_sig_rows": int(len(target_rows)),
        "fdr_cutoff": FDR_CUTOFF,
        "log2fc_cutoff": FC_CUTOFF,
        "n_sig_case_deg": int(len(sig_case)),
        "deg_contrasts": {"All_cells": audit_all, "Case": audit_case, "Control": audit_ctrl},
        "enrichment": enrichment_report,
    }
    with open(dirs['logs'] / "04_deg_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    logging.info("[DONE] Stage 04 Target Sensitivity DEG & Enrichment completed successfully!")


if __name__ == "__main__":
    main()
