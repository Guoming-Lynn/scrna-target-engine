import sys
import json
import logging
import traceback
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
from scipy import sparse as sp

# Set matplotlib backend to 'Agg' for headless rendering
matplotlib.use('Agg')

# QC plot colors come from plot_utils.
from plot_utils import COLOR_CASE as CASE_COLOR, COLOR_CONTROL as CONTROL_COLOR
from figure_recipes import ranked_violin, cell_count_bars, qc_scatter, scrublet_comparison, scrublet_histograms
from config_contract import batch_column
from input_checks import input_rejection
from pipeline_utils import resolve_target_gene


from phase01_scrublet import run_sample_scrublet
from phase01_support import (
    InputRejected,
    load_input,
    plot_qc_violin,
    robust_median_mad,
    save_fig,
    setup_logger,
    write_ingestion_audit,
)

# Main Execution
# ----------------------------------------------------------------------

def main():
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    config_path = Path(config_path)
    if not config_path.exists():
        sys.stderr.write(f"Error: Config file not found at {config_path}\n")
        sys.exit(1)

    with open(config_path, "r", encoding="utf-8-sig") as f:
        config = json.load(f)

    # 1. Parse config variables
    TARGET_GENE = config.get("TARGET_GENE", "UNKNOWN")
    PROJECT_NAME = config.get("PROJECT_NAME", "scRNA_Project")
    INPUT_PATH = Path(config.get("INPUT_PATH", "."))
    INPUT_FORMAT = config.get("INPUT_FORMAT", "h5ad")
    OUTPUT_ROOT = Path(config.get("OUTPUT_ROOT", "./output"))
    SAMPLE_MANIFEST = config.get("SAMPLE_MANIFEST", None)
    CASE_LABEL = config.get("CASE_LABEL", "Case")
    CONTROL_LABEL = config.get("CONTROL_LABEL", "Control")
    BATCH_COLUMN = batch_column(config)
    RANDOM_SEED = config.get("RANDOM_SEED", 42)
    
    QC_MT_MAX = config.get("QC_MT_MAX", 20.0)
    QC_HB_MAX = config.get("QC_HB_MAX", 5.0)
    QC_MIN_GENES = config.get("QC_MIN_GENES", 200)
    QC_MAX_GENES = config.get("QC_MAX_GENES", 6000)
    
    SCRUBLET_RATE = config.get("SCRUBLET_RATE", 0.06)
    SCRUBLET_MAX_THRESHOLD = config.get("SCRUBLET_MAX_THRESHOLD", 0.25)
    SKIP_SCRUBLET = bool(config.get("SKIP_SCRUBLET", False))
    audit = {
        "status": "started",
        "project_name": PROJECT_NAME,
        "target_gene_requested": TARGET_GENE,
        "input_path": str(INPUT_PATH),
        "input_format": INPUT_FORMAT,
        "declared_context": config.get("DATASET_CONTEXT", {}),
        "analysis_unit": "cell",
        "donor_level_inference": False,
    }

    # Set random seed
    np.random.seed(RANDOM_SEED)

    # 2. Create output directories
    dir_metadata = OUTPUT_ROOT / "00_metadata"
    dir_qc_tables = OUTPUT_ROOT / "03_tables" / "01_qc"
    dir_scrub_tables = OUTPUT_ROOT / "03_tables" / "01_qc" / "scrublet"
    dir_qc_figs = OUTPUT_ROOT / "04_figures" / "01_qc"
    dir_scrub_figs = OUTPUT_ROOT / "04_figures" / "01_qc" / "scrublet"
    dir_h5ad = OUTPUT_ROOT / "05_h5ad"
    dir_logs = OUTPUT_ROOT / "99_logs"

    for d in [dir_metadata, dir_qc_tables, dir_scrub_tables, dir_qc_figs, dir_scrub_figs, dir_h5ad, dir_logs]:
        d.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(dir_logs / f"{PROJECT_NAME}_phase01.log")
    logger.info(f"Starting Phase 1 (QC & Scrublet) for project {PROJECT_NAME}")
    logger.info(f"Target Gene: {TARGET_GENE}")
    logger.info(f"Input path: {INPUT_PATH} ({INPUT_FORMAT})")
    
    palette = {CASE_LABEL: CASE_COLOR, CONTROL_LABEL: CONTROL_COLOR}

    try:
        # 3. Load Data
        logger.info(f"Loading data from {INPUT_PATH}...")
        if not INPUT_PATH.exists():
            raise InputRejected(f"Input path does not exist: {INPUT_PATH}")
        adata = load_input(INPUT_PATH, INPUT_FORMAT, config)
            
        logger.info(f"Raw data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
        adata.var_names_make_unique()
        matrix_min = adata.X.data.min() if sp.issparse(adata.X) and adata.X.nnz else np.nanmin(adata.X)
        resolved_target = resolve_target_gene(adata.var_names, TARGET_GENE)
        rejection = input_rejection(
            n_cells=int(adata.n_obs),
            n_genes=int(adata.n_vars),
            min_cells=int(config.get("MIN_INPUT_CELLS", 20)),
            min_genes=int(config.get("MIN_INPUT_GENES", 50)),
            numeric=bool(np.issubdtype(adata.X.dtype, np.number)),
            minimum=float(matrix_min),
            target_resolved=resolved_target,
            target_requested=TARGET_GENE,
        )
        if rejection:
            raise InputRejected(rejection)
        TARGET_GENE = resolved_target
        audit["target_gene_resolved"] = resolved_target
        audit["input_matrix"] = {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "unique_gene_names_after_deduplication": bool(adata.var_names.is_unique),
            "sparse": bool(hasattr(adata.X, "tocsr")),
        }

        # 4. Load Manifest and map metadata
        if SAMPLE_MANIFEST and Path(SAMPLE_MANIFEST).exists():
            logger.info(f"Loading sample manifest from {SAMPLE_MANIFEST}")
            manifest = pd.read_csv(SAMPLE_MANIFEST)
            
            # Identify columns for merging. 
            # We assume manifest has 'barcode_suffix' and 'sample_id', 'group' etc.
            if 'barcode_suffix' in manifest.columns and INPUT_FORMAT != "h5ad":
                logger.info("Mapping cells to manifest using barcode suffixes...")
                # Extract suffix from index e.g., 'ACGT-1' -> '1'
                obs = adata.obs.copy()
                obs['_cell_id'] = obs.index.astype(str)
                obs['barcode_suffix'] = [idx.split('-')[-1] if '-' in idx else '' for idx in obs['_cell_id']]
                manifest = manifest.copy()
                manifest['barcode_suffix'] = manifest['barcode_suffix'].astype(str)
                obs['barcode_suffix'] = obs['barcode_suffix'].astype(str)
                adata.obs = obs.merge(manifest, on='barcode_suffix', how='left', suffixes=('', '_manifest'))
                adata.obs.index = adata.obs.pop('_cell_id').astype(str)
            else:
                logger.info("Merging manifest by sample ID...")
                # Try to merge on a common key, assume 'sample_id'
                if 'sample_id' in manifest.columns and 'sample_id' in adata.obs.columns:
                    # To keep original index
                    obs_idx = adata.obs.index
                    adata.obs = adata.obs.merge(manifest, on='sample_id', how='left')
                    adata.obs.index = obs_idx
                else:
                    logger.warning("Could not map manifest securely (missing matching keys). Proceeding with existing obs.")
        else:
            logger.warning("No SAMPLE_MANIFEST provided or found. Proceeding with existing metadata.")
            if 'group' not in adata.obs.columns:
                adata.obs['group'] = 'Unknown'
            if 'sample_id' not in adata.obs.columns:
                adata.obs['sample_id'] = "Sample_1"

        # Fill NAs
        adata.obs['group'] = adata.obs['group'].astype(str).replace('nan', 'Unknown')
        adata.obs['sample_id'] = adata.obs['sample_id'].astype(str).replace('nan', 'Unknown')
        audit["metadata"] = {
            "group_column": "group",
            "group_values": sorted(adata.obs['group'].unique().tolist()),
            "sample_column": "sample_id",
            "n_samples": int(adata.obs['sample_id'].nunique()),
            "has_case_label": bool(CASE_LABEL in set(adata.obs['group'])),
            "has_control_label": bool(CONTROL_LABEL in set(adata.obs['group'])),
        }
        audit["status"] = "accepted"
        write_ingestion_audit(dir_metadata / "ingestion_audit_report.json", audit)
        adata.uns["ingestion_audit"] = audit
        adata.uns["target_gene_resolved"] = resolved_target

        # 5. Set unique cell IDs & Store raw counts
        if INPUT_FORMAT != "h5ad" or 'sample_id' in adata.obs.columns:
            logger.info("Setting unique cell IDs...")
            # create sample_id_barcode
            barcodes = [idx.split('-')[0] if '-' in idx else idx for idx in adata.obs.index]
            adata.obs.index = [f"{s}_{b}" for s, b in zip(adata.obs['sample_id'], barcodes)]
            adata.obs_names_make_unique()

        logger.info("Storing raw counts in layers['counts']...")
        adata.layers["counts"] = adata.X.copy()
        
        # Save raw checkpoint
        raw_checkpoint = dir_h5ad / "01_raw_merged.h5ad"
        logger.info(f"Saving raw merged dataset to {raw_checkpoint}")
        adata.write_h5ad(raw_checkpoint)

        # 6. Compute QC metrics
        logger.info("Computing QC metrics (MT, Ribo, HB)...")
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')
        adata.var['hb'] = adata.var_names.str.match('^HB[ABDEGMQZ]')
        
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt', 'ribo', 'hb'], 
            percent_top=None, log1p=False, inplace=True
        )

        # Rename standard columns for consistency
        qc_cols_map = {
            'n_genes_by_counts': 'n_genes',
            'total_counts': 'total_counts',
            'pct_counts_mt': 'pct_mt',
            'pct_counts_ribo': 'pct_ribo',
            'pct_counts_hb': 'pct_hb'
        }
        adata.obs.rename(columns=qc_cols_map, inplace=True)
        # Keep Scanpy canonical names for downstream phases while exposing short aliases.
        adata.obs['n_genes_by_counts'] = adata.obs['n_genes']
        adata.obs['pct_counts_mt'] = adata.obs['pct_mt']
        adata.obs['pct_counts_ribo'] = adata.obs['pct_ribo']
        adata.obs['pct_counts_hb'] = adata.obs['pct_hb']
        adata.obs['log10_total_counts'] = np.log10(adata.obs['total_counts'] + 1)

        adata.uns['plot_seed'] = RANDOM_SEED
        raw_samples = sorted(adata.obs['sample_id'].astype(str).unique())
        cell_count_bars(adata.obs, dir_qc_figs / 'cell_counts_raw', CASE_LABEL, CONTROL_LABEL, title='Raw cells per sample', sample_order=raw_samples)
        scatter_limits = (0., float(adata.obs['pct_mt'].max()))
        scatter_xy = (float(adata.obs['total_counts'].max())*1.03, float(adata.obs['n_genes'].max())*1.03, max(float(adata.obs['total_counts'].min())*.9, .1))

        # 7. Pre-filter plots
        logger.info("Generating pre-filter QC plots...")
        for metric in ['n_genes', 'log10_total_counts', 'pct_mt', 'pct_ribo', 'pct_hb']:
            if metric in adata.obs.columns:
                plot_qc_violin(adata, metric, 'group', dir_qc_figs, f"pre_QC_{metric}", palette)

        qc_scatter(adata.obs, dir_qc_figs / "qc_pre_counts_vs_genes", "Pre-QC", color_limits=scatter_limits, xy_limits=scatter_xy)

        # 8. Run sample-wise adaptive MAD + hard-cap QC
        logger.info("Running sample-wise adaptive MAD QC filtering...")
        
        qc_stats = []
        adata.obs['pass_qc'] = True
        adata.obs['qc_fail_reason'] = ""

        samples = adata.obs['sample_id'].unique()
        
        for sample in samples:
            mask = adata.obs['sample_id'] == sample
            df_samp = adata.obs.loc[mask]
            
            med_g, _, mad_g = robust_median_mad(df_samp['n_genes'])
            med_mt, _, mad_mt = robust_median_mad(df_samp['pct_mt'])
            med_hb, _, mad_hb = robust_median_mad(df_samp['pct_hb'])
            med_rb, _, mad_rb = robust_median_mad(df_samp['pct_ribo'])

            # Thresholds
            min_g = max(QC_MIN_GENES, med_g - 3 * mad_g)
            max_g = min(QC_MAX_GENES, med_g + 4 * mad_g)
            max_mt = min(QC_MT_MAX, med_mt + 3 * mad_mt)
            max_hb = min(QC_HB_MAX, med_hb + 3 * mad_hb)
            max_rb = min(60.0, med_rb + 5 * mad_rb)

            pass_genes = (df_samp['n_genes'] >= min_g) & (df_samp['n_genes'] <= max_g)
            pass_mt = df_samp['pct_mt'] <= max_mt
            pass_hb = df_samp['pct_hb'] <= max_hb
            pass_rb = df_samp['pct_ribo'] <= max_rb

            pass_all = pass_genes & pass_mt & pass_hb & pass_rb
            
            adata.obs.loc[mask, 'pass_qc'] = pass_all
            
            # Record reasons
            reasons = []
            for idx in df_samp.index:
                r = []
                if not pass_genes.loc[idx]: r.append("genes")
                if not pass_mt.loc[idx]: r.append("mt")
                if not pass_hb.loc[idx]: r.append("hb")
                if not pass_rb.loc[idx]: r.append("ribo")
                reasons.append(",".join(r) if r else "pass")
            adata.obs.loc[mask, 'qc_fail_reason'] = reasons

            qc_stats.append({
                'sample_id': sample,
                'total_cells': len(df_samp),
                'pass_cells': pass_all.sum(),
                'min_genes_thresh': min_g,
                'max_genes_thresh': max_g,
                'max_mt_thresh': max_mt,
                'max_hb_thresh': max_hb,
                'max_ribo_thresh': max_rb
            })
            
        stats_df = pd.DataFrame(qc_stats)
        stats_df.to_csv(dir_qc_tables / "qc_thresholds_and_stats.csv", index=False)
        logger.info(f"Saved QC stats to {dir_qc_tables / 'qc_thresholds_and_stats.csv'}")

        # Filter
        adata = adata[adata.obs['pass_qc']].copy()
        logger.info(f"Cells remaining after QC: {adata.n_obs}")

        qc_checkpoint = dir_h5ad / "02_after_QC.h5ad"
        logger.info(f"Saving post-QC dataset to {qc_checkpoint}")
        adata.write_h5ad(qc_checkpoint)

        # 9. Post-filter plots
        logger.info("Generating post-filter QC plots...")
        for metric in ['n_genes', 'log10_total_counts', 'pct_mt', 'pct_ribo', 'pct_hb']:
            if metric in adata.obs.columns:
                plot_qc_violin(adata, metric, 'group', dir_qc_figs, f"post_QC_{metric}", palette)

        qc_scatter(adata.obs, dir_qc_figs / "qc_post_counts_vs_genes", "Post-QC", color_limits=scatter_limits, xy_limits=scatter_xy)

        # 10. Scrublet sample-wise
        logger.info("Running Scrublet sample-wise...")
        scrublet_res, distributions = run_sample_scrublet(
            adata,
            samples,
            skip=SKIP_SCRUBLET,
            rate=SCRUBLET_RATE,
            max_threshold=SCRUBLET_MAX_THRESHOLD,
            seed=int(RANDOM_SEED),
            table_dir=dir_scrub_tables,
            logger=logger,
        )
        scrublet_histograms(distributions, scrublet_res, dir_scrub_figs / 'scrublet_score_distributions')

        df_scrublet = pd.DataFrame(scrublet_res)
        df_scrublet.to_csv(dir_scrub_tables / "scrublet_stats.csv", index=False)
        
        # Persist scores before removing doublets for reproducible figure-only reruns.
        adata.obs.to_csv(dir_scrub_tables / 'post_qc_cell_scores.csv')
        adata.write_h5ad(qc_checkpoint)
        before_scrublet = adata.obs.copy()

        # Remove doublets
        adata = adata[~adata.obs['predicted_doublet']].copy()
        logger.info(f"Cells remaining after Scrublet (clean): {adata.n_obs}")
        
        clean_checkpoint = dir_h5ad / "03_clean.h5ad"
        logger.info(f"Saving clean dataset to {clean_checkpoint}")
        adata.write_h5ad(clean_checkpoint)

        # 11. Final counts and original QC-versus-clean doublet-score panels.
        cell_count_bars(adata.obs, dir_qc_figs / 'cell_counts_clean', CASE_LABEL, CONTROL_LABEL, sample_order=raw_samples)
        if not SKIP_SCRUBLET and np.isfinite(before_scrublet['doublet_score']).any():
            scrublet_comparison(before_scrublet, adata.obs, dir_scrub_figs / 'scrublet_before_after_counts_vs_genes')

        # 12. Save JSON report
        report = {
            "project_name": PROJECT_NAME,
            "target_gene": TARGET_GENE,
            "timestamp": datetime.now().isoformat(),
            "final_cells": int(adata.n_obs),
            "final_genes": int(adata.n_vars),
            "sample_stats": stats_df.to_dict(orient='records'),
            "scrublet_stats": df_scrublet.to_dict(orient='records')
        }
        
        report_path = dir_logs / "phase01_report.json"
        with open(report_path, "w") as f:
            json.dump(report, f, indent=4)
        logger.info(f"Saved run report to {report_path}")
        logger.info("Phase 1 complete!")

    except InputRejected as e:
        audit["status"] = "rejected_input"
        audit["reason"] = str(e)
        write_ingestion_audit(dir_metadata / "ingestion_audit_report.json", audit)
        logger.error(f"Input rejected: {e}")
        sys.exit(3)
    except Exception as e:
        logger.error(f"Error in Phase 1 pipeline: {str(e)}")
        logger.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
