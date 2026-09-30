"""Sample-wise Scrublet scoring used by phase 01."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scrublet as scr


def run_sample_scrublet(
    adata,
    samples,
    *,
    skip: bool,
    rate: float,
    max_threshold: float,
    table_dir: Path,
    logger: logging.Logger,
):
    """Score each sample. The caller still draws the summary histogram."""
    scrublet_res: list[dict] = []
    distributions: dict = {}
    adata.obs["doublet_score"] = np.nan
    adata.obs["scrublet_status"] = "not_evaluated"
    adata.obs["predicted_doublet"] = False
    for sample in samples:
        if skip:
            logger.warning("SKIP_SCRUBLET=true; retaining all cells after QC")
            break
        mask = adata.obs["sample_id"] == sample
        sub_adata = adata[mask]
        if len(sub_adata) < 50:
            logger.warning("Sample %s has too few cells (%s) for Scrublet. Skipping.", sample, len(sub_adata))
            continue
        scrub = scr.Scrublet(sub_adata.X, expected_doublet_rate=rate)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
        )
        auto_thresh = scrub.threshold_
        thresh = auto_thresh
        if thresh is None or not np.isfinite(thresh) or thresh > max_threshold:
            thresh = max_threshold
            predicted_doublets = doublet_scores > thresh
        adata.obs.loc[mask, "scrublet_status"] = "evaluated"
        adata.obs.loc[mask, "doublet_score"] = doublet_scores
        adata.obs.loc[mask, "predicted_doublet"] = predicted_doublets
        num_doublets = int(predicted_doublets.sum())
        scrublet_res.append(
            {
                "sample_id": sample,
                "total_cells": len(sub_adata),
                "doublets": num_doublets,
                "doublet_rate": num_doublets / len(sub_adata),
                "auto_threshold": auto_thresh,
                "effective_threshold": thresh,
                "threshold_mode": "auto" if auto_thresh == thresh else "capped",
                "threshold": thresh,
            }
        )
        simulated = np.asarray(scrub.doublet_scores_sim_)
        distributions[sample] = (doublet_scores, simulated)
        safe_sample = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(sample))
        pd.DataFrame(
            {"observed": pd.Series(doublet_scores), "simulated": pd.Series(simulated)}
        ).to_csv(table_dir / f"{safe_sample}_score_distributions.csv", index=False)
    return scrublet_res, distributions
