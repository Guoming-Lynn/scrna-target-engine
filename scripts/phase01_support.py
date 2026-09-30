"""Input loading and small QC helpers for phase 01."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.figure import Figure


class InputRejected(RuntimeError):
    """The supplied object is not a usable expression matrix for this run."""


def load_input(input_path: Path, input_format: str, config: dict):
    """Load supported formats, including ordinary GEO text matrices."""
    if input_format == "10x_mtx":
        return sc.read_10x_mtx(input_path, var_names="gene_symbols", make_unique=True)
    if input_format == "h5ad":
        return sc.read_h5ad(input_path)
    if input_format == "h5":
        return sc.read_10x_h5(input_path)
    if input_format in {"csv", "tsv", "txt", "dense_csv", "dense_tsv"}:
        return _text_matrix(input_path, input_format, config)
    raise InputRejected(f"Unsupported INPUT_FORMAT: {input_format}")


def _text_matrix(input_path: Path, input_format: str, config: dict):
    sep = "," if input_format in {"csv", "dense_csv"} else "\t"
    frame = pd.read_csv(input_path, sep=sep, index_col=0, compression="infer")
    orientation = str(config.get("TEXT_MATRIX_ORIENTATION", "genes_by_cells")).lower()
    if orientation == "cells_by_genes":
        frame = frame.T
    elif orientation == "auto":
        frame = _maybe_transpose(frame, str(config.get("TARGET_GENE", "")))
    frame = frame.apply(pd.to_numeric, errors="coerce")
    if frame.isna().any().any():
        raise InputRejected("Text matrix contains non-numeric expression values")
    return ad.AnnData(
        X=frame.to_numpy(dtype=float).T,
        obs=pd.DataFrame(index=frame.columns.astype(str)),
        var=pd.DataFrame(index=frame.index.astype(str)),
    )


def _maybe_transpose(frame: pd.DataFrame, requested: str) -> pd.DataFrame:
    row_match = requested in frame.index or requested.casefold() in {str(x).casefold() for x in frame.index}
    col_match = requested in frame.columns or requested.casefold() in {str(x).casefold() for x in frame.columns}
    if col_match and not row_match:
        return frame.T
    if not row_match and not col_match and frame.shape[0] < frame.shape[1]:
        return frame.T
    return frame


def write_ingestion_audit(path: Path, audit: dict) -> None:
    path.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")


def setup_logger(log_file: Path) -> logging.Logger:
    """Setup a basic logger to file and stdout."""
    logger = logging.getLogger("phase01_qc_scrublet")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    if logger.handlers:
        logger.handlers.clear()
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    return logger


def save_fig(fig: Figure, output_path: Path, filename_base: str) -> None:
    from plot_utils import save_publication_figure

    save_publication_figure(fig, output_path / filename_base)


def plot_qc_violin(adata, metric, group_col, output_dir, filename, palette, ylabel=""):
    """Condition-colored full-data violin with white boxes and 100 sampled dots."""
    import matplotlib.pyplot as plt
    from figure_recipes import ranked_violin
    from plot_utils import save_publication_figure

    fig, ax = plt.subplots(figsize=(9, 5))
    ranked_violin(
        adata.obs,
        group_col,
        metric,
        ax,
        palette=palette,
        seed=int(adata.uns.get("plot_seed", 42)),
    )
    if ylabel:
        ax.set_ylabel(ylabel)
    save_publication_figure(fig, Path(output_dir) / filename)


def robust_median_mad(data):
    """Return median, raw MAD, and a scaled MAD that falls back to the IQR."""
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    med = np.median(data) if data.size else np.nan
    mad_raw = np.median(np.abs(data - med)) if data.size else np.nan
    mad = 1.4826 * mad_raw if data.size else np.nan
    if not np.isfinite(mad) or mad <= 1e-9:
        q75, q25 = np.percentile(data, [75, 25])
        iqr = q75 - q25
        mad = iqr / 1.349 if iqr > 1e-9 else 1e-3
    return med, mad_raw, mad
