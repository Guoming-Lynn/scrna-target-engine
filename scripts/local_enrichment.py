"""Local gene-set enrichment. Gene lists stay on this machine."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom


class GseaNotRun(RuntimeError):
    """GSEA was not run. The message is the report text."""


def gmt_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def load_gmt(path: Path) -> dict[str, set[str]]:
    """Read a GMT file. Column 2 is always the description, never a gene."""
    sets: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 3 or not parts[0].strip():
            continue
        genes = {item.strip() for item in parts[2:] if item.strip()}
        sets[parts[0].strip()] = genes
    return sets


def local_ora_table(genes: list[str], gmt_path: Path, universe: list[str]) -> pd.DataFrame:
    """Hypergeometric overlap. The background is the supplied universe, not the GMT."""
    library = load_gmt(gmt_path)
    background = {gene.casefold() for gene in universe}
    query = {gene.casefold() for gene in genes} & background
    n_query = len(query)
    n_universe = len(background)
    rows: list[dict[str, object]] = []
    if n_query == 0 or n_universe == 0:
        return pd.DataFrame(columns=["Term", "Adjusted P-value"])
    for name, members in library.items():
        in_background = {member.casefold() for member in members} & background
        overlap = len(query & in_background)
        if overlap == 0:
            continue
        p_value = float(hypergeom.sf(overlap - 1, n_universe, len(in_background), n_query))
        rows.append({"Term": name, "Overlap": overlap, "Adjusted P-value": p_value})
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=["Term", "Adjusted P-value"])
    frame["Adjusted P-value"] = _bh([float(value) for value in frame["Adjusted P-value"]])
    return frame.sort_values("Adjusted P-value")


def local_gsea_table(
    ranked: pd.DataFrame,
    gmt_path: Path,
    permutations: int = 1000,
    seed: int = 42,
) -> pd.DataFrame:
    """Pre-ranked GSEA on a local GMT file. The same seed repeats the same table."""
    try:
        import gseapy as gp
    except ImportError as exc:
        raise GseaNotRun("未运行") from exc
    result = gp.prerank(
        rnk=ranked.iloc[:, :2],
        gene_sets=str(gmt_path),
        permutation_num=permutations,
        seed=seed,
        threads=1,
        min_size=1,
        max_size=5000,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    table = result.res2d.copy()
    kept = _gsea_columns(table)
    return kept.sort_values("FDR q-val")


def _gsea_columns(table: pd.DataFrame) -> pd.DataFrame:
    folded = {
        column.casefold().replace(" ", "").replace("_", "").replace("-", ""): column for column in table.columns
    }

    def pick(*keys: str) -> str | None:
        for key in keys:
            found = folded.get(key)
            if found is not None:
                return found
        return None

    term = pick("term", "name")
    nes = pick("nes")
    fdr = pick("fdrqval", "fdr")
    if term is None or nes is None or fdr is None:
        raise GseaNotRun("未运行")
    kept = table[[term, nes, fdr]].copy()
    kept.columns = ["Term", "NES", "FDR q-val"]
    return kept


def _bh(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=lambda index: p_values[index])
    adjusted = [1.0] * count
    running = 1.0
    for rank in range(count, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted
