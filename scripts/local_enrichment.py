"""Local gene-set enrichment. Gene lists stay on this machine."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import pandas as pd


def gmt_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def load_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = [item.strip() for item in line.split("\t") if item.strip()]
        if len(parts) < 3:
            continue
        sets[parts[0]] = set(parts[2:]) if parts[1].startswith("http") else set(parts[1:])
    return sets


def local_ora_table(genes: list[str], gmt_path: Path) -> pd.DataFrame:
    """Hypergeometric overlap against a local GMT file."""
    library = load_gmt(gmt_path)
    query = {gene.casefold() for gene in genes}
    universe: set[str] = set()
    folded: dict[str, set[str]] = {}
    for name, members in library.items():
        names = {member.casefold() for member in members}
        folded[name] = names
        universe.update(names)
    query &= universe
    n_query = len(query)
    n_universe = len(universe)
    rows: list[dict[str, object]] = []
    if n_query == 0 or n_universe == 0:
        return pd.DataFrame(columns=["Term", "Adjusted P-value"])
    for name, members in folded.items():
        overlap = len(query & members)
        if overlap == 0:
            continue
        rows.append(
            {
                "Term": name,
                "Overlap": overlap,
                "Adjusted P-value": _hypergeom_sf(overlap, n_universe, len(members), n_query),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=["Term", "Adjusted P-value"])
    frame["Adjusted P-value"] = _bh([float(value) for value in frame["Adjusted P-value"]])
    return frame.sort_values("Adjusted P-value")


def local_gsea_table(ranked: pd.DataFrame, gmt_path: Path, permutations: int = 100) -> pd.DataFrame:
    """Score local GMT sets from a ranked name/score table. No network call."""
    library = load_gmt(gmt_path)
    names = [str(value) for value in ranked.iloc[:, 0]]
    scores = [float(value) for value in ranked.iloc[:, 1]]
    order = sorted(range(len(names)), key=lambda index: scores[index], reverse=True)
    ordered_names = [names[index].casefold() for index in order]
    ordered_scores = [scores[index] for index in order]
    rows: list[dict[str, object]] = []
    for name, members in library.items():
        folded = {member.casefold() for member in members}
        observed = _running_sum(ordered_names, ordered_scores, folded)
        greater = 0
        for seed in range(permutations):
            shuffled = _shift(ordered_names, seed)
            if abs(_running_sum(shuffled, ordered_scores, folded)) >= abs(observed):
                greater += 1
        rows.append({"Term": name, "NES": observed, "FDR q-val": (greater + 1) / (permutations + 1)})
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=["Term", "NES", "FDR q-val"])
    frame["FDR q-val"] = _bh([float(value) for value in frame["FDR q-val"]])
    return frame.sort_values("FDR q-val")


def _hypergeom_sf(overlap: int, universe: int, set_size: int, query: int) -> float:
    top = min(query, set_size)
    if overlap > top or universe <= 0:
        return 1.0
    denominator = math.comb(universe, query)
    if denominator == 0:
        return 1.0
    total = 0.0
    for count in range(overlap, top + 1):
        rest = query - count
        if set_size < count or universe - set_size < rest:
            continue
        total += math.comb(set_size, count) * math.comb(universe - set_size, rest) / denominator
    return min(1.0, total)


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


def _running_sum(names: list[str], scores: list[float], members: set[str]) -> float:
    hits = [abs(score) if name in members else 0.0 for name, score in zip(names, scores)]
    hit_total = sum(hits)
    if hit_total == 0:
        return 0.0
    miss = sum(1 for name in names if name not in members)
    if miss == 0:
        return 0.0
    total = 0.0
    peak = 0.0
    for name, weight in zip(names, hits):
        total += weight / hit_total if name in members else -1.0 / miss
        peak = total if abs(total) > abs(peak) else peak
    return peak


def _shift(names: list[str], seed: int) -> list[str]:
    if not names:
        return []
    offset = seed % len(names)
    return names[offset:] + names[:offset]
