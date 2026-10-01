"""Write the studio demo matrix. Counts are synthetic and seeded."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SOURCE = "合成矩阵，由 write_demo_h5ad.py 用种子 42 生成。不是策展数据，不能当作正式分析。"
GENE = "IFITM3"
SHOWCASE_CELLS = 1600
REGRESSION_CELLS = 80


def write_demo(path: Path, n_cells: int = SHOWCASE_CELLS) -> None:
    """Write a seeded synthetic matrix. The default is large enough for figure review.

    ``80`` keeps the fast regression matrix. A larger IFITM3 gap than 12 splits
    the target cells below the subclustering minimum.
    """
    if n_cells == REGRESSION_CELLS:
        _write_regression(path)
        return
    if n_cells != SHOWCASE_CELLS:
        raise ValueError("n_cells must be 80 or 1600")
    _write_showcase(path)


def _write_regression(path: Path) -> None:
    """80 cells. Every gene is detected, so QC can keep the matrix."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(42)
    n_cells = REGRESSION_CELLS
    n_genes = 250
    counts = rng.integers(1, 6, size=(n_cells, n_genes), dtype=np.int32)
    names = [f"GENE{index:03d}" for index in range(n_genes - 1)]
    names.append(GENE)
    high = np.array([*range(20), *range(40, 60)])
    counts[high, -1] = counts[high, -1] + 12
    obs = pd.DataFrame(
        {
            "group": ["Disease"] * 40 + ["Healthy"] * 40,
            "sample_id": ["S1"] * 40 + ["S2"] * 40,
        }
    )
    data = ad.AnnData(X=counts.astype(np.float32), obs=obs, var=pd.DataFrame(index=names))
    data.uns["demo_source"] = SOURCE
    path.parent.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(path)


def _write_showcase(path: Path) -> None:
    """1600 cells in two programs. The target program has two mild subtypes.

    A stronger subtype gap splits the overview, and the phase-3 violin then
    sees only one category. QC marker genes vary so the panels are not constant.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(42)
    n_cells = SHOWCASE_CELLS
    n_genes = 250
    counts = rng.integers(1, 6, size=(n_cells, n_genes), dtype=np.int32)
    names = [f"GENE{index:03d}" for index in range(n_genes - 1)]
    names.append(GENE)
    names[0] = "MT-CO1"
    names[1] = "RPS3"
    names[2] = "HBB"
    counts[:, 0] = rng.integers(0, 7, size=n_cells)
    counts[:, 1] = rng.integers(0, 5, size=n_cells)
    counts[:, 2] = rng.integers(0, 3, size=n_cells)
    counts[0:800, 10:30] += 80
    counts[800:1600, 30:50] += 80
    counts[0:400, 50:70] += 1
    counts[400:800, 70:90] += 1
    # Exactly one quarter of the target program sits above the IFITM3 Q75 cut.
    # A few companion genes move with them so the volcano is not only the target.
    high = np.array([index for start in (0, 200, 400, 600) for index in range(start, start + 50)])
    counts[high, -1] += 12
    counts[high, 100:108] += 2
    group: list[str] = []
    sample: list[str] = []
    for _block in range(4):
        group.extend(["Disease"] * 200 + ["Healthy"] * 200)
        sample.extend(["S1"] * 200 + ["S2"] * 200)
    obs = pd.DataFrame({"group": group, "sample_id": sample})
    data = ad.AnnData(X=counts.astype(np.float32), obs=obs, var=pd.DataFrame(index=names))
    data.uns["demo_source"] = SOURCE
    path.parent.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=SHOWCASE_CELLS)
    args = parser.parse_args(argv)
    try:
        write_demo(args.out, args.cells)
    except ValueError as exc:
        sys.stderr.write(str(exc) + "\n")
        return 2
    sys.stdout.write(
        json.dumps(
            {
                "source": SOURCE,
                "gene": GENE,
                "group": "group",
                "case": "Disease",
                "control": "Healthy",
                "batch": "sample_id",
            }
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
