"""Write the studio demo matrix. Counts are synthetic and seeded."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SOURCE = "合成矩阵，由 write_demo_h5ad.py 用种子 42 生成。不是策展数据，不能当作正式分析。"
GENE = "IFITM3"


def write_demo(path: Path) -> None:
    """Write 80 cells and 250 genes. Each cell expresses every gene, so QC can keep them."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(42)
    n_cells = 80
    n_genes = 250
    counts = rng.integers(1, 6, size=(n_cells, n_genes), dtype=np.int32)
    names = [f"GENE{index:03d}" for index in range(n_genes - 1)]
    names.append(GENE)
    counts[:40, -1] = counts[:40, -1] + 30
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    write_demo(args.out)
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
