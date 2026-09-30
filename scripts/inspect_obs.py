"""List obs column names and a few unique values. The matrix is not copied out."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

UNIQUE_LIMIT = 30
H5AD_NOTE = "只有 h5ad 能列出 obs 列。没有分组时不出 case/control 对比图。"


def inspect_h5ad(path: Path, limit: int = UNIQUE_LIMIT) -> dict[str, object]:
    """Return column names from an h5ad obs table, or a note for other paths."""
    if path.suffix.casefold() != ".h5ad" or not path.is_file():
        return {"columns": [], "note": H5AD_NOTE}
    adata = ad.read_h5ad(path, backed="r")
    try:
        columns: list[dict[str, object]] = []
        for name in list(adata.obs.columns):
            values = adata.obs[name].drop_duplicates().astype(str).head(limit).tolist()
            columns.append({"name": str(name), "values": values})
        return {"columns": columns, "note": ""}
    finally:
        closer = getattr(getattr(adata, "file", None), "close", None)
        if callable(closer):
            closer()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="inspect_obs")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        payload = inspect_h5ad(args.input)
    except OSError as exc:
        payload = {"columns": [], "note": f"无法读取输入：{exc}"}
        _write(args.out, payload)
        return 1
    _write(args.out, payload)
    return 0


def _write(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
