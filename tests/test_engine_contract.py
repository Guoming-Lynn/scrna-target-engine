"""Offline enrichment, figure rejection, and config checks."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from config_contract import DEFAULT_BATCH_COLUMN, batch_column, missing_config
from figure_quality import quality_reasons
from input_checks import input_rejection
from local_enrichment import load_gmt, local_gsea_table, local_ora_table
from resolution_note import RESOLUTION_NOTE, choose_resolution

ROOT = Path(__file__).resolve().parents[1]
PHASE05 = ROOT / "scripts" / "phase05_sensitivity_deg.py"
PHASE02 = ROOT / "scripts" / "phase02_global_clustering.py"


def test_phase05_does_not_call_network_enrichment() -> None:
    text = PHASE05.read_text(encoding="utf-8").casefold()
    assert ".enrichr(" not in text
    assert "gp.enrichr" not in text
    assert "local_gmt" in text


def test_phase02_labels_heuristic_resolution() -> None:
    text = PHASE02.read_text(encoding="utf-8")
    assert RESOLUTION_NOTE in text
    assert choose_resolution([(0.2, 4), (0.4, 12), (0.8, 20)]) == 0.4
    assert choose_resolution([(0.2, 4), (0.3, 6)]) == 0.2


def test_config_rejects_missing_and_placeholder_values() -> None:
    assert missing_config({"INPUT_PATH": "a", "OUTPUT_ROOT": "b"}) == ["TARGET_GENE"]
    assert missing_config(
        {"TARGET_GENE": "YOUR_GENE", "INPUT_PATH": "a", "OUTPUT_ROOT": "b"}
    ) == ["TARGET_GENE"]
    assert missing_config({"TARGET_GENE": "IFITM3", "INPUT_PATH": "a", "OUTPUT_ROOT": "b"}) == []


@pytest.mark.parametrize(
    ("audit", "reason"),
    [
        ({"kind": "ranked_violin", "n_input": 0, "groups": []}, "empty_violin"),
        (
            {"kind": "ranked_violin", "n_input": 3, "groups": [{"density_drawn": True}]},
            "single_category_violin",
        ),
        (
            {
                "kind": "ranked_violin",
                "n_input": 3,
                "groups": [{"density_drawn": False}, {"density_drawn": False}],
            },
            "all_constant_violin",
        ),
        ({"kind": "feature", "n_cells": 0, "positive_cells": 0}, "empty_feature_panel"),
        ({"kind": "feature", "n_cells": 10, "positive_cells": 0}, "no_positive_target_expression"),
        ({"kind": "volcano", "n_plotted": 0, "n_up": 0, "n_down": 0}, "empty_volcano"),
        ({"kind": "volcano", "n_plotted": 8, "n_up": 0, "n_down": 0}, "no_significant_genes"),
    ],
)
def test_each_figure_rejection_reason(audit: dict[str, object], reason: str) -> None:
    assert reason in quality_reasons([audit])


def test_input_rejection_paths() -> None:
    assert input_rejection(
        n_cells=10,
        n_genes=100,
        min_cells=20,
        min_genes=50,
        numeric=True,
        minimum=0,
        target_resolved="IFITM3",
        target_requested="IFITM3",
    ) == "Expression matrix has too few cells (10)"
    assert "negative" in str(
        input_rejection(
            n_cells=100,
            n_genes=100,
            min_cells=20,
            min_genes=50,
            numeric=True,
            minimum=-1,
            target_resolved="IFITM3",
            target_requested="IFITM3",
        )
    )
    assert input_rejection(
        n_cells=100,
        n_genes=100,
        min_cells=20,
        min_genes=50,
        numeric=True,
        minimum=0,
        target_resolved=None,
        target_requested="NOGENE",
    )


def test_local_ora_uses_only_the_gmt_file(tmp_path: Path) -> None:
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "HALLMARK_A\tdesc\tGENE1\tGENE2\tGENE3\tGENE4\tGENE5\n"
        "HALLMARK_B\tdesc\tOTHER1\tOTHER2\tOTHER3\n",
        encoding="utf-8",
    )
    universe = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "OTHER1", "OTHER2", "OTHER3", "BG1"]
    table = local_ora_table(["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"], gmt, universe)
    assert table.iloc[0]["Term"] == "HALLMARK_A"
    assert "Adjusted P-value" in table.columns


def test_gmt_description_is_not_a_gene(tmp_path: Path) -> None:
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("SET\tNA\tGENE1\tGENE2\n", encoding="utf-8")
    members = {gene.casefold() for gene in load_gmt(gmt)["SET"]}
    assert members == {"gene1", "gene2"}


def test_ora_p_value_matches_scipy_on_the_supplied_universe(tmp_path: Path) -> None:
    from scipy.stats import hypergeom

    gmt = tmp_path / "sets.gmt"
    genes = [f"G{index}" for index in range(20)]
    gmt.write_text("SET\tNA\t" + "\t".join(genes[:10]) + "\n", encoding="utf-8")
    query = genes[:5]
    table = local_ora_table(query, gmt, genes)
    expected = float(hypergeom.sf(4, 20, 10, 5))
    assert table.iloc[0]["Adjusted P-value"] == pytest.approx(expected)
    wider = local_ora_table(query, gmt, genes + ["EXTRA"])
    assert float(wider.iloc[0]["Adjusted P-value"]) != pytest.approx(expected)


def test_gsea_repeats_with_the_same_seed(tmp_path: Path) -> None:
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("SET_A\tNA\tA\tB\tC\nSET_B\tNA\tX\tY\tZ\n", encoding="utf-8")
    ranked = pd.DataFrame(
        {"names": ["A", "B", "C", "X", "Y", "Z"], "score": [3, 2, 1, -1, -2, -3]}
    )
    first = local_gsea_table(ranked, gmt, permutations=50, seed=42)
    second = local_gsea_table(ranked, gmt, permutations=50, seed=42)
    assert first.reset_index(drop=True).equals(second.reset_index(drop=True))
    assert "FDR q-val" in first.columns


def test_phase05_does_not_call_online_enrichment_method() -> None:
    text = PHASE05.read_text(encoding="utf-8")
    assert ".enrichr(" not in text


def test_batch_column_default_is_defined_once() -> None:
    assert batch_column({}) == DEFAULT_BATCH_COLUMN
    assert DEFAULT_BATCH_COLUMN == "sample_id"
    assert batch_column({"BATCH_COLUMN": "donor"}) == "donor"
    assert batch_column({"BATCH_COLUMN": ""}) == "sample_id"
    for name in (
        "phase01_qc_scrublet.py",
        "phase02_global_clustering.py",
        "phase03_subclustering.py",
    ):
        text = (ROOT / "scripts" / name).read_text(encoding="utf-8")
        assert 'config.get("BATCH_COLUMN"' not in text
        assert "config.get('BATCH_COLUMN'" not in text
        assert "batch_column(config)" in text


def test_inspect_h5ad_returns_obs_columns(tmp_path: Path) -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from inspect_obs import main

    path = tmp_path / "tiny.h5ad"
    ad.AnnData(
        X=np.zeros((2, 1)),
        obs=pd.DataFrame({"group": ["Disease", "Healthy"], "sample_id": ["S1", "S2"]}),
        var=pd.DataFrame(index=["IFITM3"]),
    ).write_h5ad(path)
    out = tmp_path / "inspect.json"
    assert main(["--input", str(path), "--out", str(out)]) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    columns = {item["name"]: item["values"] for item in payload["columns"]}
    assert columns["group"] == ["Disease", "Healthy"]
    assert columns["sample_id"] == ["S1", "S2"]


def test_inspect_non_h5ad_lists_no_columns(tmp_path: Path) -> None:
    from inspect_obs import inspect_h5ad

    path = tmp_path / "counts.csv"
    path.write_text("gene,a\nIFITM3,1\n", encoding="utf-8")
    payload = inspect_h5ad(path)
    assert payload["columns"] == []
    assert "没有分组" in str(payload["note"])


def test_demo_matrix_is_synthetic_and_names_its_source(tmp_path: Path) -> None:
    import anndata as ad
    from write_demo_h5ad import SOURCE, main

    path = tmp_path / "demo.h5ad"
    assert main(["--out", str(path)]) == 0
    data = ad.read_h5ad(path)
    assert data.n_obs == 80
    assert data.n_vars == 250
    assert "IFITM3" in list(data.var_names)
    assert float(data.X.min()) >= 0
    assert data.uns["demo_source"] == SOURCE
    assert "策展" in SOURCE


def test_scrublet_call_passes_configured_random_state() -> None:
    source = (ROOT / "scripts" / "phase01_scrublet.py").read_text(encoding="utf-8")
    call = source.split("Scrublet(", 1)[1].split(")", 1)[0]
    assert "random_state=" in call
    phase01 = (ROOT / "scripts" / "phase01_qc_scrublet.py").read_text(encoding="utf-8")
    assert "seed=int(RANDOM_SEED)" in phase01


def test_ora_bh_includes_zero_overlap_sets_in_the_universe(tmp_path: Path) -> None:
    from scipy.stats import hypergeom

    gmt = tmp_path / "sets.gmt"
    genes = [f"G{index}" for index in range(20)]
    gmt.write_text(
        "HIT\tNA\t" + "\t".join(genes[:10]) + "\n"
        "MISS\tNA\t" + "\t".join(genes[10:15]) + "\n"
        "OUTSIDE\tNA\tNOT_IN_UNIVERSE\n",
        encoding="utf-8",
    )
    table = local_ora_table(genes[:5], gmt, genes).set_index("Term")
    raw = float(hypergeom.sf(4, 20, 10, 5))
    assert set(table.index) == {"HIT", "MISS"}
    assert float(table.loc["MISS", "Adjusted P-value"]) == pytest.approx(1.0)
    assert float(table.loc["HIT", "Adjusted P-value"]) == pytest.approx(min(1.0, raw * 2))


def test_synthetic_10x_h5_loads_gene_symbols(tmp_path: Path) -> None:
    import h5py
    import numpy as np
    from phase01_support import load_input

    path = tmp_path / "filtered_feature_bc_matrix.h5"
    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("barcodes", data=[b"AAAC-1", b"AAAT-1"])
        matrix.create_dataset("data", data=np.array([1, 3], dtype=np.int32))
        matrix.create_dataset("indices", data=np.array([0, 1], dtype=np.int32))
        matrix.create_dataset("indptr", data=np.array([0, 1, 2], dtype=np.int32))
        matrix.create_dataset("shape", data=np.array([2, 2], dtype=np.int32))
        features = matrix.create_group("features")
        features.create_dataset("id", data=[b"ENSG1", b"ENSG2"])
        features.create_dataset("name", data=[b"IFITM3", b"ACTB"])
        features.create_dataset("feature_type", data=[b"Gene Expression", b"Gene Expression"])
        features.create_dataset("genome", data=[b"GRCh38", b"GRCh38"])
    loaded = load_input(path, "h5", {})
    assert loaded.n_obs == 2
    assert list(loaded.var_names) == ["IFITM3", "ACTB"]
    assert float(loaded.X[0, 0]) == 1.0


def test_harmony_cluster_count_never_asks_for_one_cluster() -> None:
    from pipeline_utils import harmony_cluster_count

    assert harmony_cluster_count(40) == 2
    assert harmony_cluster_count(80) == 3
    assert harmony_cluster_count(10_000) == 100


def test_local_gsea_returns_fdr_column(tmp_path: Path) -> None:
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("SET_A\tdesc\tA\tB\tC\nSET_B\tdesc\tX\tY\tZ\n", encoding="utf-8")
    ranked = pd.DataFrame({"names": ["A", "B", "C", "X", "Y", "Z"], "score": [3, 2, 1, -1, -2, -3]})
    table = local_gsea_table(ranked, gmt, permutations=5)
    assert "FDR q-val" in table.columns
    assert set(table["Term"]) == {"SET_A", "SET_B"}
