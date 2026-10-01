"""The seeded demo matrix runs all four quicklook phases without a network."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "scripts" / "run_pipeline.py"
DEMO = ROOT / "scripts" / "write_demo_h5ad.py"
_GUARD = """
import socket
_connect = socket.socket.connect

def _connect_offline(self, address):
    if isinstance(address, tuple) and address:
        host = str(address[0])
        if host not in {"127.0.0.1", "::1", "localhost"}:
            raise OSError(101, "network disabled: %s" % host)
    return _connect(self, address)

socket.socket.connect = _connect_offline
"""


def test_demo_matrix_finishes_four_phases_offline(tmp_path: Path) -> None:
    guard = tmp_path / "guard"
    guard.mkdir()
    (guard / "sitecustomize.py").write_text(_GUARD, encoding="utf-8")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(guard) + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONNOUSERSITE"] = "1"
    h5ad = tmp_path / "demo.h5ad"
    wrote = subprocess.run(
        [sys.executable, str(DEMO), "--out", str(h5ad), "--cells", "80"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert wrote.returncode == 0, wrote.stderr
    import anndata as ad

    names = [str(name) for name in ad.read_h5ad(h5ad).var_names]
    gmt = tmp_path / "local.gmt"
    gmt.write_text(
        "SET_A\tNA\t" + "\t".join(names[:8]) + "\n" + "SET_B\tNA\t" + "\t".join(names[8:16]) + "\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    out.mkdir()
    config = {
        "TARGET_GENE": "IFITM3",
        "INPUT_PATH": str(h5ad),
        "INPUT_FORMAT": "h5ad",
        "OUTPUT_ROOT": str(out),
        "QUICK_MODE": True,
        "RANDOM_SEED": 42,
        "FILTER_LOW_INFORMATION_FIGURES": True,
        "BUILD_PORTFOLIO": True,
        "GROUP_COLUMN": "group",
        "CASE_LABEL": "Disease",
        "CONTROL_LABEL": "Healthy",
        "BATCH_COLUMN": "sample_id",
        "LOCAL_GMT": str(gmt),
    }
    (out / "config.json").write_text(json.dumps(config), encoding="utf-8")
    ran = subprocess.run(
        [sys.executable, str(PIPELINE), "--config", str(out / "config.json")],
        cwd=str(ROOT),
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    status_path = out / "99_logs" / "pipeline_status.json"
    assert ran.returncode == 0, status_path.read_text(encoding="utf-8") if status_path.is_file() else ran.stderr
    status = json.loads(status_path.read_text(encoding="utf-8"))
    assert status["outcome"] == "completed_analysis"
    assert status["analysis_unit"] == "cell"
    assert status["donor_level_inference"] is False
    assert [phase["status"] for phase in status["phases"].values()] == [
        "ok",
        "quick_provisional",
        "quick_provisional",
        "ok",
    ]
    for folder in ("01_qc", "02_global", "03_target", "04_deg"):
        assert list((out / "04_figures" / folder).glob("*.png")), folder
    for name in (
        "volcano_All_cells_IFITM3_high_vs_low.png",
        "volcano_Disease_Case_IFITM3_high_vs_low.png",
        "volcano_Healthy_Control_IFITM3_high_vs_low.png",
        "gsea_case_pathway_nes.png",
    ):
        assert (out / "04_figures" / "04_deg" / name).is_file()
    report = json.loads((out / "99_logs" / "04_deg_report.json").read_text(encoding="utf-8"))
    assert report["deg_contrasts"]["Case"]["status"] == "completed"
    assert report["deg_contrasts"]["Control"]["status"] == "completed"
    assert report["enrichment"]["gsea_status"] == "completed"
    quality = json.loads(
        (out / "99_logs" / "figure_quality_dataset_report.json").read_text(encoding="utf-8")
    )
    assert quality["status"] == "reject"
    assert "low_information_figure_present" in quality["reasons"]
