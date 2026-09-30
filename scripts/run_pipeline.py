#!/usr/bin/env python3
"""Unattended runner for the rapid exploratory scRNA target-gene pipeline.

Each phase is an isolated process so a failed optional phase cannot corrupt
earlier checkpoints. Required phase failures stop dependent work.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

PHASES = [
    ("phase01", "phase01_qc_scrublet.py", ("03_clean.h5ad",)),
    ("phase02", "phase02_global_clustering.py", ("04_global_clustered.h5ad",)),
    ("phase03", "phase03_subclustering.py", ("05_target_subclustered_annotated.h5ad",)),
    ("phase04_sensitivity", "phase05_sensitivity_deg.py", ("07_target_sensitivity_DE_ready.h5ad",)),
]
REJECTED_INPUT_CODE = 3


def main() -> int:
    ap = argparse.ArgumentParser(description="Run all scRNA target-gene phases")
    ap.add_argument("--config", required=True, help="JSON configuration path")
    ap.add_argument("--force", action="store_true", help="do not skip completed phases")
    ap.add_argument("--stop-after", choices=[p[0] for p in PHASES])
    args = ap.parse_args()
    config_path = Path(args.config).resolve()
    cfg = json.loads(config_path.read_text(encoding="utf-8-sig"))
    from config_contract import missing_config
    missing = missing_config(cfg)
    if missing:
        raise SystemExit(f"Missing required config values: {', '.join(missing)}")
    for key in ('INPUT_PATH', 'OUTPUT_ROOT', 'SAMPLE_MANIFEST'):
        if cfg.get(key):
            cfg[key] = str((config_path.parent / Path(cfg[key])).resolve())
    out = Path(cfg['OUTPUT_ROOT'])
    logs = out / "99_logs"
    logs.mkdir(parents=True, exist_ok=True)
    status_path = logs / "pipeline_status.json"
    if status_path.exists() and not args.force:
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
        except Exception:
            status = {"phases": {}}
    else:
        status = {"phases": {}}
    fingerprint = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode())
    for source in sorted(Path(__file__).resolve().parent.glob('*.py')):
        fingerprint.update(source.read_bytes())
    cfg_hash = fingerprint.hexdigest()
    if status.get("config_hash") != cfg_hash:
        status = {"phases": {}}
    status.update({"config": str(config_path), "config_hash": cfg_hash, "started_at": time.strftime("%Y-%m-%dT%H:%M:%S")})
    scripts_dir = Path(__file__).resolve().parent
    resolved_config = logs / 'resolved_config.json'
    resolved_config.write_text(json.dumps(cfg, indent=2), encoding='utf-8')
    import os
    env = os.environ.copy()
    cache = logs / 'runtime_cache'
    cache.mkdir(exist_ok=True)
    env.setdefault('NUMBA_CACHE_DIR', str(cache / 'numba'))
    env.setdefault('MPLCONFIGDIR', str(cache / 'matplotlib'))
    # Windows conda environments may need MKL DLLs discoverable when
    # launched via an absolute python path.
    prefix = Path(sys.prefix)
    libbin = prefix / 'Library' / 'bin'
    if libbin.exists():
        env['PATH'] = str(libbin) + os.pathsep + env.get('PATH', '')
    rerun_downstream = args.force
    for name, script, checkpoints in PHASES:
        entry = status["phases"].setdefault(name, {})
        entry["script"] = script
        entry["started_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        checkpoint_ok = any((out / "05_h5ad" / f).exists() for f in checkpoints)
        if not rerun_downstream and entry.get('returncode') == 0 and checkpoint_ok:
            if name in ('phase02', 'phase03'):
                previous_status = entry.get('status')
                result = subprocess.run([sys.executable, str(scripts_dir / 'human_review.py'), '--checkpoint', str(out / '05_h5ad' / checkpoints[0]), '--phase', name, '--config', str(resolved_config)], env=env)
                entry['status'] = 'needs_review' if result.returncode == 2 else (('quick_provisional' if cfg.get('QUICK_MODE', False) else 'reviewed') if result.returncode == 0 else 'failed')
                status_path.write_text(json.dumps(status, indent=2), encoding='utf-8')
                if result.returncode:
                    return result.returncode
                # Accepted names/cleaned populations invalidate dependent analyses.
                rerun_downstream = previous_status not in ('reviewed', 'quick_provisional')
            if args.stop_after == name:
                break
            continue
        rerun_downstream = True
        cmd = [sys.executable, str(scripts_dir / script), str(resolved_config)]
        try:
            with (logs / f'{name}.log').open('w', encoding='utf-8') as log_file:
                proc = subprocess.run(cmd, cwd=str(scripts_dir.parent), env=env, stdout=log_file, stderr=subprocess.STDOUT, check=False)
            entry["returncode"] = int(proc.returncode)
            checkpoint_ok = any((out / "05_h5ad" / f).exists() for f in checkpoints)
            entry["checkpoint_ok"] = checkpoint_ok
            if proc.returncode == REJECTED_INPUT_CODE:
                entry["status"] = "rejected_input"
            else:
                entry["status"] = "ok" if proc.returncode == 0 and checkpoint_ok else "failed"
        except Exception as exc:
            entry["returncode"] = -1
            entry["status"] = "failed"
            entry["error"] = repr(exc)
        entry["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8")
        if entry["status"] == "rejected_input":
            status["outcome"] = "rejected_input"
            status["analysis_unit"] = "cell"
            status["donor_level_inference"] = False
            status["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8")
            return REJECTED_INPUT_CODE
        if entry["status"] == "failed":
            status["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8")
            return 1
        if name in ('phase02', 'phase03'):
            result = subprocess.run([sys.executable, str(scripts_dir / 'human_review.py'), '--checkpoint', str(out / '05_h5ad' / checkpoints[0]), '--phase', name, '--config', str(resolved_config)], env=env)
            entry['status'] = 'needs_review' if result.returncode == 2 else (('quick_provisional' if cfg.get('QUICK_MODE', False) else 'reviewed') if result.returncode == 0 else 'failed')
            status_path.write_text(json.dumps(status, indent=2), encoding='utf-8')
            if result.returncode:
                return result.returncode
        if args.stop_after == name:
            break
    portfolio = {"status": "disabled"}
    dataset_quality = {"status": "disabled"}
    if cfg.get("FILTER_LOW_INFORMATION_FIGURES", True):
        quality_log = logs / "figure_quality.log"
        quality_cmd = [
            sys.executable,
            str(scripts_dir / "figure_quality.py"),
            "--output-root",
            str(out),
            "--dataset-level",
        ]
        try:
            with quality_log.open("w", encoding="utf-8") as log_file:
                quality_proc = subprocess.run(
                    quality_cmd,
                    cwd=str(scripts_dir.parent),
                    env=env,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            dataset_quality = {
                "status": "ok" if quality_proc.returncode == 0 else ("rejected_dataset" if quality_proc.returncode == 4 else "failed_optional"),
                "returncode": int(quality_proc.returncode),
                "report": str(logs / "figure_quality_dataset_report.json"),
            }
            status["figure_quality"] = dataset_quality
        except Exception as exc:
            dataset_quality = {"status": "failed_optional", "error": repr(exc)}
            status["figure_quality"] = dataset_quality
    if cfg.get("BUILD_PORTFOLIO", True) and dataset_quality.get("status") != "rejected_dataset":
        portfolio_log = logs / "portfolio_board.log"
        portfolio_cmd = [
            sys.executable,
            str(scripts_dir / "portfolio_board.py"),
            "--output-root",
            str(out),
            "--target-gene",
            str(cfg.get("TARGET_GENE", "scRNA")),
        ]
        try:
            with portfolio_log.open("w", encoding="utf-8") as log_file:
                board_proc = subprocess.run(
                    portfolio_cmd,
                    cwd=str(scripts_dir.parent),
                    env=env,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            portfolio = {
                "status": "ok" if board_proc.returncode == 0 else "failed_optional",
                "returncode": int(board_proc.returncode),
                "output_root": str(out / "06_portfolio"),
                "log": str(portfolio_log),
            }
        except Exception as exc:
            portfolio = {"status": "failed_optional", "error": repr(exc), "log": str(portfolio_log)}
    status["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    outcome = f"stopped_after_{args.stop_after}" if args.stop_after else "completed_analysis"
    degradation_reasons = []
    deg_report = logs / "04_deg_report.json"
    if deg_report.exists():
        try:
            deg_info = json.loads(deg_report.read_text(encoding="utf-8"))
            enrichment_status = deg_info.get("enrichment", {}).get("status")
            if enrichment_status in ("no_significant_output", "not_run"):
                degradation_reasons.append(f"enrichment_{enrichment_status}")
            strata = deg_info.get("target_strata_counts", {})
            if int(strata.get("High", 0)) < 5 or int(strata.get("Low", 0)) < 5:
                degradation_reasons.append("insufficient_target_strata_for_some_cell_level_contrasts")
        except Exception:
            degradation_reasons.append("deg_report_unreadable")
    audit_report = logs.parent / "00_metadata" / "ingestion_audit_report.json"
    if audit_report.exists():
        try:
            audit_info = json.loads(audit_report.read_text(encoding="utf-8"))
            if not audit_info.get("metadata", {}).get("has_case_label", False) or not audit_info.get("metadata", {}).get("has_control_label", False):
                degradation_reasons.append("missing_case_or_control_label")
        except Exception:
            pass
    if degradation_reasons and not args.stop_after:
        outcome = "completed_with_degradation"
    status["outcome"] = outcome
    status["degradation_reasons"] = sorted(set(degradation_reasons))
    status["analysis_unit"] = "cell"
    status["donor_level_inference"] = False
    status["portfolio"] = portfolio
    status_path.write_text(json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
