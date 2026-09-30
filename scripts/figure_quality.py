#!/usr/bin/env python3
"""Reject entire datasets whose figure set contains low-information panels.

The filter is intentionally evidence-based. It acts on plot audit sidecars,
never on expression matrices or statistical tables. The dataset-level command
moves a rejected project as one unit so the curated area never mixes a good
panel from a bad run with panels from another run.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from figure_output import write_figure_index


def _audits(path: Path) -> list[dict]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return []
    if isinstance(value, dict):
        return [value]
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    return []


def quality_reasons(audits: list[dict]) -> list[str]:
    """Return deterministic rejection reasons for an audited figure bundle."""
    reasons: list[str] = []
    for audit in audits:
        kind = str(audit.get("kind", ""))
        if kind == "ranked_violin":
            groups = audit.get("groups") or []
            n_input = int(audit.get("n_input") or 0)
            if n_input <= 0 or not groups:
                reasons.append("empty_violin")
            elif len(groups) < 2:
                reasons.append("single_category_violin")
            elif not any(bool(group.get("density_drawn")) for group in groups):
                reasons.append("all_constant_violin")
        elif kind == "feature":
            n_cells = int(audit.get("n_cells") or 0)
            positive = int(audit.get("positive_cells") or 0)
            if n_cells <= 0:
                reasons.append("empty_feature_panel")
            elif positive <= 0:
                reasons.append("no_positive_target_expression")
        elif kind == "volcano":
            plotted = int(audit.get("n_plotted") or 0)
            significant = int(audit.get("n_up") or 0) + int(audit.get("n_down") or 0)
            if plotted <= 0:
                reasons.append("empty_volcano")
            elif significant <= 0:
                reasons.append("no_significant_genes")
    return sorted(set(reasons))


def _bundle_files(sidecar: Path) -> list[Path]:
    name = sidecar.name
    if not name.endswith(".plot.json"):
        return [sidecar]
    stem = name[: -len(".plot.json")]
    extensions = (
        f"{stem}.png", f"{stem}.pdf", f"{stem}.svg", f"{stem}.plot.json",
        f"{stem}.plot_data.csv", f"{stem}.layout.json", f"{stem}.caption.txt",
        f"{stem}.csv", f"{stem}.json",
    )
    return [sidecar.parent / item for item in extensions if (sidecar.parent / item).exists()]


def filter_output(output_root: Path, quarantine: bool = True) -> dict:
    output_root = output_root.resolve()
    figure_root = output_root / "04_figures"
    if not figure_root.exists():
        return {"status": "no_figures", "rejected": []}
    rejected: list[dict] = []
    for sidecar in sorted(figure_root.rglob("*.plot.json")):
        reasons = quality_reasons(_audits(sidecar))
        if not reasons:
            continue
        rel = sidecar.relative_to(figure_root)
        files = _bundle_files(sidecar)
        target_dir = output_root / "07_figure_quarantine" / rel.parent if quarantine else None
        if target_dir:
            target_dir.mkdir(parents=True, exist_ok=True)
        moved: list[str] = []
        for source in files:
            if target_dir:
                destination = target_dir / source.name
                if destination.exists():
                    destination.unlink()
                shutil.move(str(source), str(destination))
                moved.append(str(destination.relative_to(output_root)))
            else:
                source.unlink()
                moved.append(str(source.relative_to(output_root)))
        rejected.append({
            "figure": str(rel.with_suffix("")),
            "source_sidecar": str(rel),
            "reasons": reasons,
            "files_moved_or_removed": moved,
        })

    for stage_dir in sorted(p for p in figure_root.iterdir() if p.is_dir()):
        write_figure_index(stage_dir)
    write_figure_index(figure_root)
    report = {
        "status": "filtered",
        "output_root": str(output_root),
        "figure_root": str(figure_root),
        "quarantine": str(output_root / "07_figure_quarantine") if quarantine else None,
        "rejected_count": len(rejected),
        "rejected": rejected,
        "rules": {
            "single_category_violin": "reject a violin panel with fewer than two displayed categories",
            "empty_condition_panel": "reject a composite figure containing an empty violin or feature panel",
            "all_constant_violin": "reject a multi-category violin with no density-bearing category",
            "no_positive_target_expression": "reject a target feature panel with zero positive cells",
            "no_significant_genes": "reject a volcano with no significant up/down genes",
        },
    }
    (output_root / "99_logs" / "figure_quality_report.json").parent.mkdir(parents=True, exist_ok=True)
    (output_root / "99_logs" / "figure_quality_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return report


def assess_dataset(output_root: Path) -> dict:
    """Assess the whole project; one bad audited figure rejects the project."""
    output_root = output_root.resolve()
    figure_root = output_root / "04_figures"
    sidecars = sorted(figure_root.rglob("*.plot.json")) if figure_root.exists() else []
    rejected: list[dict] = []
    for sidecar in sidecars:
        reasons = quality_reasons(_audits(sidecar))
        if reasons:
            rejected.append({
                "figure": str(sidecar.relative_to(output_root)),
                "reasons": reasons,
            })
    reasons = []
    if not figure_root.exists():
        reasons.append("no_figure_directory")
    elif not sidecars:
        reasons.append("no_audited_figures")
    if rejected:
        reasons.append("low_information_figure_present")
    return {
        "dataset": output_root.name,
        "output_root": str(output_root),
        "status": "pass" if not reasons else "reject",
        "reasons": sorted(set(reasons)),
        "audited_figure_count": len(sidecars),
        "bad_figures": rejected,
    }


def quarantine_dataset(output_root: Path, quarantine_root: Path, assessment: dict) -> dict:
    """Move the complete project directory, preserving its internal files."""
    output_root = output_root.resolve()
    quarantine_root = quarantine_root.resolve()
    quarantine_root.mkdir(parents=True, exist_ok=True)
    destination = quarantine_root / f"{output_root.parent.name}__{output_root.name}"
    if destination.exists():
        shutil.rmtree(destination)
    shutil.move(str(output_root), str(destination))
    assessment = dict(assessment)
    assessment["quarantine_destination"] = str(destination)
    (destination / "figure_quality_dataset_report.json").write_text(
        json.dumps(assessment, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return assessment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--no-quarantine", action="store_true", help="delete rejected render bundles")
    parser.add_argument("--dataset-level", action="store_true", help="assess the entire dataset as one deliverable")
    parser.add_argument("--quarantine-dataset-to", type=Path, help="move a rejected dataset directory here")
    args = parser.parse_args()
    if args.dataset_level:
        assessment = assess_dataset(args.output_root)
        if assessment["status"] == "reject" and args.quarantine_dataset_to:
            assessment = quarantine_dataset(args.output_root, args.quarantine_dataset_to, assessment)
        report_path = Path(assessment.get("quarantine_destination", args.output_root)) / "99_logs" / "figure_quality_dataset_report.json"
        if assessment.get("quarantine_destination"):
            report_path = Path(assessment["quarantine_destination"]) / "figure_quality_dataset_report.json"
        else:
            report_path = args.output_root / "99_logs" / "figure_quality_dataset_report.json"
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(json.dumps(assessment, indent=2, ensure_ascii=False), encoding="utf-8")
        print(json.dumps(assessment, ensure_ascii=False))
        return 0 if assessment["status"] == "pass" else 4
    report = filter_output(args.output_root, quarantine=not args.no_quarantine)
    print(json.dumps({"status": report["status"], "rejected_count": len(report.get("rejected", [])),
                      "report": str((args.output_root / "99_logs" / "figure_quality_report.json").resolve())}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
