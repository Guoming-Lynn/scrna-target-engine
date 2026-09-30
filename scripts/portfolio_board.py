#!/usr/bin/env python3
"""Build a dense, stage-organized portfolio board from pipeline PNG outputs.

The board is an overview artifact. It preserves each source PNG unchanged,
uses a fixed physical canvas and records every included source in a manifest.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import textwrap
from pathlib import Path

from PIL import Image
from PIL import ImageDraw, ImageFont


BG = "#EEF1F4"
HEADER = "#D9DEE3"
INK = "#19232D"
WHITE = "#FFFFFF"
DEFAULT_DPI = 180
CANVAS_WIDTH_IN = 1960 / DEFAULT_DPI

STAGE_SPECS = [
    ("01_qc", "01  QC / SCRUBLET", [
        "cell_counts_clean", "cell_counts_raw", "qc_post_counts_vs_genes",
        "qc_post_log10_total_counts", "qc_post_n_genes", "qc_post_pct_hb",
        "qc_post_pct_mt", "qc_post_pct_ribo", "qc_pre_counts_vs_genes",
        "qc_pre_log10_total_counts", "qc_pre_n_genes", "qc_pre_pct_hb",
        "qc_pre_pct_mt", "qc_pre_pct_ribo", "scrublet_before_after_counts_vs_genes",
        "scrublet_score_distributions",
    ]),
    ("02_global", "02  GLOBAL CLUSTERING", [
        "cluster_to_final_name", "composition_cell_types_by_condition", "markers_top_dotplot",
        "qc_metrics_by_cluster", "target_.*_by_group_violin", "target_.*_case_.*_feature_violin",
        "target_.*_control_.*_feature_violin", "target_.*_small_n_summary", "target_.*_umap",
        "target_.*_violin", "umap_by_condition", "umap_by_sample", "umap_cell_labels_on_map",
        "umap_cell_labels", "umap_module_scores", "umap_preliminary_clusters",
    ]),
    ("03_target", "03  TARGET CELL ANALYSIS", [
        "cluster_to_final_name", "composition_cell_types_by_condition", "markers_top_dotplot",
        "qc_metrics_by_cluster", "target_.*_by_group_violin", "target_.*_case_.*_feature_violin",
        "target_.*_control_.*_feature_violin", "target_.*_small_n_summary", "target_.*_umap",
        "target_.*_violin", "umap_by_condition", "umap_by_sample", "umap_cell_labels_on_map",
        "umap_cell_labels", "umap_module_scores", "umap_preliminary_clusters",
    ]),
    ("04_deg", "04  DEG / ENRICHMENT", [
        "composition_target_expression_strata", "gsea_.*_pathway_nes", "ora_.*_down_pathways",
        "ora_.*_up_pathways", "volcano_.*_case_.*_high_vs_low", "volcano_.*_all_cells_.*_high_vs_low",
        "volcano_.*_control_.*_high_vs_low",
    ]),
]


def _natural_key(path: Path) -> tuple:
    return tuple(int(part) if part.isdigit() else part.casefold() for part in re.split(r"(\d+)", path.stem))


def _ordered_pngs(stage_dir: Path, patterns: list[str]) -> list[Path]:
    candidates = sorted(stage_dir.rglob("*.png"), key=_natural_key)
    by_stem = {p.stem: p for p in candidates}
    selected: list[Path] = []
    used: set[Path] = set()
    for pattern in patterns:
        rx = re.compile(rf"^{pattern}$", re.IGNORECASE)
        matches = sorted((p for p in candidates if rx.match(p.stem)), key=_natural_key)
        for path in matches:
            if path not in used:
                selected.append(path)
                used.add(path)
    for path in candidates:
        if path not in used:
            selected.append(path)
    return selected


def _label(path: Path) -> str:
    text = path.stem.replace("_", " ")
    return "\n".join(textwrap.wrap(text, width=28, break_long_words=False))


def _read_thumbnail(path: Path, max_size: tuple[int, int] = (900, 620)):
    with Image.open(path) as source:
        image = source.convert("RGB")
        image.thumbnail(max_size, Image.Resampling.LANCZOS)
        return image.copy()


def build_board(output_root: Path, output_dir: Path | None = None, target_gene: str | None = None) -> dict:
    figure_root = output_root / "04_figures"
    output_dir = output_dir or (output_root / "06_portfolio")
    output_dir.mkdir(parents=True, exist_ok=True)

    config_target = target_gene
    for candidate in (output_root / "99_logs" / "resolved_config.json", output_root / "config.json"):
        if config_target:
            break
        if candidate.exists():
            try:
                config_target = json.loads(candidate.read_text(encoding="utf-8-sig")).get("TARGET_GENE")
            except Exception:
                pass
    target_gene = str(config_target or "scRNA").strip()

    sections: list[tuple[str, list[Path]]] = []
    for stage, title, patterns in STAGE_SPECS:
        stage_dir = figure_root / stage
        if stage_dir.exists():
            paths = _ordered_pngs(stage_dir, patterns)
            if paths:
                sections.append((title, paths))
    if not sections:
        raise FileNotFoundError(f"No PNG figures found under {figure_root}")

    columns = 6
    tile_h = 1.27
    image_h = 1.05
    label_h = tile_h - image_h
    section_header_h = 0.38
    section_gap = 0.08
    top_h = 0.52
    bottom_h = 0.22
    total_h = top_h + bottom_h + sum(
        section_header_h + math.ceil(len(paths) / columns) * tile_h + section_gap
        for _, paths in sections
    )

    canvas_w = 1960
    px_per_in = DEFAULT_DPI
    canvas_h = int(round(total_h * px_per_in))
    canvas = Image.new("RGB", (canvas_w, canvas_h), BG)
    draw = ImageDraw.Draw(canvas)
    font_candidates = [
        r"C:\Windows\Fonts\Arial.ttf",
        r"C:\Windows\Fonts\Aptos.ttf",
        r"C:\Windows\Fonts\DejaVuSans.ttf",
    ]
    bold_candidates = [
        r"C:\Windows\Fonts\Arialbd.ttf",
        r"C:\Windows\Fonts\Aptos Display.ttf",
        r"C:\Windows\Fonts\DejaVuSans-Bold.ttf",
    ]
    def _font(candidates, size):
        for candidate in candidates:
            if Path(candidate).exists():
                return ImageFont.truetype(candidate, size=size)
        return ImageFont.load_default()
    title_font = _font(bold_candidates, 27)
    section_font = _font(bold_candidates, 16)
    label_font = _font(font_candidates, 8)
    draw.rectangle((0, 0, canvas_w, int(top_h * px_per_in)), fill=HEADER)
    draw.text((24, 21), f"{target_gene} | single-cell analysis portfolio", fill=INK, font=title_font)

    left = 20
    right = 20
    gap = 16
    tile_w = int((canvas_w - left - right - (columns - 1) * gap) / columns)
    cursor_y = int(top_h * px_per_in)
    section_header_px = int(section_header_h * px_per_in)
    tile_h_px = int(tile_h * px_per_in)
    image_h_px = int(image_h * px_per_in)
    label_h_px = tile_h_px - image_h_px
    manifest_sections = []

    for title, paths in sections:
        draw.rectangle((0, cursor_y, canvas_w, cursor_y + section_header_px), fill=BG)
        draw.text((20, cursor_y + 11), f"{title}  ({len(paths)} figures)", fill=INK, font=section_font)
        cursor_y += section_header_px + 8
        rows = math.ceil(len(paths) / columns)
        section_manifest = {"title": title, "figures": []}
        for idx, path in enumerate(paths):
            row, col = divmod(idx, columns)
            x = left + col * (tile_w + gap)
            y = cursor_y + row * tile_h_px
            draw.rectangle((x, y, x + tile_w, y + image_h_px), fill=WHITE)
            thumb = _read_thumbnail(path, (tile_w - 8, image_h_px - 8))
            tx = x + (tile_w - thumb.width) // 2
            ty = y + (image_h_px - thumb.height) // 2
            canvas.paste(thumb, (tx, ty))
            draw.multiline_text((x, y + image_h_px + 5), _label(path), fill="#33404C", font=label_font, spacing=0)
            section_manifest["figures"].append({
                "label": path.stem,
                "source": str(path.resolve()),
                "relative_source": str(path.relative_to(output_root)),
            })
        cursor_y += rows * tile_h_px + int(section_gap * px_per_in)
        manifest_sections.append(section_manifest)

    stem = f"{target_gene}_all_figures_dense_v2"
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    canvas.save(png_path, dpi=(DEFAULT_DPI, DEFAULT_DPI), optimize=True)
    canvas.save(pdf_path, "PDF", resolution=DEFAULT_DPI)

    manifest = {
        "title": f"{target_gene} | single-cell analysis portfolio",
        "output_png": str(png_path.resolve()),
        "output_pdf": str(pdf_path.resolve()),
        "dpi": DEFAULT_DPI,
        "canvas_width_px": 1960,
        "source_figure_count": sum(len(s["figures"]) for s in manifest_sections),
        "sections": manifest_sections,
        "accessibility_note": "Overview board; each panel retains its source PNG and filename label.",
    }
    (output_dir / f"{stem}.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--target-gene")
    args = parser.parse_args()
    manifest = build_board(args.output_root.resolve(), args.output_dir.resolve() if args.output_dir else None, args.target_gene)
    print(json.dumps({"status": "ok", "source_figure_count": manifest["source_figure_count"],
                      "output_png": manifest["output_png"], "output_pdf": manifest["output_pdf"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
