# QC / stage 2 / stage 3 plotting contract

Use `scripts/figure_recipes.py` for executable recipes. The purpose is to show
retained-cell quality, expression distribution and localization; these panels
do not establish a direction of disease effect. Python/matplotlib only.

- Condition colors: **Case #E69F00**, **Control #0072B2** (Okabe–Ito).
  Use these only when marks directly encode disease/control. A panel filtered
  to Case still uses subtype colors if its violins represent subtypes. Both
  rows of the expression violin and both condition feature/violin pairs must match the same
  subtype UMAP. Resolve colors from categorical labels plus `uns[key_colors]`,
  never from the sorted display positions; preserve an existing UMAP palette.
  Observed/simulated doublets also use blue/orange, with explicit labels.
  Multiclass cell types keep stable categorical colors. Volcano retains the
  separately requested red-up/blue-down convention and excludes the target gene.
- A category is plotted only if **n > 100 (at least 101)** finite cells are
  available in the current panel. Case and Control are evaluated separately.
  Omit its violin, box, dots and tick if n <= 100; this is not a cell filter.
  Record omitted categories/counts in .plot.json. If none qualify, show a clear
  no-eligible-categories message. The 100 jitter dots are a separate sampling cap.
- Violin body uses all finite observations, KDE cut at the observed minimum and
  maximum. Overlay an opaque white IQR box, red median (`#D62728`), 1.5-IQR whiskers.
  Only scatter points are sampled: at most **100 cells per displayed category**,
  without replacement, `RANDOM_SEED`, jitter ±0.20, alpha 0.35. No significance
  stars are inferred by these plotting helpers.
- Each panel orders categories by **median descending**, ties by mean descending
  then name. Panels can therefore have different category orders; labels are
  repeated, not shared. Color mapping stays attached to category identity.
  Constant/all-zero groups show boxes/medians and points without fake density.
- Long names wrap; do not truncate names into ambiguous prefixes. Bar positions
  and tick positions are identical integers. Sample bars are not hue-dodged.
- `target_GENE_by_group_violin`: **2 rows × 1 column**, case above control; shared y
  scale, independently ranked categories and full-data boxes in both panels.
- Condition feature/violin pairs use two files: `target_GENE_case_LABEL_feature_violin` and
  `target_GENE_control_LABEL_feature_violin`. Each contains only that condition,
  UMAP left and violin right. Both use the same existing embedding, coordinate
  bounds and expression-color range; no condition-specific re-embedding.
- Gene maps draw zero cells gray, positive cells low-to-high. Default color
  limits are p1/p99 **of positive cells in the entire checkpoint**, shared by
  case/control. `FEATURE_PERCENTILES: [1, 99]` is configurable. This saturates
  only the display scale; expression values/cells/statistics remain unchanged.
- QC pre/post genes-versus-counts scatters share axes and MT color scale.
  Scrublet before/after scatters share axes and doublet-score p1/p99 scale.
- White background, no decorative grids, left/bottom axes only, consistent
  sans-serif fonts, explicit units/group/cell counts. Deliver PNG at 300 dpi
  and editable-text PDF; rasterize dense marks, keep labels as vector text.
- Sidecar `.plot.json` records distribution sample sizes, seed, order, medians,
  missing-value exclusions and color saturation counts. It is a display audit,
  not a statistical test report.

## Figure-only revision

Use `scripts/replot_checkpoint.py --config ... --checkpoint ... --groupby
cell_type --output NEW_EMPTY_DIR` for global figures, or `--groupby subtype`
for finalized target cells. It regenerates review/marker/expression/QC panels
from the existing checkpoint without changing annotations, clustering or h5ad.
For initial unresolved checkpoints use `leiden` or `target_leiden_selected`.

Human deletion and final naming remain required. On acceptance, preliminary
stage figures move into that review round's `pre_review_stage_figures`; final
figures and naming tables are published back into normal stage directories.

Volcano labels: rank by adjusted significance within each direction, keep a spaced label column near the selected points, and order ties at saturated boundaries to avoid crossing leaders. Inspect leaders as well as text overlap in every rendered volcano. Composition plots have no decorative grid. Violin captions must say red median, matching the rendered line.
