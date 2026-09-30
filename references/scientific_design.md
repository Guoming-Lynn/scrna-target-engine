# Scientific plotting principles adapted for rapid scRNA reports

This package incorporates selected rules from the installed nature-figure,
scipilot-figure-skill and nature-statistics skills. It is self-contained and does
not require those skills, their local paths, SciencePlots or an image API at run
time. User-approved styles and exploratory scope take priority over generic
journal defaults.

**Agent requirement:** portability is not permission to ignore installed skills.
Before plotting, follow SKILL.md's mandatory discovery/read/apply/visual-review
instruction. When relevant plotting skills are installed, consult their actual
current requirements in addition to this distilled reference and record their
application. The user's explicit figure rules take priority.

| Source principle | Pipeline implementation | Boundary |
|---|---|---|
| nature-figure: semantic color consistency | UMAP category labels bind violin colors permanently; blue/orange only for condition encoding | Do not replace the approved palette with another journal theme |
| nature-figure: shared comparison scales | Same expression percentile limits and embedding coordinates across case/control | Percentiles affect color only, never data |
| SciPilot: distribution before summary | Full-data white-box violins plus 100 sampled points; n shown | n<=100 groups are omitted from that panel only; eligible constant groups omit KDE |
| SciPilot: full labels and readable density | Wrap long names, separate condition panels, rasterize dense points | Do not crop names or remove cells to improve appearance |
| SciPilot: render and inspect | Shared exporter writes tick-overlap audit; agent inspects actual PNGs | Automated audit does not certify every title, legend or statistical claim |
| nature-statistics: define plotted n/uncertainty | n denotes cells, box=IQR, whiskers=1.5 IQR, median line; sidecar caption states sampling and absence of tests | Do not reinterpret sampled dots as independent donors or add invented stars |
| nature-figure/SciPilot: disclose axis transformations | Volcano saturation uses visible boundary triangles and plot-data CSV | Full DE tables are unchanged and the configured target is plot-excluded only |
| nature-statistics: distinguish score from significance | GSEA shows signed NES, and hatching identifies FDR>=0.05 | Positive/negative NES does not mean case/control identity |

## Figure questions and data fields

- QC: what cell-quality distributions changed? Cell-level UMI, detected genes,
  MT/ribo/HB percentages; pre/post comparison uses matched axes where applicable.
- UMAP: where are annotated cells/target expression localized? Use the existing
  checkpoint embedding, categorical color map or numeric expression colorbar.
- Violin: what is the distribution within each displayed category? Rank by
  median while keeping label identity and color; do not infer a biological order.
- Composition: which denominator defines 100%? The filename distinguishes
  cell types within condition from conditions within cell type.
- Volcano: which genes satisfy the declared FDR AND effect-size thresholds?
  Show signs, sample contrast in title, and disclose saturated display values.
- ORA/GSEA: which pathways were tested/enriched? Full pathway names remain
  accessible; the primary results table remains the source of exact q values.

## Visual review before delivery

Inspect UMAP, a long-label violin, composition, QC scatter and a volcano at the
intended group-meeting size. Check category-color agreement, text collisions,
colorbars, missing categories, full name wrapping, panel order and shared scales.
Read .layout.json warnings and resolve visible overlap. Do not copy a journal's
small print font size into a slide by default. PNG 300 dpi + PDF is the rapid
report default; optional SVG is available if a downstream editor needs it.

No promises of 'publication quality' follow automatically from a theme or a
validator. Preserve the manual cluster-removal/final-naming gates. This revision
does not silently introduce donor models, new cell filters or network analysis.
