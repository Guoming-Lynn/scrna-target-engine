# Source-stage figure inventory and current implementation

Original snapshots are under `hadha_source/`; do not execute them as dependencies.
Reusable code: `../scripts/figure_recipes.py`. Current style: `plotting_style.md`.

| Original source | Original output families | Current handling |
|---|---|---|
| run_step01.py.txt, 244–367 | raw counts, five QC violins, pre-filter genes/counts MT scatter | Phase01: centered counts, white-box ranked violins, pre scatter |
| Same, 524–542 | post-filter MT scatter | Restored with shared pre/post axes and color limits |
| Same, 617–630 | observed/simulated doublet histogram per sample | Density grid, auto/effective threshold, distribution CSVs |
| Same, 653–743 | final counts, QC/final doublet-score scatter | Both restored |
| run_step02_round2.py.txt | metadata/QC/resolution UMAP, QC boxes (filename says violin), marker features/dotplot, target feature/violin | Metadata/resolution UMAP, marker dotplot, shared QC/expression suite; marker feature grid remains a reference candidate |
| run_step03_round1 through round3 | metadata/QC/resolution plots, module UMAP, composition, target feature/split violin, QC boxes | Subtype UMAP, QC/expression suite, composition, configured marker dotplot; full module/metadata grids are not yet wired |
| update_all_requested_figures.py.txt, 159–242 | feature plus split inner-box violin | Adapted to separate case/control feature + ranked white-box violin pairs |
| sensitivity_notebook_code.py.txt, 1146–1602 and 2274–2933 | density, jitter, white box, ranked variants, reference stars | ranked_violin: full-data density/box, median ranking, 100 sampled dots. No automatically inferred reference stars |
| plot_round3_diagnostics.py.txt | diagnostic marker violin/box grid and feature/dotplots | Optional reference for requested contamination markers; never an automatic deletion rule |
| update_round3_figures_and_naming.py.txt | final named UMAP, diagnostic dots, composition | Human review publishes named UMAP/marker/composition/QC/expression into stage directory |
| Sensitivity notebook volcano/enrichment section | target-excluded volcano and enrichment suite | Current stage 04: red/blue target-excluded volcano, ORA summaries, GSEA NES; curve/concordance remain follow-up work |

## Current naming

Use output_naming.md and figure_output.py for current stage/name mapping. The following semantic families are retained.

## QC families

`cell_counts_raw`, `pre_QC_METRIC`, `post_QC_METRIC`, `qc_pre_counts_vs_genes`,
`qc_post_counts_vs_genes`, `cell_counts_clean`; under scrublet:
`scrublet_score_distributions`, `scrublet_before_after_counts_vs_genes`.
When Scrublet is skipped, do not fabricate scores or simulated distributions.

## Shared stage 2 / stage 3 / final-review names

- `qc_metrics_by_cluster`: metrics in rows, white boxes, descending medians.
- `target_GENE_umap`: gray zeros; positive-cell p1/p99 expression scale.
- `target_GENE_violin`: all cells by cluster/final name.
- `target_GENE_by_group_violin`: case above control, independently ranked.
- `target_GENE_case_LABEL_feature_violin`: case UMAP left, violin right.
- `target_GENE_control_LABEL_feature_violin`: control pair, shared scales.

All are PNG + PDF. Sidecars record sampling and color saturation. Missing genes
skip target panels rather than fabricating all-zero expression. Preliminary plots
are archived in the review round when final named plots replace stage outputs.
