# Output naming

Delivery uses four consecutive stage directories under both 03_tables and
04_figures: 01_qc, 02_global, 03_target, 04_deg. Scrublet is a subfolder of 01_qc.
Figures use descriptive names rather than inherited project numbers. Figure
order is not a biological claim; optional outputs need not leave numerical gaps.

| Previous name | Current name |
|---|---|
| 25_UMAP_legend | umap_cell_labels |
| 25_UMAP_on_data | umap_cell_labels_on_map |
| qc_metrics_violin | qc_metrics_by_cluster |
| 27_group_composition | composition_cell_types_by_condition |
| 16_marker_dotplot | markers_top_dotplot |
| umap_HADHA | target_HADHA_umap |
| violin_HADHA | target_HADHA_violin |
| violin_HADHA_by_group | target_HADHA_by_group_violin |
| 14_HADHA_case_AC_feature_and_violin | target_HADHA_case_AC_feature_violin |
| 14_HADHA_control_PA_feature_and_violin | target_HADHA_control_PA_feature_violin |
| final_clean_cell_counts | cell_counts_clean |
| post_QC_scatter | qc_post_counts_vs_genes |
| scrublet_histograms_grid | scrublet_score_distributions |
| 05_Volcano_AC_Case_HADHA_high_vs_low | volcano_AC_Case_HADHA_high_vs_low |
| 08_Enrichment_Case_UP_barplot | ora_case_up_pathways |
| 10_GSEA_top_pathways_NES_barplot | gsea_case_pathway_nes |

The mapping is implemented by figure_output.figure_path, not by renaming user
files after analysis. Direct Scanpy exports must also use this mapping/index.
A new figure should use a descriptive basename directly, with no arbitrary
prefix. The registry remains idempotent so repeated calls never add prefixes.

Every folder has an HTML gallery and CSV index, based on actual PNG/PDF files.
.plot.json records data/display details; .layout.json records tick-overlap
warnings; .caption.txt explains boxes/scales where relevant. These are sidecars,
not additional figures. A volcano .plot_data.csv preserves pre-display values.

Legacy executable and checkpoint names remain stable to allow old checkpoints
to be used. Old analysis directories are not silently renamed. Replot into a
fresh directory for a clean current-naming delivery. Original source snapshots
under references/hadha_source keep their historical names for traceability.

The root and stage figure_index.csv files recursively include nested PNGs with relative paths (including 01_qc/scrublet); rebuild after publishing reviewed figures. Subtype marker tables use subtype_markers_all.csv.
