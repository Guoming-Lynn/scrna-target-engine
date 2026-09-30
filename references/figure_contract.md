# Rapid figure contract

Read [figure_stage_reference.md](figure_stage_reference.md) for source mapping,
implemented outputs and remaining reference candidates. Read
[plotting_style.md](plotting_style.md) for executable recipes and exact styles.

Required for this QC / global / target-cell revision:

1. QC pre/post distributions and matched genes-versus-counts scatters.
2. Scrublet observed/simulated densities, threshold lines, final counts and
   matched post-QC/final doublet-score scatters (unless explicitly skipped).
3. Cluster/final-name QC distributions: white boxes, 100 points, median ranking.
4. Target UMAP, overall violin, case-over-control violin, separate condition
   feature/violin pairs with shared coordinates and color scale.
5. Final named UMAP and marker/composition panels regenerated after human review.

An original reference inventory is not an implementation guarantee. Pseudotime
is not part of the default run. Human decisions control deletion and names.
Plotting never changes counts, expression, annotations or statistical tables.
