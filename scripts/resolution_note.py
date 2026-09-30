"""Label for the quicklook resolution heuristic."""

RESOLUTION_NOTE = "分辨率由启发式自动选择"
SELECTION_METHOD = "heuristic_first_10_to_15_clusters"


def choose_resolution(counts: list[tuple[float, int]]) -> float:
    """Pick the first resolution with 10 to 15 clusters, otherwise the first one."""
    if not counts:
        raise ValueError("resolution list is empty")
    for resolution, n_clusters in counts:
        if 10 <= n_clusters <= 15:
            return resolution
    return counts[0][0]
