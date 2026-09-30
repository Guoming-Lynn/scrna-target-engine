"""Input rejection reasons shared by the loader."""


def input_rejection(
    *,
    n_cells: int,
    n_genes: int,
    min_cells: int,
    min_genes: int,
    numeric: bool,
    minimum: float,
    target_resolved: str | None,
    target_requested: str,
) -> str | None:
    if n_cells < min_cells:
        return f"Expression matrix has too few cells ({n_cells})"
    if n_genes < min_genes:
        return f"Expression matrix has too few features ({n_genes})"
    if not numeric:
        return "Expression matrix is not numeric"
    if minimum < 0:
        return "Expression matrix contains negative values"
    if target_resolved is None:
        return (
            f"Target gene '{target_requested}' was not found "
            "(exact or unambiguous case-insensitive match)"
        )
    return None
