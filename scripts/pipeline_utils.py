"""Small shared helpers for robust, single-cell exploratory runs."""
from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd


def resolve_target_gene(var_names: Iterable[str], requested: str) -> Optional[str]:
    """Resolve an exact gene name, then an unambiguous case-insensitive match."""
    names = pd.Index([str(x) for x in var_names])
    requested = str(requested or "").strip()
    if not requested:
        return None
    if requested in names:
        return requested
    matches = names[names.str.casefold() == requested.casefold()]
    return str(matches[0]) if len(matches) == 1 else None


def is_target_name(value: object, target: str) -> bool:
    return str(value).casefold() == str(target).casefold()
