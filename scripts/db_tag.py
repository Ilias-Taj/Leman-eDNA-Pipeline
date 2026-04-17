"""Helper to derive a short, filesystem-safe DB tag from one or more .udb paths.

Used by scripts 5 and 7 to place DB-specific outputs into
`taxonomy/{tag}/` and `taxonomy_summary/{tag}/` subfolders so multiple
databases can coexist without re-running the whole pipeline.

Examples:
    derive_tag(["refs/silva_18s_v123.udb"])                       -> "silva"
    derive_tag(["refs/midori2_COI.udb"])                          -> "midori2"
    derive_tag(["refs/eKOI_COI.udb"])                             -> "ekoi"
    derive_tag(["refs/silva_18s_v123.udb", "refs/eKOI_COI.udb"])  -> "silva-ekoi"
    derive_tag(["refs/eKOI_COI.udb", "refs/eKOI_COI.udb"])        -> "ekoi"
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Optional

_STRIP_PATTERNS = [
    re.compile(r"_v\d+$", re.IGNORECASE),
    re.compile(r"_18s$", re.IGNORECASE),
    re.compile(r"_coi$", re.IGNORECASE),
    re.compile(r"_jedi$", re.IGNORECASE),
    re.compile(r"_rrna$", re.IGNORECASE),
    re.compile(r"_ssu$", re.IGNORECASE),
    re.compile(r"_nr$", re.IGNORECASE),
]


def label_from_path(db_path: str) -> str:
    """Return short label from a .udb filename (lowercased, marker/version stripped)."""
    stem = Path(db_path).stem.lower()
    changed = True
    while changed:
        changed = False
        for pat in _STRIP_PATTERNS:
            new_stem = pat.sub("", stem)
            if new_stem != stem:
                stem = new_stem
                changed = True
    return stem or "db"


def derive_tag(db_paths: Iterable[Optional[str]]) -> str:
    """Derive a combined tag from a sequence of DB paths.

    - None/empty entries are ignored.
    - Duplicate labels are collapsed while preserving order.
    - Joined with '-'.
    - Returns 'default' if no DBs provided.
    """
    labels: list[str] = []
    for p in db_paths:
        if not p:
            continue
        lab = label_from_path(p)
        if lab and lab not in labels:
            labels.append(lab)
    return "-".join(labels) if labels else "default"
