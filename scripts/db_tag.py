"""Helper to derive per-marker DB labels for taxonomy output paths.

New layout: taxonomy/{MARKER}/{db_label}/taxonomy_{MARKER}.txt
Each marker's taxonomy is stored independently per DB, so running with
different DB combinations doesn't duplicate files.

Examples:
    marker_db_label("refs/pr2_18s_v511.udb")  -> "pr2"
    marker_db_label("refs/midori2_COI.udb")    -> "midori2"
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

# Hardcoded mapping from known .udb filenames to short labels
_DB_LABELS = {
    "silva_18s_v123":  "silva",
    "pr2_18s_v511":    "pr2",
    "midori2_coi":     "midori2",
    "ekoi_coi":        "ekoi",
    "porter_coi_v51":  "porter",
}


def label_from_path(db_path: str) -> str:
    """Return short label for a known .udb file, or the stem as fallback."""
    stem = Path(db_path).stem.lower()
    return _DB_LABELS.get(stem, stem)
