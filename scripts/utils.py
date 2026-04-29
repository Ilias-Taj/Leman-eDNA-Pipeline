#!/usr/bin/env python3
"""Shared utilities for the eDNA pipeline scripts.

Centralises common functions used across multiple pipeline steps:
- Tool discovery (vsearch, filtlong)
- FASTQ I/O helpers
- Directory traversal helpers
- Timing context manager
"""

import shutil
import sys
import time
from pathlib import Path
from typing import Generator, Optional, Tuple

# ── Constants ────────────────────────────────────────────────────────────────

# Directories created by the pipeline that should be skipped when iterating
# over barcode sample folders.
SKIP_DIRS = frozenset({
    "logs", "validation", "merged", "temp_clustering",
    "taxonomy", "taxonomy_summary", "blast_results",
})


# ── Tool discovery ───────────────────────────────────────────────────────────

def find_tool(name: str, extra_candidates: Optional[list] = None) -> str:
    """Locate a command-line tool by checking common paths then PATH.

    Args:
        name: Tool name (e.g. "vsearch", "filtlong").
        extra_candidates: Additional absolute paths to check first.

    Returns:
        Absolute path to the tool.

    Raises:
        SystemExit: If the tool is not found anywhere.
    """
    candidates = [
        Path(f"./env/bin/{name}"),
        Path(f"/opt/homebrew/bin/{name}"),
        Path(f"/usr/local/bin/{name}"),
    ]
    if extra_candidates:
        candidates = [Path(p) for p in extra_candidates] + candidates

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    found = shutil.which(name)
    if found:
        return found

    print(f"ERROR: '{name}' not found in PATH or common locations.", file=sys.stderr)
    print(f"  Install with: conda install -c bioconda {name}", file=sys.stderr)
    sys.exit(1)


def find_vsearch() -> str:
    """Shorthand for find_tool('vsearch')."""
    return find_tool("vsearch")


# ── FASTQ helpers ────────────────────────────────────────────────────────────

def iter_fastq(fileobj) -> Generator[Tuple[str, str, str, str], None, None]:
    """Yield (header, seq, plus, qual) tuples from an open FASTQ file handle.

    Each element is a raw line *with* its trailing newline so it can be
    written directly to an output file without modification.

    Args:
        fileobj: An open text-mode file handle (plain or gzip.open('rt')).

    Yields:
        Tuple of (header_line, sequence_line, plus_line, quality_line).
    """
    buf = []
    for line in fileobj:
        buf.append(line)
        if len(buf) == 4:
            yield buf[0], buf[1], buf[2], buf[3]
            buf = []


# ── Directory helpers ────────────────────────────────────────────────────────

def iter_barcode_dirs(input_dir: Path):
    """Yield sorted barcode directories, skipping pipeline output folders.

    Args:
        input_dir: Root run directory (e.g. out/Water_eDNA_18S_COI_14_01_26).

    Yields:
        Path objects for each sample/barcode subdirectory.
    """
    for d in sorted(input_dir.iterdir()):
        if d.is_dir() and d.name not in SKIP_DIRS:
            yield d


# ── Timing ───────────────────────────────────────────────────────────────────

class timer:
    """Context manager that times a block and prints elapsed seconds.

    Usage::

        with timer("Clustering 18S") as t:
            do_work()
        print(f"took {t.elapsed:.1f}s")
    """

    def __init__(self, label: str):
        self.label = label
        self.elapsed = 0.0

    def __enter__(self):
        self._start = time.time()
        return self

    def __exit__(self, *exc):
        self.elapsed = time.time() - self._start
        print(f"[TIME] {self.label}: {self.elapsed:.1f}s")
        return False
