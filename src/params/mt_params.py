"""Parse the ViennaRNA-format MT parameter file (rna_DirksPierce09.par).

Reads the stacking, mismatch, loop, dangle and other Mathews-Turner sections.
These parameters are fixed inputs to the energy function; they are not estimated.

Unit convention: all integer values in the file are kcal/mol x 100.  The
``MTParams`` object stores them as-is (integers); callers should divide by 100
when they need kcal/mol floats.

Pair-type indexing (ViennaRNA convention used throughout this file):
    0 = NN (undefined / no pair)
    1 = CG
    2 = GC
    3 = GU
    4 = UG
    5 = AU
    6 = UA

The ``stack`` table is indexed [closing_pair][opening_pair], where "closing"
is the outer pair (i, j) and "opening" is the inner pair (i+1, j-1).
Both indices run from 0 to 6 (NN included).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Sentinel for "not stacked" / infinity entries in the .par file.
_NST = 100_000  # large integer; never a real energy (kcal/mol x 100)

# Pair-type labels in the order they appear in the file (index 0 is NN).
PAIR_TYPES: tuple[str, ...] = ("NN", "CG", "GC", "GU", "UG", "AU", "UA")
PAIR_INDEX: dict[str, int] = {p: i for i, p in enumerate(PAIR_TYPES)}

# Base labels used in mismatch / dangle tables.
BASE_TYPES: tuple[str, ...] = ("N", "A", "C", "G", "U")
BASE_INDEX: dict[str, int] = {b: i for i, b in enumerate(BASE_TYPES)}


def _parse_int(token: str) -> int:
    """Convert a single token from the .par file to an integer.

    ``NST`` and ``INF`` are mapped to the internal sentinel ``_NST``.
    """
    if token in ("NST", "INF"):
        return _NST
    return int(token)


def _strip_comments(line: str) -> str:
    """Remove inline C-style comments (/* ... */) and trailing whitespace."""
    line = re.sub(r"/\*.*?\*/", "", line)
    return line.strip()


def _tokenize(line: str) -> list[str]:
    """Return the numeric/sentinel tokens on a single (comment-stripped) line."""
    return _strip_comments(line).split()


@dataclass
class MTParams:
    """Mathews-Turner parameters loaded from a ViennaRNA-format .par file.

    Only the ``stack`` section is fully parsed on construction; other sections
    are stored as raw text blocks for future expansion.

    Attributes
    ----------
    stack : np.ndarray, shape (7, 7), dtype int32
        Stacking free energies in kcal/mol x 100.  Indexed
        ``stack[closing_pair][opening_pair]`` using the PAIR_INDEX mapping.
        Undefined entries (NST) are stored as ``_NST`` (100_000).
    raw_sections : dict[str, list[str]]
        All other section names mapped to their raw lines, for future parsing.
    """

    stack: np.ndarray  # shape (7, 7), kcal/mol x 100
    raw_sections: dict[str, list[str]] = field(default_factory=dict)


def load_mt_params(path: str | Path) -> MTParams:
    """Parse a ViennaRNA-format parameter file and return an ``MTParams`` object.

    Only the ``stack`` section is fully parsed; all other sections are kept as
    raw line lists so they can be parsed on demand later.

    Parameters
    ----------
    path :
        Path to the ``.par`` file (e.g. ``params/rna_DirksPierce09.par``).

    Returns
    -------
    MTParams
        Parsed MT parameters.  ``stack`` is a 7×7 integer array.

    Raises
    ------
    ValueError
        If the ``stack`` section is not found or has unexpected dimensions.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    # Split file into named sections.
    sections: dict[str, list[str]] = {}
    current_section: str | None = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("##"):
            # File header — ignore.
            continue
        if stripped.startswith("# "):
            current_section = stripped[2:].strip()
            sections[current_section] = []
        elif current_section is not None:
            sections[current_section].append(line)

    if "stack" not in sections:
        raise ValueError("'stack' section not found in parameter file.")

    stack = _parse_stack(sections.pop("stack"))
    return MTParams(stack=stack, raw_sections=sections)


def _parse_stack(lines: list[str]) -> np.ndarray:
    """Parse the ``stack`` section into a 7×7 numpy array.

    The file lists 7 rows (CG, GC, GU, UG, AU, UA, NN) with 7 values each.
    Index 0 is reserved for NN; the file rows are mapped to indices 1-6 then 0.

    Returns
    -------
    np.ndarray, shape (7, 7), dtype int32
        ``arr[closing_pair_index][opening_pair_index]``.
    """
    # File row order: CG GC GU UG AU UA NN  → PAIR_INDEX 1 2 3 4 5 6 0
    file_row_order = [1, 2, 3, 4, 5, 6, 0]

    rows: list[list[int]] = []
    for line in lines:
        tokens = _tokenize(line)
        if not tokens:
            continue
        rows.append([_parse_int(t) for t in tokens])

    if len(rows) != 7:
        raise ValueError(f"Expected 7 rows in 'stack' section, got {len(rows)}.")
    for i, row in enumerate(rows):
        if len(row) != 7:
            raise ValueError(
                f"Expected 7 values in stack row {i}, got {len(row)}: {row}"
            )

    arr = np.full((7, 7), _NST, dtype=np.int32)
    for file_idx, pair_idx in enumerate(file_row_order):
        for file_col, col_pair_idx in enumerate(file_row_order):
            arr[pair_idx, col_pair_idx] = rows[file_idx][file_col]

    return arr


def stack_energy(mt: MTParams, closing: str, opening: str) -> int:
    """Look up the stacking energy for a base-pair stack in kcal/mol x 100.

    Parameters
    ----------
    mt :
        Loaded MT parameters.
    closing :
        The outer (closing) pair type, e.g. ``"CG"``.
    opening :
        The inner (opening) pair type, e.g. ``"AU"``.

    Returns
    -------
    int
        Stacking energy in kcal/mol x 100.  Returns ``_NST`` (100_000) if
        either pair type is unknown or the combination is undefined (NST).
    """
    ci = PAIR_INDEX.get(closing, 0)
    oi = PAIR_INDEX.get(opening, 0)
    return int(mt.stack[ci, oi])
