"""Filtering helpers for RNA sequence/structure records.

Each function takes a list of RNARecord objects, applies one filter criterion,
and returns the surviving records plus a count of how many were removed.
Filters are composable: pipe the output of one into the next.

Canonical RNA bases: A, C, G, U.
"""

from __future__ import annotations

from src.data.parse import RNARecord

# Default alphabet accepted as canonical.
CANONICAL_BASES: frozenset[str] = frozenset("ACGU")


def filter_length_mismatch(
    records: list[RNARecord],
) -> tuple[list[RNARecord], int]:
    """Remove records where the sequence and structure lengths differ.

    Parameters
    ----------
    records :
        Input record list.

    Returns
    -------
    (kept, n_removed)
        ``kept`` — records whose sequence and structure have equal length.
        ``n_removed`` — number of records dropped.
    """
    kept = [r for r in records if len(r.sequence) == len(r.structure)]
    return kept, len(records) - len(kept)


def filter_noncanonical(
    records: list[RNARecord],
    alphabet: frozenset[str] = CANONICAL_BASES,
) -> tuple[list[RNARecord], int]:
    """Remove records that contain bases outside the allowed alphabet.

    The check is case-sensitive; convert sequences to upper-case before
    calling if the source data uses mixed case.

    Parameters
    ----------
    records :
        Input record list.
    alphabet :
        Set of accepted single-character base labels.  Defaults to
        ``{'A', 'C', 'G', 'U'}``.

    Returns
    -------
    (kept, n_removed)
        ``kept`` — records whose entire sequence is drawn from ``alphabet``.
        ``n_removed`` — number of records dropped.
    """
    kept = [r for r in records if all(b in alphabet for b in r.sequence)]
    return kept, len(records) - len(kept)


def filter_length(
    records: list[RNARecord],
    min_len: int | None = None,
    max_len: int | None = None,
) -> tuple[list[RNARecord], int]:
    """Remove records whose sequence length is outside [min_len, max_len].

    Either bound may be ``None`` to skip that check.

    Parameters
    ----------
    records :
        Input record list.
    min_len :
        Minimum allowed sequence length (inclusive).  ``None`` means no lower bound.
    max_len :
        Maximum allowed sequence length (inclusive).  ``None`` means no upper bound.

    Returns
    -------
    (kept, n_removed)
        ``kept`` — records within the length bounds.
        ``n_removed`` — number of records dropped.
    """
    def _in_bounds(r: RNARecord) -> bool:
        n = len(r.sequence)
        if min_len is not None and n < min_len:
            return False
        if max_len is not None and n > max_len:
            return False
        return True

    kept = [r for r in records if _in_bounds(r)]
    return kept, len(records) - len(kept)
