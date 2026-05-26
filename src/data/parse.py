"""Parse dot-bracket structure strings into base-pair sets.

Supports multiple bracket types as independent LIFO stacks:
    ()  []  {}  <>  and letter pairs Aa .. Zz (uppercase = open, lowercase = close)

The parser is intentionally pure: no minimum-loop-length filtering, no sequence
awareness. Filtering belongs in clean.py.
"""

from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Bracket alphabet
# ---------------------------------------------------------------------------

# Map each opening bracket to its closing counterpart.
_OPEN_TO_CLOSE: dict[str, str] = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
    **{c: c.lower() for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"},
}

_CLOSE_TO_OPEN: dict[str, str] = {v: k for k, v in _OPEN_TO_CLOSE.items()}

_OPEN_CHARS: frozenset[str] = frozenset(_OPEN_TO_CLOSE)
_CLOSE_CHARS: frozenset[str] = frozenset(_CLOSE_TO_OPEN)


# ---------------------------------------------------------------------------
# Record type shared across the data sub-package
# ---------------------------------------------------------------------------


@dataclass
class RNARecord:
    """A single RNA sequence/structure pair from the dataset.

    Attributes
    ----------
    name :
        Identifier (e.g. PDB ID or dataset label).
    sequence :
        RNA sequence string (expected alphabet: A, C, G, U).
    structure :
        Dot-bracket structure string, same length as ``sequence``.
    """

    name: str
    sequence: str
    structure: str


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def parse_structure(s: str) -> set[tuple[int, int]]:
    """Parse a dot-bracket string into a set of 0-indexed base pairs (i, j), i < j.

    Each bracket type is matched on its own independent LIFO stack, so
    pseudoknots encoded with different bracket types (e.g. ``()`` and ``[]``)
    are handled correctly.

    Parameters
    ----------
    s :
        Dot-bracket structure string.  Valid characters: ``.``, ``(``, ``)``,
        ``[``, ``]``, ``{``, ``}``, ``<``, ``>``, and the 26 letter pairs
        ``Aa``–``Zz``.

    Returns
    -------
    set of (int, int)
        Unordered set of base pairs; each pair has the smaller index first.

    Raises
    ------
    ValueError
        If the string contains unmatched brackets or unrecognised characters.
    """
    stacks: dict[str, list[int]] = {ch: [] for ch in _OPEN_CHARS}
    pairs: set[tuple[int, int]] = set()

    for pos, ch in enumerate(s):
        if ch == ".":
            continue
        if ch in _OPEN_CHARS:
            stacks[ch].append(pos)
        elif ch in _CLOSE_CHARS:
            open_ch = _CLOSE_TO_OPEN[ch]
            if not stacks[open_ch]:
                raise ValueError(
                    f"Unmatched closing '{ch}' at position {pos} in: {s!r}"
                )
            open_pos = stacks[open_ch].pop()
            pairs.add((open_pos, pos))
        else:
            raise ValueError(
                f"Unrecognised character '{ch}' at position {pos} in: {s!r}"
            )

    for open_ch, stack in stacks.items():
        if stack:
            raise ValueError(
                f"Unmatched opening '{open_ch}' at positions {stack} in: {s!r}"
            )

    return pairs
