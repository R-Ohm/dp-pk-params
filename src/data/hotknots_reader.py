"""Reader for the HotKnots S-Train / S-Test file format.

Each record is three consecutive lines:
    1. Header: > PSEUDOBASE_ID=PKB00240; ABBREVIATION=BChV ; ORGANISM=...
    2. Sequence: ACGUACGU...
    3. Structure: (((...)))...

The reader parses the header into a name (PSEUDOBASE_ID) and stores the
abbreviation as secondary metadata, then applies the standard cleaning filters.
"""

from __future__ import annotations

from pathlib import Path

from src.data.clean import filter_length, filter_length_mismatch, filter_noncanonical
from src.data.parse import RNARecord


def _parse_header(line: str) -> tuple[str, str]:
    """Return (name, abbreviation) from a header line.

    Handles two formats:
    - Key-value:  > PSEUDOBASE_ID=PKB00240; ABBREVIATION=BChV ; ...
    - Bare name:  > AKV-MuLV
    Keys and values are stripped of whitespace; missing keys return empty string.
    """
    body = line.lstrip(">").strip()
    fields: dict[str, str] = {}
    for part in body.split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            fields[k.strip()] = v.strip()

    if fields:
        name = fields.get("PSEUDOBASE_ID") or fields.get("SSTRAND_ID", "")
        abbrev = fields.get("ABBREVIATION", "")
    else:
        # Bare-name header: > AKV-MuLV
        name = body
        abbrev = ""
    return name, abbrev


def normalize_sequence(seq: str) -> str:
    """Normalize an RNA sequence to uppercase canonical bases.

    Converts lowercase to uppercase and substitutes T -> U (DNA notation).
    """
    return seq.upper().replace("T", "U")


def load_hotknots_file(
    path: str | Path,
    min_len: int | None = None,
    max_len: int | None = None,
) -> tuple[list[RNARecord], dict[str, int]]:
    """Parse a HotKnots-format file into a cleaned list of RNARecord.

    Applies cleaning filters in order:
        1. length mismatch  (len(sequence) != len(structure))
        2. non-canonical bases  (anything outside A, C, G, U)
        3. length bounds  (optional min/max)

    Parameters
    ----------
    path :
        Path to the .txt file (e.g. ``data/TRA_HotKnots_max200_S80P.txt``).
    min_len, max_len :
        Optional inclusive length bounds passed to ``filter_length``.

    Returns
    -------
    (records, stats)
        records — cleaned RNARecord list.
        stats   — counts at each stage: raw, after each filter, and n_removed
                  by each filter.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    raw: list[RNARecord] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith(">"):
            pk_id, abbrev = _parse_header(line)
            name = pk_id or abbrev or f"record_{len(raw)}"
            seq = normalize_sequence(lines[i + 1].strip() if i + 1 < len(lines) else "")
            struct = lines[i + 2].strip() if i + 2 < len(lines) else ""
            raw.append(RNARecord(name=name, sequence=seq, structure=struct))
            i += 3
        else:
            i += 1

    after_mismatch, n_mismatch = filter_length_mismatch(raw)
    after_noncanonical, n_noncanonical = filter_noncanonical(after_mismatch)
    after_length, n_length = filter_length(after_noncanonical, min_len, max_len)

    stats = {
        "raw": len(raw),
        "after_mismatch": len(after_mismatch),
        "after_noncanonical": len(after_noncanonical),
        "after_length": len(after_length),
        "n_mismatch": n_mismatch,
        "n_noncanonical": n_noncanonical,
        "n_length": n_length,
    }
    return after_length, stats
