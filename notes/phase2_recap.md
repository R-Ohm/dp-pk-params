# Phase 2 Recap — Structure Parsing

**Goal:** Turn a dot-bracket structure string into a set of base pairs (pseudoknots included), provide filtering helpers to clean raw dataset records, and load the actual training/test data files into memory.

---

## What was built

**`src/data/parse.py`** — dot-bracket parser and shared record type

- `RNARecord` dataclass holding `(name, sequence, structure)` — the shared data type used across the pipeline.
- `parse_structure(s: str) -> set[tuple[int, int]]` — parses a dot-bracket string into 0-indexed base pairs `(i, j)` with `i < j`. Each bracket type (`()`, `[]`, `{}`, `<>`, and letter pairs `Aa`–`Zz`) maintains its own independent LIFO stack, so pseudoknots encoded with mixed bracket types are handled correctly. The parser is intentionally pure — no minimum loop length filtering, no sequence awareness.

**`src/data/clean.py`** — record-level filtering helpers

- `filter_length_mismatch(records)` — removes records where `len(sequence) != len(structure)`.
- `filter_noncanonical(records, alphabet)` — removes records containing bases outside the allowed alphabet (default: `{A, C, G, U}`).
- `filter_length(records, min_len, max_len)` — removes records outside a length window; both bounds are optional and inclusive.
- All three functions return `(kept_list, n_removed)` so they can be chained and logged.

**`src/data/hotknots_reader.py`** — file reader for the HotKnots dataset format (added as a gap fix)

- `load_hotknots_file(path, min_len, max_len)` — parses a three-line-per-record `.txt` file (header `>`, sequence, structure) into a cleaned list of `RNARecord`, then applies all three cleaning filters in sequence and returns the records plus a stats dict.
- `normalize_sequence(seq)` — uppercases the sequence and substitutes `T→U` (DNA notation), applied before filtering so these records are kept rather than dropped.
- Handles two header formats: key-value (`PSEUDOBASE_ID=...`, `SSTRAND_ID=...`) and bare-name (`> AKV-MuLV`).
- Verified load stats:
  - **S-Train** (`TRA_HotKnots_max200_S80P.txt`): 1807 records loaded, 0 dropped, all `len(sequence) == len(structure)`.
  - **S-Test / ShPK** (`TES_pkshort_max200_RenJab-S20.txt`): 78 records loaded, 0 dropped, all lengths consistent.

**`tests/test_parse.py`** — 25 unit tests, all passing

- Hairpin `((((....))))` → 4 nested pairs.
- H-type pseudoknot `(((...[[[)))...]]]` → 6 pairs across two independent stacks.
- All-dot string → empty set.
- Triple-bracket pseudoknot `([{...}])` using `()`, `[]`, `{}` → 3 pairs, one per stack.
- Letter bracket pairs (`A..a`, `ABba`).
- Error cases: unmatched open, unmatched close, mismatched bracket types, invalid character.
- Full coverage of all three `clean.py` filters including edge cases (empty input, custom alphabet, inclusive bounds).

---

## Status

Phase 2 complete. Ready to move to Phase 3: feature counting (band/pseudoloop detection and the 11 DP pseudoknot feature counts). This is the highest-risk component — requires the Rastegari & Condon (2007) definitions to be confirmed before writing any code.
