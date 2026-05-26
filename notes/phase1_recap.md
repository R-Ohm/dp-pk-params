# Phase 1 Recap — Parameter I/O

**Goal:** Load the two parameter files into clean in-memory representations so all downstream code has a single, trusted source of truth for the parameters.

---

## What was built

**`src/params/pk_params.py`** — loader for the 11 DP pseudoknot parameters

- Defines a `PKParams` dataclass holding all 11 parameters in their native stored units: the 9 additive parameters as kcal/mol×100 integers, and the 2 multiplicative parameters (`pk_stack_x`, `pk_internal_x`) as plain decimals.
- `load_pk_params()` reads the collaborator's JSON and validates that all 11 keys are present.
- `PKParams.to_kcal()` converts the 9 additive parameters to kcal/mol floats while leaving the multiplicative ones unchanged.
- `save_pk_params()` writes a parameter set back to JSON in the same format — needed later when we produce the re-estimated parameter file.

**`src/params/mt_params.py`** — parser for the ViennaRNA-format MT parameter file

- Parses `rna_DirksPierce09.par` and extracts the 7×7 stacking energy table (pair types: CG, GC, GU, UG, AU, UA, NN) into a numpy array.
- All other sections (mismatch, dangle, loop penalties, etc.) are stored as raw text blocks so they can be parsed on demand in later phases without re-reading the file.
- These MT parameters are fixed inputs — they are not estimated.

**`tests/test_pk_params.py`** — 16 unit tests, all passing

- Checks all 11 keys load correctly, that the kcal/mol conversion is exact (e.g. −138 → −1.38, 246 → 2.46), that multiplicative parameters are not scaled, and that a save/load round-trip is lossless.
- Checks the stack table is the right shape and dtype, and spot-checks several entries against known values from the file (CG-CG = −152, GC-GC = −216, AU-AU = −74).

---

## Status

Phase 1 complete. Ready to move to Phase 2: structure parsing (dot-bracket → base-pair set, multi-bracket support, cleaning helpers).
