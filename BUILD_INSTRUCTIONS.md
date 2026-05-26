# BUILD_INSTRUCTIONS.md

Instructions for Claude Code to build the `dp-pk-params` project from scratch.
This project is self-contained: build everything fresh here. Do not import or copy
code from any previous project.

Read `CLAUDE.md` first — it has the parameter table, the energy model, the unit
conventions, and the one open question. This file is the ordered build plan.

## Ground rules

- Work in small, tested increments. After each step, run the tests and stop for the
  user to review before continuing. Do not build multiple phases in one turn.
- Show code before writing files. Wait for the user to approve.
- No deep-learning framework. The estimation is constrained optimization (scipy).
- Two highest-risk components: feature counting and the energy function. They get
  unit tests against hand-computed examples — this is mandatory, not optional.
- Keep the energy unit convention explicit everywhere: additive parameters are
  stored as kcal/mol x 100 integers; the two multiplicative parameters
  (pk_stack_x, pk_internal_x) are plain decimals.
- Do NOT build Phase 5 (estimation) until the user confirms the open question about
  competing-structure generation has been answered by their supervisor.

## Phase 0 — Scaffold

1. Run `setup.sh` (already provided) to create the directory tree, `.gitignore`,
   `pyproject.toml`, `README.md`, and stub files.
2. `git init`, then an initial commit.
3. Set up the virtual environment and `pip install -e ".[dev]"`. If anything fails
   to install, report which package and stop — do not work around it silently.
4. Confirm `pytest` runs (the stub tests should pass).

Stop. Let the user confirm the environment is clean before Phase 1.

## Phase 1 — Parameter I/O

Goal: load the two parameter files into clean in-memory representations.

1. `src/params/pk_params.py`: load `params/rna_pk_DirksPierce09_HotKnotsV2.json`
   into a dataclass or dict of 11 named parameters. Provide:
   - a loader that returns the 11 values,
   - a function converting the stored integers (kcal/mol x 100) to kcal/mol floats,
     leaving the two multiplicative parameters as-is,
   - a saver that writes a parameter set back to JSON in the same format.
2. `src/params/mt_params.py`: parse the ViennaRNA-format `rna_DirksPierce09.par`.
   At minimum, parse the `stack` section (the 6x6 stacking table). Other sections
   (mismatch, dangle, loop, etc.) can be parsed as needed later; for now, parse
   `stack` and leave clear hooks for the rest.
3. Tests in `tests/test_pk_params.py`: load the JSON, assert all 11 keys are
   present, assert the kcal/mol conversion (e.g. pk_in_ext -138 -> -1.38, band 246
   -> 2.46), assert pk_stack_x stays 0.89.

Stop for review.

## Phase 2 — Structure parsing

Goal: turn a dot-bracket structure string into a set of base pairs, pseudoknots
included.

1. `src/data/parse.py`: a function `parse_structure(s: str) -> set[tuple[int,int]]`
   returning 0-indexed pairs (i, j) with i < j. Handle multiple bracket types as
   independent LIFO stacks: `()`, `[]`, `{}`, `<>`, and letter pairs `Aa`..`Zz`.
   No minimum-loop-length filtering here — the parser is pure.
2. `src/data/clean.py`: filtering helpers — drop non-canonical sequence characters,
   length filters, sequence/structure length-consistency check. Each filter is its
   own function, returning a count of what it removed.
3. Tests in `tests/test_parse.py`:
   - hairpin `((((....))))` -> 4 nested pairs,
   - H-type pseudoknot `(((...[[[)))...]]]` -> 6 pairs across two stacks,
   - empty `......` -> empty set,
   - a triple-level case using `{}` as well.
   Remember bracket matching is LIFO within each stack.

Stop for review.

## Phase 3 — Feature counting (HIGHEST RISK)

Goal: given a parsed structure, count the 11 DP pseudoknot features.

1. `src/features/bands.py`: detect bands, pseudoloops, closed regions and span-band
   loops, following Rastegari & Condon (2007). The user has that paper; ask them for
   specific definitions if anything is ambiguous rather than guessing.
2. `src/features/pk_features.py`: a function that takes a parsed structure and
   returns counts for the 11 features named in CLAUDE.md (pk_in_ext, band,
   unpaired_in_pk, cr_in_pk, pk_mloop_init, pk_mloop_bp, pk_mloop_unpaired, and the
   contexts for the two multiplicative features). Output: a count vector aligned to
   the 11-parameter order.
3. Tests in `tests/test_pk_features.py`: build the H-type pseudoknot from
   supplementary Figure 4 by hand and assert the feature counts match what the
   figure shows (1 external pseudoknot, 2 bands, 8 unpaired-in-pk bases, etc.).
   This test is mandatory and must pass before Phase 4.

Stop for review. Feature counting errors are silent and catastrophic — do not
proceed until the Figure 4 counts are verified.

## Phase 4 — Energy function

Goal: compute the DP free energy of a structure.

1. `src/energy/dp_energy.py`: implement
   `dG = MT_contributions + pseudoknot_contributions`.
   - pseudoknot part: sum over the 9 additive features of count * parameter, plus
     the 2 multiplicative parameters scaling band-spanning stacked-pair and
     internal-loop energies.
   - MT part: use the stacking parameters from `mt_params.py`. For the Figure 4
     check, the MT stacked-pair energies are given in the figure itself
     (-1.52, -2.16, etc.), so the check can be done with those values directly.
2. THE KEY TEST, in `tests/test_dp_energy.py`: reproduce supplementary Figure 4 —
   the H-type pseudoknot must total **-8.02 kcal/mol** with DP09 parameters. If it
   does not match exactly, stop and debug; do not proceed.
3. Add the Figure 6 complex (3-band) pseudoknot as a second check if its numbers
   can be read off the figure.

Stop for review.

## CHECKPOINT — resolve the open question

Before Phase 5, the user must confirm with their supervisor how competing
structures are generated for training:
  (a) a pseudoknot folding algorithm (HotKnots) in the loop, or
  (b) perturbing reference structures.
Do not build Phase 5 until the user states the answer. Ask them; do not assume.

## Phase 5 — Estimation method (BLOCKED until checkpoint passes)

Goal: estimate the 9 additive parameters (with a strategy for the 2 multiplicative
ones) from structural data.

1. `src/estimation/`: implement max-margin estimation. For each training sequence,
   the reference structure should have lower energy than its competing structures;
   since energy is linear in the 9 additive parameters, each such comparison is a
   linear inequality. Solve the resulting max-margin problem with scipy.
2. For the 2 multiplicative parameters: either hold them fixed at DP09 values, or
   grid-search them in an outer loop — follow whatever the user decided with the
   supervisor.
3. `scripts/estimate.py`: wire it end to end; output a new parameter JSON in
   `params/` (e.g. `dp_pk_estimated_v1.json`).

Stop for review.

## Phase 6 — Evaluation

Goal: compare the estimated parameters against DP09.

1. `src/evaluation/fmeasure.py`: sensitivity, PPV, harmonic-mean F-measure over
   predicted vs reference base pairs.
2. `scripts/evaluate.py`: fold the 78-structure ShPK set with a given parameter set
   and report F-measure. Reference points: DP03 = 0.616, DP09 = 0.817.
3. Produce a comparison: estimated parameters vs DP09, overall and per structure.

Stop for review.

## After Phase 6

Write a short results summary and clean up the repo. The deliverable is the
estimated parameter JSON plus the evaluation showing whether it matches or beats
DP09's 0.817 F-measure on ShPK.

## Reminders for Claude Code

- One phase at a time. Stop for review after each. Show code before writing.
- The two worked examples (Figure 4 -> -8.02 kcal/mol; Figure 6) are the ground
  truth for correctness. Tests against them are mandatory.
- Do not build Phase 5 before the user confirms the open question is resolved.
- No deep-learning framework. scipy is the optimization tool.
- If a definition from Rastegari & Condon (2007) is unclear, ask the user — they
  have the paper. Do not guess at band/pseudoloop definitions.
