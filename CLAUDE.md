# CLAUDE.md

## Project Overview

This project uses machine learning to re-estimate the **11 pseudoknot-specific
free-energy parameters of the Dirks-Pierce (DP) energy model**, aiming to improve on
the published DP09 values (Andronescu, Pop & Condon 2010).

The deliverable is a small JSON file of 11 numbers — a better-performing replacement
for the collaborator's current `rna_pk_DirksPierce09_HotKnotsV2.json`. This is NOT a
structure-prediction project and uses NO neural network for structure prediction. The
energy model stays explicit and physical; ML only estimates 11 coefficients.

**Timeline:** ~8 weeks, ending 23 July 2026. Supervisor: Dr. Hosna Jabbari, U. Alberta.
Weekly Friday one-on-one + group meeting.

## The 11 Parameters (the entire estimation target)

The collaborator's file `rna_pk_DirksPierce09_HotKnotsV2.json` holds the current
values. Each maps to a DP pseudoknot feature (Table 6, Andronescu et al. 2010).
Additive parameters are stored as kcal/mol x 100 (integer); the two multiplicative
parameters are stored as plain decimals.

| JSON key          | Feature                                          | DP09 value |
|-------------------|--------------------------------------------------|------------|
| pk_in_ext         | initiate an external pseudoknot                  | -1.38      |
| pk_in_mloop       | initiate a pseudoknot inside a multiloop         | 10.07      |
| pk_in_pk          | initiate a pseudoknot inside a pseudoknot        | 15.00      |
| band              | initiate a band (pseudoknotted stem)             | 2.46       |
| unpaired_in_pk    | unpaired base in a pseudoknot                    | 0.06       |
| cr_in_pk          | nested substructure (closed region) in a pk      | 0.96       |
| pk_stack_x        | MULTIPLICATIVE: stacked pair spanning a band     | 0.89       |
| pk_internal_x     | MULTIPLICATIVE: internal loop spanning a band    | 0.74       |
| pk_mloop_init     | initiate a multiloop that spans a band           | 3.41       |
| pk_mloop_bp       | branch in a multiloop that spans a band          | 0.56       |
| pk_mloop_unpaired | unpaired base in a multiloop that spans a band   | 0.12       |

`pk_stack_x` and `pk_internal_x` are multiplicative: the energy of a band-spanning
stacked pair / internal loop is its MT energy TIMES this factor. The other 9 are
additive penalties. This distinction matters: the DP energy is linear in the 9
additive parameters but NOT linear in the 2 multiplicative ones.

## The Energy Model

DP free energy of a structure:  dG  =  (MT contributions)  +  (pseudoknot contributions)

- MT contributions: stacking, loops, mismatches, dangles for the nested parts.
  Computed from the MT09 parameter file `rna_DirksPierce09.par` (ViennaRNA format).
  These are FIXED in this project; they are not re-estimated.
- Pseudoknot contributions: a sum over pseudoknot features, each = (feature count) x
  (parameter). This is the part governed by the 11 parameters above.

**Reference worked example (use as the primary unit test):** Figure 4 of the
supplementary material computes the DP09 energy of an H-type pseudoknot by hand and
gets **-8.02 kcal/mol**. The implemented energy function MUST reproduce this exact
value on that structure before it is trusted. Figure 6 gives a complex (3-band)
worked example as a second check.

## Training Signal: Structural Data

There is almost no thermodynamic (measured-energy) data for pseudoknots, so this
project trains on STRUCTURAL data: sequences with known reference structures, no
energy labels. The estimation principle (constraint generation, after Andronescu et
al.): tune the parameters so that, for each sequence, the reference structure has
LOWER free energy than competing structures for the same sequence.

Because the DP energy is linear in the 9 additive parameters, "reference beats
competitor" is a linear inequality on those parameters; estimation becomes a
max-margin / structured-SVM problem (solvable with scipy). The 2 multiplicative
parameters make the full problem non-linear; options are to (a) hold the 2
multiplicative parameters fixed at DP09 values and estimate only the 9 additive ones,
or (b) alternate / grid-search over the 2 multiplicative ones. Decide with supervisor.

## OPEN QUESTION — resolve before building the estimation method

**How are competing structures generated?** Two options:
  (a) run a pseudoknot folding algorithm (e.g. HotKnots) with current parameters and
      take its predicted structures as competitors — faithful but heavy, puts a
      folder in the training loop;
  (b) generate competitors by perturbing the reference structure — lighter, but
      lower-quality competitors.
This determines the estimation pipeline. ASK THE SUPERVISOR. Do not build the
estimation method (Phase 4 below) until this is answered. Phases 1-3 are unaffected
and should proceed.

## Data

- **Collaborator's parameter file** `rna_pk_DirksPierce09_HotKnotsV2.json` — the 11
  current values; baseline and initialization.
- **MT parameter file** `rna_DirksPierce09.par` — ViennaRNA-format MT09 parameters
  for the nested-structure energy. Fixed input, not estimated.
- **Structural data** — pseudoknotted RNA sequence/structure pairs. Source: the
  PseudoRNA dataset at
  https://www.cs.ubc.ca/labs/algorithms/Publications/PaperMaterials/PseudoRNA/
  (the S-Train / S-Test sets from Andronescu et al.). Download and inspect first.
- **ShPK evaluation set** — 78 short pseudoknotted structures named in Tables 1-2 of
  the supplementary material; DP09 scores F-measure 0.817 on it (DP03: 0.616).

## Evaluation

Primary: structure-prediction F-measure (sensitivity, PPV, harmonic mean) on the
ShPK set, computed by folding with the estimated parameters and comparing to
reference structures. Target: match or beat DP09's 0.817. Report DP03 (0.616) and
DP09 (0.817) as the two reference points.

Secondary: inspect the 11 estimated values for physical plausibility against DP09.

## Build Order (Phases)

1. **Data + parameter I/O.** Download the PseudoRNA structural data. Port the
   pseudoknot-aware structure parser/cleaner from the previous project. Write loaders
   for the JSON pseudoknot parameters and the .par MT file.
2. **Feature counting.** Given a parsed structure, identify bands, pseudoloops,
   span-band loops (per Rastegari & Condon 2007) and count the 11 DP pseudoknot
   features. HIGHEST-RISK COMPONENT — a wrong count corrupts every energy.
3. **Energy function.** Implement dG = MT + pseudoknot contributions. MUST reproduce
   the Figure 4 worked example (-8.02 kcal/mol) and the Figure 6 example before being
   trusted.
4. **STOP — resolve the competing-structures open question with the supervisor.**
5. **Estimation method.** Max-margin estimation of the 9 additive parameters (and a
   strategy for the 2 multiplicative ones), trained on structural data.
6. **Evaluation.** Fold ShPK with estimated parameters; compare F-measure to DP09.

## Coding Standards

- Python 3.10+. numpy / scipy / pandas. NO deep-learning framework — the estimation
  is max-margin / constrained optimization, not a neural network. Do not add torch.
- Black formatting, ruff linting, type hints on public functions.
- Per-function docstrings on all non-trivial functions.
- Do not reuse variable names across separate branches of a function.
- Energy parameters are kcal/mol x 100 integers for additive features; keep the unit
  convention explicit in code and documented at every boundary.
- Unit-test feature counting and the energy function against the hand-computed
  worked examples. These are not optional.
- Never commit datasets; gitignore data/.

## Repository Structure

```
dp-pk-params/
├── CLAUDE.md
├── README.md
├── pyproject.toml
├── .gitignore
├── data/
│   └── README.md
├── params/                  # the JSON and .par parameter files live here (small,
│                            #   may be committed unlike bulk data)
├── src/
│   ├── data/                # structure parsing & cleaning (ported, proven)
│   ├── params/              # load/save JSON pk params and .par MT params
│   ├── features/            # band/pseudoloop detection, 11-feature counting
│   ├── energy/              # dG = MT + pseudoknot energy function
│   ├── estimation/          # max-margin parameter estimation (Phase 5)
│   └── evaluation/          # F-measure on ShPK vs DP09 baseline
├── scripts/
│   ├── download_data.py
│   ├── count_features.py
│   ├── estimate.py          # built only after the open question is resolved
│   └── evaluate.py
├── configs/
├── tests/
│   └── fixtures/            # the Figure 4 / Figure 6 worked examples
├── notebooks/
└── runs/
```

## What Claude Code Should Do

- Build Phases 1-3 now; they are unaffected by the open question.
- Treat feature counting (Phase 2) as the highest-risk component. Unit-test it hard.
- Before trusting the energy function, reproduce the Figure 4 example (-8.02 kcal/mol)
  exactly. If it does not match, stop and debug.
- Do NOT build the estimation method until the competing-structures question is
  answered. If asked to, point back to this file.
- Reuse the parser/cleaner from the previous project; do not rewrite them.
- Keep the kcal/mol x 100 integer convention explicit at every boundary.

## What Claude Code Should Avoid

- Any neural-network or structure-prediction model.
- Adding a deep-learning framework.
- Re-estimating the MT parameters or the Cao-Chen model — out of scope. Only the 11
  DP pseudoknot parameters.
- Guessing the answer to the open question instead of flagging it.

## References

- Andronescu, M. S., Pop, C., & Condon, A. E. (2010). Improved free energy parameters
  for RNA pseudoknotted secondary structure prediction. RNA, 16(1), 26-42.
  PRIMARY REFERENCE — Table 6 (the 11 parameters), supplementary Figures 4 & 6
  (worked energy examples), Tables 1-2 (the ShPK evaluation set).
- Dirks, R. M., & Pierce, N. A. (2003). A partition function algorithm for nucleic
  acid secondary structure including pseudoknots.
- Rastegari, B., & Condon, A. (2007). Parsing nucleic acid pseudoknotted secondary
  structure. — band / pseudoloop / span-band definitions used for feature counting.
- Ren, J., Rastegari, B., Condon, A., & Hoos, H. H. (2005). HotKnots. — the folding
  algorithm referenced for competing-structure generation.

## Contacts

- Supervisor: Dr. Hosna Jabbari (jabbari@ualberta.ca).
- Collaborator: Yahya — maintains the HotKnots-based pseudoknot program and the
  parameter JSON.
