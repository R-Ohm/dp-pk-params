# dp-pk-params

Machine-learning re-estimation of the 11 pseudoknot-specific free-energy parameters
of the Dirks-Pierce (DP) energy model, aiming to improve on the published DP09 set
(Andronescu, Pop & Condon 2010).

The deliverable is a small JSON file of 11 numbers: a better-performing replacement
for the collaborator's current pseudoknot parameter file. This is a
parameter-estimation project, not a structure-prediction project.

See CLAUDE.md for the full description, the parameter table, the energy model, the
build order, and the open question that must be resolved before the estimation
method is built.

## Setup
```bash
python -m venv .venv
# Windows: .\.venv\Scripts\Activate.ps1
source .venv/bin/activate
pip install -e ".[dev]"
```

## Status
Early. Phases 1-3 (data, feature counting, energy function) proceed now; the
estimation method is blocked on a supervisor question (see CLAUDE.md).
