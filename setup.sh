#!/usr/bin/env bash
# setup.sh — initialize the dp-pk-params repository structure. Run once.

set -e
echo "Setting up dp-pk-params project structure..."

mkdir -p src/{data,params,features,energy,estimation,evaluation}
mkdir -p scripts configs tests/fixtures notebooks runs data params docs

touch src/__init__.py
for d in data params features energy estimation evaluation; do
  touch "src/$d/__init__.py"
done
touch tests/__init__.py

cat > .gitignore <<'EOF'
__pycache__/
*.py[cod]
*.egg-info/
.pytest_cache/
.ruff_cache/
.mypy_cache/
venv/
.venv/
env/

# Bulk data — not committed
data/raw/
data/processed/
*.parquet

# Run outputs
runs/
*.log

.ipynb_checkpoints/
.vscode/
.idea/
.DS_Store

# NOTE: params/ is intentionally NOT ignored — the parameter files are small
# and are inputs/outputs worth version-controlling.
EOF

cat > data/README.md <<'EOF'
# Data

Bulk data is gitignored. Download structural data here.

## PseudoRNA structural dataset
Pseudoknotted RNA sequence/structure pairs (the S-Train / S-Test sets from
Andronescu et al. 2010). Source:
https://www.cs.ubc.ca/labs/algorithms/Publications/PaperMaterials/PseudoRNA/
Download, inspect the format, then parse with src/data/.

## ShPK evaluation set
78 short pseudoknotted structures, named in Tables 1-2 of the Andronescu
supplementary material. Used for the F-measure evaluation against DP09.
EOF

cat > params/README.md <<'EOF'
# Parameter files

Small enough to version-control; kept here rather than in data/.

- rna_pk_DirksPierce09_HotKnotsV2.json  -- the 11 DP pseudoknot parameters
  (collaborator's file; baseline DP09 values; estimation initialization).
- rna_DirksPierce09.par                 -- ViennaRNA-format MT09 parameters for the
  nested-structure energy. Fixed input; not re-estimated.

Place both files here. Estimated parameter sets produced by this project are also
written here, with descriptive names (e.g. dp_pk_estimated_v1.json).
EOF

cat > pyproject.toml <<'EOF'
[project]
name = "dp-pk-params"
version = "0.1.0"
description = "ML re-estimation of the 11 Dirks-Pierce pseudoknot energy parameters"
requires-python = ">=3.10"
dependencies = [
    "numpy>=1.24",
    "scipy>=1.11",
    "pandas>=2.0",
    "matplotlib>=3.7",
    "pyyaml>=6.0",
    "tqdm>=4.65",
]
# No deep-learning framework: estimation is max-margin / constrained optimization.

[project.optional-dependencies]
dev = ["pytest>=7.4", "ruff>=0.1", "black>=23.0", "mypy>=1.5"]

[tool.ruff]
line-length = 100
target-version = "py310"

[tool.black]
line-length = 100

[tool.pytest.ini_options]
testpaths = ["tests"]
EOF

cat > README.md <<'EOF'
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
EOF

# --- source stubs with intent docstrings ---
cat > src/params/pk_params.py <<'EOF'
"""Load and save the 11 DP pseudoknot parameters.

The collaborator's JSON stores additive parameters as kcal/mol x 100 (integers) and
the two multiplicative parameters (pk_stack_x, pk_internal_x) as plain decimals.
This module provides a consistent in-memory representation and keeps the unit
convention explicit. See the parameter table in CLAUDE.md.
"""
EOF

cat > src/params/mt_params.py <<'EOF'
"""Parse the ViennaRNA-format MT parameter file (rna_DirksPierce09.par).

Reads the stacking, mismatch, loop, dangle and other Mathews-Turner sections.
These parameters are fixed inputs to the energy function; they are not estimated.
"""
EOF

cat > src/features/pk_features.py <<'EOF'
"""Detect bands / pseudoloops / span-band loops and count the 11 DP pseudoknot
features for a given structure.

Follows the loop/band classification of Rastegari & Condon (2007). HIGHEST-RISK
component of the project: a wrong count silently corrupts every energy. Unit-test
against the hand-computed worked examples (supplementary Figures 4 and 6).
"""
EOF

cat > src/energy/dp_energy.py <<'EOF'
"""DP free-energy function: dG = MT contributions + pseudoknot contributions.

The pseudoknot part sums feature_count * parameter over the 11 features; the two
multiplicative parameters scale band-spanning stacked-pair and internal-loop
energies. MUST reproduce the supplementary Figure 4 worked example (-8.02 kcal/mol)
before being trusted.
"""
EOF

cat > src/estimation/README.md <<'EOF'
Phase 5: max-margin estimation of the 9 additive parameters (plus a strategy for the
2 multiplicative ones). DO NOT IMPLEMENT until the competing-structures open
question in CLAUDE.md is resolved with the supervisor.
EOF

cat > src/evaluation/fmeasure.py <<'EOF'
"""F-measure evaluation on the ShPK set.

Sensitivity, PPV and their harmonic mean over predicted vs reference base pairs.
Baselines for comparison: DP03 = 0.616, DP09 = 0.817 (Andronescu et al. 2010).
"""
EOF

cat > tests/test_pk_features.py <<'EOF'
"""Tests for the 11-feature counting. Verify against the hand-computed worked
examples in the Andronescu supplementary material (Figures 4 and 6)."""

def test_placeholder():
    assert True
EOF

cat > tests/test_dp_energy.py <<'EOF'
"""Energy-function tests. The key test: the Figure 4 H-type pseudoknot must compute
to -8.02 kcal/mol with DP09 parameters. Figure 6 is a second worked example."""

def test_placeholder():
    assert True
EOF

cat > scripts/download_data.py <<'EOF'
"""Download the PseudoRNA structural dataset (S-Train / S-Test).

Source: https://www.cs.ubc.ca/labs/algorithms/Publications/PaperMaterials/PseudoRNA/
Inspect the file format before parsing.
"""

def main() -> None:
    raise NotImplementedError("Implement in Phase 1.")

if __name__ == "__main__":
    main()
EOF

cat > scripts/count_features.py <<'EOF'
"""Produce 11-feature count vectors for a set of parsed structures."""

def main() -> None:
    raise NotImplementedError("Implement after src/features/pk_features.py (Phase 2).")

if __name__ == "__main__":
    main()
EOF

cat > scripts/estimate.py <<'EOF'
"""Run max-margin estimation of the DP pseudoknot parameters.

DO NOT IMPLEMENT until the competing-structures open question in CLAUDE.md is
resolved with the supervisor (Phase 4).
"""

def main() -> None:
    raise NotImplementedError("Blocked on the CLAUDE.md open question.")

if __name__ == "__main__":
    main()
EOF

cat > scripts/evaluate.py <<'EOF'
"""Fold the ShPK set with a parameter set and report F-measure vs DP09."""

def main() -> None:
    raise NotImplementedError("Implement in Phase 6.")

if __name__ == "__main__":
    main()
EOF

echo "Done. Repository structure created."
echo ""
echo "Next steps:"
echo "  1. git init && git add . && git commit -m 'Initial scaffold'"
echo "  2. python -m venv .venv && activate && pip install -e '.[dev]'"
echo "  3. Place rna_pk_DirksPierce09_HotKnotsV2.json and rna_DirksPierce09.par"
echo "     into params/."
echo "  4. Read CLAUDE.md — note the parameter table and the one open question."
echo "  5. Phase 1: download the PseudoRNA data; port the parser from the old repo."
