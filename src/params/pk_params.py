"""Load and save the 11 DP pseudoknot parameters.

The collaborator's JSON stores additive parameters as kcal/mol x 100 (integers) and
the two multiplicative parameters (pk_stack_x, pk_internal_x) as plain decimals.
This module provides a consistent in-memory representation and keeps the unit
convention explicit. See the parameter table in CLAUDE.md.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

# Ordered list of the 9 additive parameter keys (stored as kcal/mol x 100 integers).
ADDITIVE_KEYS: tuple[str, ...] = (
    "pk_in_ext",
    "pk_in_mloop",
    "pk_in_pk",
    "band",
    "unpaired_in_pk",
    "cr_in_pk",
    "pk_mloop_init",
    "pk_mloop_bp",
    "pk_mloop_unpaired",
)

# The two multiplicative parameters are plain decimals (no unit scaling).
MULTIPLICATIVE_KEYS: tuple[str, ...] = (
    "pk_stack_x",
    "pk_internal_x",
)

ALL_KEYS: tuple[str, ...] = ADDITIVE_KEYS + MULTIPLICATIVE_KEYS


@dataclass
class PKParams:
    """The 11 DP pseudoknot free-energy parameters.

    Additive parameters are stored in kcal/mol x 100 (integer hundredths); use
    ``to_kcal()`` to get plain kcal/mol floats. Multiplicative parameters are
    plain decimals and are NOT scaled.
    """

    # 9 additive parameters (kcal/mol x 100, integers)
    pk_in_ext: int
    pk_in_mloop: int
    pk_in_pk: int
    band: int
    unpaired_in_pk: int
    cr_in_pk: int
    pk_mloop_init: int
    pk_mloop_bp: int
    pk_mloop_unpaired: int

    # 2 multiplicative parameters (plain decimals)
    pk_stack_x: float
    pk_internal_x: float

    def to_kcal(self) -> dict[str, float]:
        """Return all parameters in kcal/mol.

        Additive parameters are divided by 100; multiplicative parameters are
        returned unchanged.
        """
        result: dict[str, float] = {}
        for key in ADDITIVE_KEYS:
            result[key] = getattr(self, key) / 100.0
        for key in MULTIPLICATIVE_KEYS:
            result[key] = getattr(self, key)
        return result


def load_pk_params(path: str | Path) -> PKParams:
    """Load the 11 DP pseudoknot parameters from a JSON file.

    The JSON must have a ``pseudoknot_parameters`` top-level key containing all
    11 parameter fields. Additive parameters must be integers (kcal/mol x 100);
    multiplicative parameters must be floats.

    Parameters
    ----------
    path:
        Path to the parameter JSON file (e.g. ``params/rna_pk_DirksPierce09_HotKnotsV2.json``).

    Returns
    -------
    PKParams
        All 11 parameters in their native stored units.

    Raises
    ------
    KeyError
        If any of the 11 required keys is missing from the JSON.
    ValueError
        If the file cannot be parsed as valid JSON.
    """
    path = Path(path)
    with path.open("r") as fh:
        raw = json.load(fh)

    params_raw: dict = raw["pseudoknot_parameters"]

    missing = [k for k in ALL_KEYS if k not in params_raw]
    if missing:
        raise KeyError(f"Missing required parameter keys: {missing}")

    return PKParams(
        pk_in_ext=int(params_raw["pk_in_ext"]),
        pk_in_mloop=int(params_raw["pk_in_mloop"]),
        pk_in_pk=int(params_raw["pk_in_pk"]),
        band=int(params_raw["band"]),
        unpaired_in_pk=int(params_raw["unpaired_in_pk"]),
        cr_in_pk=int(params_raw["cr_in_pk"]),
        pk_mloop_init=int(params_raw["pk_mloop_init"]),
        pk_mloop_bp=int(params_raw["pk_mloop_bp"]),
        pk_mloop_unpaired=int(params_raw["pk_mloop_unpaired"]),
        pk_stack_x=float(params_raw["pk_stack_x"]),
        pk_internal_x=float(params_raw["pk_internal_x"]),
    )


def save_pk_params(params: PKParams, path: str | Path, name: str = "") -> None:
    """Write a PKParams object to JSON in the collaborator's file format.

    The output mirrors the structure of ``rna_pk_DirksPierce09_HotKnotsV2.json``:
    a top-level ``pseudoknot_parameters`` dict with additive values as integers
    and multiplicative values as floats.

    Parameters
    ----------
    params:
        The parameter set to save.
    path:
        Destination file path.
    name:
        Optional descriptive name stored inside the JSON.
    """
    path = Path(path)
    raw = asdict(params)
    payload: dict = {"pseudoknot_parameters": {"name": name, **raw}}
    with path.open("w") as fh:
        json.dump(payload, fh, indent=2)
