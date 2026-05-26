"""Tests for src/params/pk_params.py and src/params/mt_params.py."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from src.params.pk_params import (
    ALL_KEYS,
    ADDITIVE_KEYS,
    MULTIPLICATIVE_KEYS,
    PKParams,
    load_pk_params,
    save_pk_params,
)
from src.params.mt_params import (
    PAIR_INDEX,
    MTParams,
    load_mt_params,
    stack_energy,
    _NST,
)

# Locate parameter files relative to the repo root.
REPO_ROOT = Path(__file__).parent.parent
PK_JSON = REPO_ROOT / "params" / "rna_pk_DirksPierce09_HotKnotsV2.json"
MT_PAR = REPO_ROOT / "params" / "rna_DirksPierce09.par"


# ---------------------------------------------------------------------------
# PKParams tests
# ---------------------------------------------------------------------------


class TestLoadPKParams:
    def test_all_keys_present(self):
        params = load_pk_params(PK_JSON)
        for key in ALL_KEYS:
            assert hasattr(params, key), f"Missing attribute: {key}"

    def test_additive_values_match_json(self):
        """Raw stored values must match the collaborator's JSON exactly."""
        params = load_pk_params(PK_JSON)
        assert params.pk_in_ext == -138
        assert params.pk_in_mloop == 1007
        assert params.pk_in_pk == 1500
        assert params.band == 246
        assert params.unpaired_in_pk == 6
        assert params.cr_in_pk == 96
        assert params.pk_mloop_init == 341
        assert params.pk_mloop_bp == 56
        assert params.pk_mloop_unpaired == 12

    def test_multiplicative_values_match_json(self):
        params = load_pk_params(PK_JSON)
        assert params.pk_stack_x == pytest.approx(0.89)
        assert params.pk_internal_x == pytest.approx(0.74)


class TestToKcal:
    def test_additive_conversion(self):
        """Additive parameters must be divided by 100 to get kcal/mol."""
        params = load_pk_params(PK_JSON)
        kcal = params.to_kcal()
        assert kcal["pk_in_ext"] == pytest.approx(-1.38)
        assert kcal["band"] == pytest.approx(2.46)
        assert kcal["unpaired_in_pk"] == pytest.approx(0.06)
        assert kcal["pk_mloop_unpaired"] == pytest.approx(0.12)

    def test_multiplicative_passthrough(self):
        """Multiplicative parameters must NOT be scaled."""
        params = load_pk_params(PK_JSON)
        kcal = params.to_kcal()
        assert kcal["pk_stack_x"] == pytest.approx(0.89)
        assert kcal["pk_internal_x"] == pytest.approx(0.74)

    def test_all_keys_present_in_kcal(self):
        params = load_pk_params(PK_JSON)
        kcal = params.to_kcal()
        for key in ALL_KEYS:
            assert key in kcal, f"Missing key in to_kcal() output: {key}"


class TestSavePKParams:
    def test_roundtrip(self):
        """save then load must reproduce the original PKParams exactly."""
        original = load_pk_params(PK_JSON)
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            save_pk_params(original, tmp_path, name="test-roundtrip")
            reloaded = load_pk_params(tmp_path)
            for key in ADDITIVE_KEYS:
                assert getattr(original, key) == getattr(reloaded, key), key
            for key in MULTIPLICATIVE_KEYS:
                assert getattr(original, key) == pytest.approx(getattr(reloaded, key)), key
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_saved_json_structure(self):
        """The JSON must have the pseudoknot_parameters top-level key."""
        params = load_pk_params(PK_JSON)
        with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            save_pk_params(params, tmp_path)
            with tmp_path.open() as fh:
                data = json.load(fh)
            assert "pseudoknot_parameters" in data
            pk = data["pseudoknot_parameters"]
            for key in ALL_KEYS:
                assert key in pk, f"Missing key in saved JSON: {key}"
        finally:
            tmp_path.unlink(missing_ok=True)

    def test_missing_key_raises(self):
        """A JSON missing a required key must raise KeyError."""
        bad_json = {"pseudoknot_parameters": {"pk_in_ext": -138}}
        with tempfile.NamedTemporaryFile(
            suffix=".json", mode="w", delete=False
        ) as tmp:
            json.dump(bad_json, tmp)
            tmp_path = Path(tmp.name)
        try:
            with pytest.raises(KeyError):
                load_pk_params(tmp_path)
        finally:
            tmp_path.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# MTParams tests
# ---------------------------------------------------------------------------


class TestLoadMTParams:
    def test_stack_shape(self):
        """Stack table must be 7x7."""
        mt = load_mt_params(MT_PAR)
        assert mt.stack.shape == (7, 7)

    def test_stack_dtype(self):
        mt = load_mt_params(MT_PAR)
        assert mt.stack.dtype == np.int32

    def test_stack_cg_cg(self):
        """CG-CG stacking energy must be -152 (kcal/mol x 100)."""
        mt = load_mt_params(MT_PAR)
        assert stack_energy(mt, "CG", "CG") == -152

    def test_stack_gc_gc(self):
        """GC-GC stacking energy must be -216 (kcal/mol x 100)."""
        mt = load_mt_params(MT_PAR)
        assert stack_energy(mt, "GC", "GC") == -216

    def test_stack_au_au(self):
        """AU-AU stacking energy must be -74 (kcal/mol x 100)."""
        mt = load_mt_params(MT_PAR)
        assert stack_energy(mt, "AU", "AU") == -74

    def test_nst_entries_are_sentinel(self):
        """NN row/column entries must be the NST sentinel."""
        mt = load_mt_params(MT_PAR)
        # NN-CG is NST in the file
        assert stack_energy(mt, "NN", "CG") == _NST

    def test_other_sections_present(self):
        """Other sections must be stored as raw text blocks."""
        mt = load_mt_params(MT_PAR)
        assert "mismatch_hairpin" in mt.raw_sections
        assert "mismatch_interior" in mt.raw_sections
