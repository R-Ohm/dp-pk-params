"""Tests for src/data/parse.py and src/data/clean.py."""

import pytest

from src.data.parse import RNARecord, parse_structure
from src.data.clean import (
    CANONICAL_BASES,
    filter_length,
    filter_length_mismatch,
    filter_noncanonical,
)


# ---------------------------------------------------------------------------
# parse_structure tests
# ---------------------------------------------------------------------------


class TestParseStructure:
    def test_hairpin(self):
        """((((....)))): 4 nested Watson-Crick pairs."""
        pairs = parse_structure("((((....))))")
        assert pairs == {(0, 11), (1, 10), (2, 9), (3, 8)}

    def test_h_type_pseudoknot(self):
        """(((...[[[)))...]]] : 6 pairs across two independent stacks."""
        # ( stack: 0,1,2  [ stack: 6,7,8
        # ) at 9->2, 10->1, 11->0   ] at 15->8, 16->7, 17->6
        pairs = parse_structure("(((...[[[)))...]]]")
        assert pairs == {(0, 11), (1, 10), (2, 9), (6, 17), (7, 16), (8, 15)}

    def test_empty_string(self):
        """All-dot string yields an empty pair set."""
        assert parse_structure("......") == set()

    def test_triple_bracket_pseudoknot(self):
        """Three-level pseudoknot using (), [], {}: 3 pairs, one from each stack."""
        # ( at 0, [ at 2, { at 4, } at 6, ] at 8, ) at 10
        s = "([{...}])"
        pairs = parse_structure(s)
        assert (0, 8) in pairs   # ( ) pair
        assert (1, 7) in pairs   # [ ] pair
        assert (2, 6) in pairs   # { } pair
        assert len(pairs) == 3

    def test_single_pair(self):
        assert parse_structure("(.)")  == {(0, 2)}

    def test_letter_bracket_uppercase_opens(self):
        """Uppercase letter opens, lowercase closes: A..a -> pair (0, 3)."""
        pairs = parse_structure("A..a")
        assert pairs == {(0, 3)}

    def test_multiple_letter_types_independent(self):
        """Different letter types maintain separate stacks."""
        # A at 0, B at 1, b at 2 (closes B), a at 3 (closes A)
        pairs = parse_structure("ABba")
        assert pairs == {(0, 3), (1, 2)}

    def test_unmatched_open_raises(self):
        with pytest.raises(ValueError, match="Unmatched opening"):
            parse_structure("((.")

    def test_unmatched_close_raises(self):
        with pytest.raises(ValueError, match="Unmatched closing"):
            parse_structure(")..")

    def test_mismatched_bracket_types_raise(self):
        """A closing ] should not match an opening (."""
        with pytest.raises(ValueError):
            parse_structure("(...]")

    def test_invalid_character_raises(self):
        # X is a valid letter-bracket; use ! which is outside the entire alphabet.
        with pytest.raises(ValueError, match="Unrecognised character"):
            parse_structure("(.!.)")

    def test_angle_brackets(self):
        """<> bracket type is supported."""
        pairs = parse_structure("<.>")
        assert pairs == {(0, 2)}

    def test_nested_same_bracket(self):
        """Multiple pairs of the same bracket type nest correctly (LIFO)."""
        pairs = parse_structure("(())")
        assert pairs == {(0, 3), (1, 2)}


# ---------------------------------------------------------------------------
# clean.py tests
# ---------------------------------------------------------------------------


def _records(*specs: tuple[str, str, str]) -> list[RNARecord]:
    """Helper: build RNARecord list from (name, seq, struct) tuples."""
    return [RNARecord(name=n, sequence=s, structure=t) for n, s, t in specs]


class TestFilterLengthMismatch:
    def test_removes_mismatched(self):
        recs = _records(
            ("ok", "ACGU", "...."),
            ("bad", "ACG", "...."),   # 3 vs 4
        )
        kept, n = filter_length_mismatch(recs)
        assert n == 1
        assert len(kept) == 1
        assert kept[0].name == "ok"

    def test_all_valid(self):
        recs = _records(("a", "ACGU", "...."), ("b", "GG", ".."))
        kept, n = filter_length_mismatch(recs)
        assert n == 0
        assert len(kept) == 2

    def test_empty_input(self):
        kept, n = filter_length_mismatch([])
        assert kept == []
        assert n == 0


class TestFilterNoncanonical:
    def test_removes_noncanonical(self):
        recs = _records(
            ("ok", "ACGU", "...."),
            ("bad", "ACGN", "...."),  # N is not canonical
            ("bad2", "acgu", "...."), # lowercase not canonical
        )
        kept, n = filter_noncanonical(recs)
        assert n == 2
        assert kept[0].name == "ok"

    def test_custom_alphabet(self):
        recs = _records(("t", "ACGT", "...."))
        kept, n = filter_noncanonical(recs, alphabet=frozenset("ACGT"))
        assert n == 0

    def test_all_canonical(self):
        recs = _records(("a", "ACGU", "...."), ("b", "GGCC", "...."))
        kept, n = filter_noncanonical(recs)
        assert n == 0
        assert len(kept) == 2

    def test_empty_input(self):
        kept, n = filter_noncanonical([])
        assert kept == []
        assert n == 0


class TestFilterLength:
    def test_min_len(self):
        recs = _records(("short", "AC", ".."), ("long", "ACGU", "...."))
        kept, n = filter_length(recs, min_len=3)
        assert n == 1
        assert kept[0].name == "long"

    def test_max_len(self):
        recs = _records(("short", "AC", ".."), ("long", "ACGU", "...."))
        kept, n = filter_length(recs, max_len=3)
        assert n == 1
        assert kept[0].name == "short"

    def test_both_bounds(self):
        recs = _records(
            ("tiny", "A", "."),
            ("ok", "ACG", "..."),
            ("big", "ACGUA", "....."),
        )
        kept, n = filter_length(recs, min_len=2, max_len=4)
        assert n == 2
        assert kept[0].name == "ok"

    def test_no_bounds(self):
        recs = _records(("a", "A", "."), ("b", "ACGU", "...."))
        kept, n = filter_length(recs)
        assert n == 0
        assert len(kept) == 2

    def test_inclusive_bounds(self):
        """min_len and max_len are inclusive."""
        recs = _records(("exact", "ACG", "..."))
        kept, n = filter_length(recs, min_len=3, max_len=3)
        assert n == 0
        assert len(kept) == 1
