"""Load and save the 11 DP pseudoknot parameters.

The collaborator's JSON stores additive parameters as kcal/mol x 100 (integers) and
the two multiplicative parameters (pk_stack_x, pk_internal_x) as plain decimals.
This module provides a consistent in-memory representation and keeps the unit
convention explicit. See the parameter table in CLAUDE.md.
"""
