"""Detect bands / pseudoloops / span-band loops and count the 11 DP pseudoknot
features for a given structure.

Follows the loop/band classification of Rastegari & Condon (2007). HIGHEST-RISK
component of the project: a wrong count silently corrupts every energy. Unit-test
against the hand-computed worked examples (supplementary Figures 4 and 6).
"""
