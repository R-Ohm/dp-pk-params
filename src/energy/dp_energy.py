"""DP free-energy function: dG = MT contributions + pseudoknot contributions.

The pseudoknot part sums feature_count * parameter over the 11 features; the two
multiplicative parameters scale band-spanning stacked-pair and internal-loop
energies. MUST reproduce the supplementary Figure 4 worked example (-8.02 kcal/mol)
before being trusted.
"""
