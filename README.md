This repository contains Integration_Patterns, an R function designed to quantify and compare morphological integration patterns across groups using two-block Partial Least Squares (PLS) and within-group covariance structures.

The implementation is accompanied by an example dataset of hominoid cranial and C5 vertebra landmark coordinates, available on Zenodo:

Dataset DOI: https://doi.org/10.5281/zenodo.19203426

The dataset is provided as a reproducible example for testing and illustrating the Integration_Patterns function.

Method overview

The method builds upon the PLS framework implemented in the R package Morpho and introduces a geometric approach to compare integration patterns across groups (e.g., species).

The workflow is as follows:

1. A pooled two-block PLS is computed using Morpho::pls2B on two matched datasets (e.g., cranial and cervical vertebrae coordinates).
2. The first PLS axis (PLS1) is extracted for each block, defining a shared PLS1–PLS1 space.
3. For each group:
   - A covariance matrix is computed from the PLS1 scores.
   - The first eigenvector of this matrix is taken as the dominant axis of covariation (integration axis).
4. Differences between groups are quantified as angular distances between these axes (0–90°).
5. Statistical significance is assessed using permutation tests that randomly reassign individuals between groups.

By default, the function performs pooled within-group PLS (pooled_within = TRUE), removing group means prior to PLS estimation and focusing on within-group covariance structure.