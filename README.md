# Integration Patterns

This repository contains `Integration_Patterns`, an R implementation designed to quantify and compare morphological integration patterns across species using two-block Partial Least Squares (PLS) and within-species covariance structures.

The method extends the PLS approach implemented in the R package **Morpho** by:

1. Extracting the first PLS axis between two morphological modules (e.g., cranial and cervical vertebrae) using `Morpho::pls2B`.  
2. Computing species-specific covariance matrices of the PLS scores (2D space).  
3. Identifying the main axis of covariation (first eigenvector) for each species.  
4. Quantifying the difference between species as the angle between their principal covariation axes.  
5. Assessing the significance of these differences via permutation tests.

## Installation

To use this function, download or clone the repository and source the main script in R:

```R
source("Integration_Patterns.R")
```

This script depends on:

- `Morpho` (for PLS); install with `install.packages("Morpho")`

## Usage Example (Simulated 3D Coordinates Data with Five Species)

```R
library(Morpho)
set.seed(42)

simulate_landmarks <- function(n_ind, n_land = 30) {
  array(rnorm(n_land * 3 * n_ind), dim = c(n_land, 3, n_ind))
}

# Simulate two modules for 100 individuals
X <- simulate_landmarks(100)
Y <- simulate_landmarks(100)

# Define species vector (five species with 20 individuals each)
species <- rep(c(
  "Homo_sapiens",
  "Pan_troglodytes",
  "Gorilla_gorilla",
  "Pongo_pygmaeus",
  "Hylobates_lar"
), each = 20)

# Run the analysis
res <- Integration_Patterns(
  X = X,
  Y = Y,
  species = species,
  rounds_pls = 0,    # no permutation inside pls2B
  n_perm = 200,      # number of permutations for angle comparison
  min_n = 3          # minimum individuals per species
)

# View angle matrix and p-value matrix
print(res$angle_matrix)
print(res$pval_matrix)
```

## Output

The function returns a list containing:

- `pls`: the Morpho PLS object.  
- `scores`: a dataframe of PLS scores for each individual along with species labels.  
- `data_by_species`: a list of score matrices by species.  
- `cov_matrices`: within-species covariance matrices.  
- `eigen_results`: eigenvectors and eigenvalues for each species.  
- `angle_matrix`: matrix of angles (degrees) between species' principal covariation axes.  
- `pval_matrix`: permutation-based significance matrix for the angle differences.  

## Citation

If you use this software, please cite it as:

Arlegi, M. (2025). *Integration Patterns: Quantifying between-species differences in morphological integration using two-block PLS* (Version 1.0.1). Zenodo. https://doi.org/10.5281/zenodo.17623869

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17623869.svg)](https://doi.org/10.5281/zenodo.17623869)
```

## Licence

This code is released under the GPL-3 License. See the `LICENSE` file for details.
