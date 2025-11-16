###############################################
# Integration Patterns: Comparing morphological 
# integration patterns across species using 
# two-block PLS (Morpho::pls2B)
# Author: Mikel Arlegi
# License: GPL-3
###############################################

library(Morpho)

###############################################
# Helper function: angle between principal axes 
###############################################
Integration_Patterns_angle <- function(cov_matrix1, cov_matrix2) { 
  
  if (!all(dim(cov_matrix1) == c(2, 2)) || !all(dim(cov_matrix2) == c(2, 2))) {
    stop("Covariance matrices must be 2x2.")
  }
  
  eigen1 <- eigen(cov_matrix1)
  eigen2 <- eigen(cov_matrix2)
  
  v1 <- eigen1$vectors[, 1]
  v2 <- eigen2$vectors[, 1]
  
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- v2 / sqrt(sum(v2^2))
  
  cosine_similarity <- sum(v1 * v2)
  cosine_similarity <- max(min(cosine_similarity, 1), -1)
  
  angle <- acos(abs(cosine_similarity)) * (180 / pi)
  
  list(
    angle_degrees = angle,
    cos_similarity = cosine_similarity,
    eigenvector1  = v1,
    eigenvector2  = v2,
    eigenvalues1  = eigen1$values,
    eigenvalues2  = eigen2$values
  )
}

###############################################
# Permutation test for the angle between species
###############################################
Integration_Patterns_permutation <- function(data1, data2, n_perm = 1000) {
  
  if (nrow(data1) < 3 || nrow(data2) < 3) {
    warning("Low sample size in one species. Permutation test might be unreliable.")
  }
  
  cov1_obs <- cov(data1)
  cov2_obs <- cov(data2)
  res_obs  <- Integration_Patterns_angle(cov1_obs, cov2_obs)
  angle_obs <- res_obs$angle_degrees
  
  pooled <- rbind(data1, data2)
  n1 <- nrow(data1)
  n_total <- nrow(pooled)
  
  angles_perm <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    idx <- sample(n_total)
    group1 <- pooled[idx[1:n1], ]
    group2 <- pooled[idx[(n1 + 1):n_total], ]
    
    angles_perm[i] <- Integration_Patterns_angle(cov(group1), cov(group2))$angle_degrees
  }
  
  p_value <- mean(angles_perm >= angle_obs)
  
  list(
    angle_observed = angle_obs,
    angles_perm = angles_perm,
    p_value = p_value
  )
}

###############################################
# Main pipeline function
###############################################
Integration_Patterns <- function(
  X, Y,
  species,
  rounds_pls = 0,
  n_perm = 1000,
  min_n = 3
) {
  
  if (length(dim(X)) != 3 || length(dim(Y)) != 3) {
    stop("X and Y must be 3D arrays (landmarks x dimensions x individuals).")
  }
  
  if (dim(X)[3] != dim(Y)[3]) {
    stop("X and Y must contain the same individuals in the same order.")
  }
  
  n_ind <- dim(X)[3]
  
  if (length(species) != n_ind) {
    stop("Length of species vector must match number of individuals.")
  }
  
  pls_res <- Morpho::pls2B(X, Y, rounds = rounds_pls)
  
  scores_X <- pls_res$Xscores[, 1]
  scores_Y <- pls_res$Yscores[, 1]
  
  df <- data.frame(
    species = species,
    score_X = scores_X,
    score_Y = scores_Y,
    stringsAsFactors = FALSE
  )
  
  species_levels <- sort(unique(species))
  
  data_by_species <- list()
  cov_matrices <- list()
  eigen_results <- list()
  
  for (sp in species_levels) {
    sp_data <- subset(df, species == sp, select = c("score_X", "score_Y"))
    m <- as.matrix(sp_data)
    data_by_species[[sp]] <- m
    
    if (nrow(m) < min_n) {
      cov_matrices[[sp]] <- NA
      eigen_results[[sp]] <- NA
    } else {
      cov_m <- cov(m)
      cov_matrices[[sp]] <- cov_m
      e <- eigen(cov_m)
      v1 <- e$vectors[, 1] / sqrt(sum(e$vectors[, 1]^2))
      eigen_results[[sp]] <- list(
        eigenvector1 = v1,
        eigenvalues = e$values
      )
    }
  }
  
  n_sp <- length(species_levels)
  angle_matrix <- matrix(NA, n_sp, n_sp, dimnames = list(species_levels, species_levels))
  pval_matrix  <- matrix(NA, n_sp, n_sp, dimnames = list(species_levels, species_levels))
  
  for (i in seq_len(n_sp)) {
    for (j in seq_len(n_sp)) {
      sp_i <- species_levels[i]
      sp_j <- species_levels[j]
      
      if (i == j) {
        angle_matrix[i, j] <- 0
        pval_matrix[i, j]  <- NA
        next
      }
      
      if (is.na(cov_matrices[[sp_i]])[1] || is.na(cov_matrices[[sp_j]])[1]) {
        angle_matrix[i, j] <- NA
        pval_matrix[i, j]  <- NA
        next
      }
      
      angle_res <- Integration_Patterns_angle(cov_matrices[[sp_i]], cov_matrices[[sp_j]])
      angle_matrix[i, j] <- angle_res$angle_degrees
      
      perm_res <- Integration_Patterns_permutation(
        data_by_species[[sp_i]],
        data_by_species[[sp_j]],
        n_perm = n_perm
      )
      pval_matrix[i, j] <- perm_res$p_value
    }
  }
  
  angle_matrix[lower.tri(angle_matrix)] <- t(angle_matrix)[lower.tri(angle_matrix)]
  pval_matrix[lower.tri(pval_matrix)]  <- t(pval_matrix)[lower.tri(pval_matrix)]
  
  list(
    pls  = pls_res,
    scores = df,
    data_by_species = data_by_species,
    cov_matrices = cov_matrices,
    eigen_results = eigen_results,
    angle_matrix = angle_matrix,
    pval_matrix  = pval_matrix,
    params = list(rounds_pls = rounds_pls, n_perm = n_perm, min_n = min_n)
  )
}
