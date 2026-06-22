###############################################
# Integration_Patterns
# Comparing morphological integration patterns across groups
# Author: Mikel Arlegi
# License: GPL-3
###############################################

#' Angular comparison of group-specific integration patterns from two-block PLS
#'
#' Performs a comparative analysis of morphological integration patterns among groups
#' (e.g., species) using two-block partial least squares (PLS). The function estimates
#' group-specific axes of covariation in a shared PLS space and quantifies their
#' divergence using angular distances and permutation-based significance tests.
#'
#' The procedure first performs a pooled two-block PLS (\code{Morpho::pls2B}) on two
#' matched data blocks (\code{X} and \code{Y}), which may be Procrustes shape coordinates
#' (3D arrays, p x k x n) or other multivariate descriptors (matrices, n x variables).
#' Scores from a selected PLS axis are extracted for each block and used to define a
#' common PLS score space. By default, the first PLS axis is used
#' (\code{pls_axis = 1}), corresponding to the dominant mode of between-block
#' covariation. For each group, the covariance matrix of the selected PLS scores is
#' computed and its first eigenvector is taken as the dominant axis of covariation.
#' Pairwise angular distances between these axes represent differences in integration
#' patterns among groups for the selected PLS mode.
#'
#' By default, the function performs pooled within-group PLS (\code{pooled_within = TRUE}),
#' in which group means are removed prior to PLS estimation. This configuration removes
#' between-group mean differences and focuses inference on within-group covariance
#' structure. Alternatively, raw (non-centred) PLS can be performed
#' (\code{pooled_within = FALSE}), so the PLS space reflects combined within- and
#' between-group variation.
#'
#' Specimen correspondence between \code{X} and \code{Y} is required. If specimen
#' identifiers are provided (dimnames on the 3rd dimension for arrays; rownames for
#' matrices), \code{Y} is automatically reordered to match \code{X}. If identifiers are
#' absent, the function assumes both blocks are already in identical order.
#'
#' Angular differences are assessed via permutation tests that pool individuals from
#' each pair of groups and randomly reallocate them while preserving sample sizes,
#' generating a null distribution under the hypothesis of no group-specific difference
#' in integration axis orientation.
#'
#' @param X A 3D array (p x k x n) of Procrustes coordinates for the first module, or a
#'   matrix (n x variables) of other multivariate descriptors (rows = specimens).
#' @param Y A 3D array (p x k x n) for the second module, or a matrix (n x variables).
#'   Must contain the same specimens as \code{X}. If specimen identifiers are present,
#'   \code{Y} is reordered to match \code{X}.
#' @param groups A factor or character vector of length \code{n} indicating group
#'   membership for each specimen.
#' @param pooled_within Logical. If \code{TRUE} (default), performs pooled within-group
#'   PLS by centring variables within groups prior to PLS.
#' @param pls_axis Integer. PLS axis to analyse. Default is \code{1}, corresponding to
#'   the dominant PLS axis.
#' @param n_perm Number of permutations used to assess angular differences.
#' @param min_n Minimum sample size required per group to estimate covariance matrices.
#'
#' @details
#' The function implements a geometric comparison of integration patterns by evaluating
#' the orientation of group-specific covariation axes in the score space of a selected
#' PLS axis. For each group, the covariance matrix of the selected X and Y PLS scores is
#' computed and its first eigenvector is taken as the principal axis of covariation.
#' Because eigenvectors represent axes rather than oriented vectors, sign is ignored
#' and angular distances are computed as acute angles in the interval 0–90 degrees.
#'
#' Two analytical modes are available:
#' \itemize{
#'   \item When \code{pooled_within = TRUE} (default), group means are removed prior to PLS
#'   estimation. This produces a pooled within-group PLS space reflecting the average
#'   within-group covariation pattern, and angular comparisons test whether within-group
#'   integration patterns are conserved across taxa independently of mean divergence.
#'   \item When \code{pooled_within = FALSE}, PLS is performed on raw data without group
#'   centring. The resulting space reflects combined within- and between-group variation,
#'   and angular comparisons quantify differences in group-specific integration axes within
#'   this mixed-variance space.
#' }
#'
#' The function also returns diagnostics for the selected PLS axis. These include the
#' singular value associated with the axis, the percentage of total between-block
#' covariation accounted for by that axis, and the correlation between the corresponding
#' X and Y PLS scores. These diagnostics help assess whether the selected PLS axis
#' represents a sufficiently dominant and interpretable mode of covariation for angular
#' comparison.
#'
#' @return An object of class \code{"Integration_Patterns"} containing:
#'   \itemize{
#'     \item \code{pls}: The PLS result returned by \code{Morpho::pls2B}.
#'     \item \code{pls_diagnostics}: Data frame reporting the selected PLS axis, its
#'     singular value, the percentage of total between-block covariation accounted for
#'     by that axis, and the correlation between the corresponding X and Y PLS scores.
#'     \item \code{scores}: Data frame of selected PLS scores with group labels.
#'     \item \code{data_by_groups}: List of selected PLS score matrices per group.
#'     \item \code{cov_matrices}: Within-group covariance matrices in the selected PLS
#'     score space.
#'     \item \code{eigen_results}: Eigenvectors and eigenvalues for each group.
#'     \item \code{angle_matrix}: Matrix of pairwise angular distances (degrees).
#'     \item \code{pval_matrix}: Matrix of permutation-based p-values.
#'     \item \code{params}: List of function parameters used.
#'     \item \code{call}: The matched call.
#'   }
#'
#' @export
#' @keywords multivariate morphometrics integration PLS
#'
#' @references
#' Rohlf, F.J. & Corti, M. (2000). The use of two-block partial least-squares
#' to study covariation in shape. Systematic Biology 49:740–753.
#'
#' Klingenberg, C.P. & Marugán-Lobón, J. (2013). Evolutionary covariation in
#' geometric morphometric data. Systematic Biology 62:591–610.
#'
#' Schlager, S. (2017). Morpho and Rvcg – Shape analysis in R.
#' In: Zheng, G., Li, S., Székely, G. (eds) Statistical Shape and Deformation Analysis.
#' Academic Press.
#'
#' @author Mikel Arlegi
#'
#' @examples
#' \dontrun{
#' Hominoids <- readRDS("Hominoids.rds")
#'
#' res <- Integration_Patterns(
#'   X = Hominoids$crania,
#'   Y = Hominoids$C5_vertebrae,
#'   groups = Hominoids$species,
#'   pooled_within = TRUE,
#'   pls_axis = 1,
#'   n_perm = 999
#' )
#'
#' res$pls_diagnostics
#' res$angle_matrix
#' res$pval_matrix
#' }

Integration_Patterns <- function(
    X, Y,
    groups,
    pooled_within = TRUE,
    pls_axis = 1,
    n_perm = 999,
    min_n = 3
) {
  
  if (!requireNamespace("Morpho", quietly = TRUE)) {
    stop("Package 'Morpho' is required but not installed.", call. = FALSE)
  }
  
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  
  is_arrayX <- (length(dim(X)) == 3)
  is_arrayY <- (length(dim(Y)) == 3)
  is_matX   <- is.matrix(X)
  is_matY   <- is.matrix(Y)
  
  if (!(is_arrayX || is_matX) || !(is_arrayY || is_matY)) {
    stop("X and Y must be either 3D arrays (p x k x n) or matrices (n x variables).", call. = FALSE)
  }
  
  n_ind_X <- if (is_arrayX) dim(X)[3] else nrow(X)
  n_ind_Y <- if (is_arrayY) dim(Y)[3] else nrow(Y)
  
  if (n_ind_X != n_ind_Y) {
    stop("X and Y must contain the same number of specimens (same n).", call. = FALSE)
  }
  
  n_ind <- n_ind_X
  if (length(groups) != n_ind) {
    stop("Length of groups must match the number of specimens (n).", call. = FALSE)
  }
  
  if (any(is.na(groups))) {
    stop("'groups' contains missing values.", call. = FALSE)
  }
  
  if (length(unique(groups)) < 2) {
    stop("'groups' must contain at least two groups.", call. = FALSE)
  }
  
  if (length(pls_axis) != 1 || is.na(pls_axis) || pls_axis < 1 || pls_axis != as.integer(pls_axis)) {
    stop("'pls_axis' must be a single positive integer.", call. = FALSE)
  }
  
  groups <- as.character(groups)
  group_levels <- sort(unique(groups))
  
  ### Reorder by specimen identifiers
  if (is_arrayX && !is.null(dimnames(X)[[3]]) && !is.null(dimnames(Y)[[3]])) {
    
    id_x <- trimws(dimnames(X)[[3]])
    id_y <- trimws(dimnames(Y)[[3]])
    
    if (anyDuplicated(id_x) > 0 || anyDuplicated(id_y) > 0) {
      stop("Duplicate specimen identifiers detected in X or Y (3rd dimension).", call. = FALSE)
    }
    if (!setequal(id_x, id_y)) {
      stop("Mismatched specimen identifiers between X and Y (3rd dimension).", call. = FALSE)
    }
    
    Y <- Y[, , match(id_x, id_y), drop = FALSE]
    dimnames(Y)[[3]] <- dimnames(X)[[3]]
    
  } else if (is_matX && !is.null(rownames(X)) && !is.null(rownames(Y))) {
    
    id_x <- trimws(rownames(X))
    id_y <- trimws(rownames(Y))
    
    if (anyDuplicated(id_x) > 0 || anyDuplicated(id_y) > 0) {
      stop("Duplicate specimen identifiers detected in X or Y (rownames).", call. = FALSE)
    }
    if (!setequal(id_x, id_y)) {
      stop("Mismatched specimen identifiers between X and Y (rownames).", call. = FALSE)
    }
    
    Y <- Y[match(id_x, id_y), , drop = FALSE]
    rownames(Y) <- rownames(X)
  }
  
  X_use <- X
  Y_use <- Y
  
  ### Pooled within-group
  if (isTRUE(pooled_within)) {
    
    for (gr in group_levels) {
      idx <- which(groups == gr)
      if (length(idx) == 0) next
      
      if (is_arrayX) {
        mean_X <- apply(X[,,idx, drop = FALSE], c(1, 2), mean)
        for (ii in idx) X_use[,,ii] <- X[,,ii] - mean_X
      } else {
        mu_X <- colMeans(X[idx, , drop = FALSE])
        X_use[idx, ] <- sweep(X[idx, , drop = FALSE], 2, mu_X, "-")
      }
      
      if (is_arrayY) {
        mean_Y <- apply(Y[,,idx, drop = FALSE], c(1, 2), mean)
        for (ii in idx) Y_use[,,ii] <- Y[,,ii] - mean_Y
      } else {
        mu_Y <- colMeans(Y[idx, , drop = FALSE])
        Y_use[idx, ] <- sweep(Y[idx, , drop = FALSE], 2, mu_Y, "-")
      }
    }
  }
  
  pls_res <- Morpho::pls2B(X_use, Y_use, rounds = 0)
  
  ### Check selected PLS axis
  if (pls_axis > ncol(pls_res$Xscores) || pls_axis > ncol(pls_res$Yscores)) {
    stop("'pls_axis' must correspond to an available PLS axis.", call. = FALSE)
  }
  
  ### PLS diagnostics for the selected axis
  get_pls_diagnostics <- function(pls_res, X_use, Y_use, pls_axis) {
    
    singular_values <- NULL
    
    if (!is.null(pls_res$svd) && !is.null(pls_res$svd$d)) {
      singular_values <- pls_res$svd$d
    } else if (!is.null(pls_res$singular.values)) {
      singular_values <- pls_res$singular.values
    } else if (!is.null(pls_res$sv)) {
      singular_values <- pls_res$sv
    }
    
    ### Fallback: compute singular values directly from the cross-covariance matrix
    if (is.null(singular_values)) {
      
      X_mat <- X_use
      Y_mat <- Y_use
      
      if (length(dim(X_mat)) == 3) {
        X_mat <- Morpho::two.d.array(X_mat)
      }
      if (length(dim(Y_mat)) == 3) {
        Y_mat <- Morpho::two.d.array(Y_mat)
      }
      
      Xc <- scale(X_mat, center = TRUE, scale = FALSE)
      Yc <- scale(Y_mat, center = TRUE, scale = FALSE)
      
      Cxy <- crossprod(Xc, Yc) / (nrow(Xc) - 1)
      singular_values <- svd(Cxy)$d
    }
    
    if (pls_axis > length(singular_values)) {
      singular_value_axis <- NA_real_
      percent_total_covariation <- NA_real_
    } else {
      singular_value_axis <- singular_values[pls_axis]
      percent_total_covariation <- 100 * singular_value_axis^2 / sum(singular_values^2)
    }
    
    score_correlation <- cor(
      pls_res$Xscores[, pls_axis],
      pls_res$Yscores[, pls_axis],
      use = "complete.obs"
    )
    
    data.frame(
      pls_axis = pls_axis,
      singular_value = singular_value_axis,
      percent_total_covariation = percent_total_covariation,
      score_correlation = score_correlation
    )
  }
  
  pls_diagnostics <- get_pls_diagnostics(
    pls_res = pls_res,
    X_use = X_use,
    Y_use = Y_use,
    pls_axis = pls_axis
  )
  
  df <- data.frame(
    groups = groups,
    score_X = pls_res$Xscores[, pls_axis],
    score_Y = pls_res$Yscores[, pls_axis],
    stringsAsFactors = FALSE
  )
  
  data_by_groups <- vector("list", length(group_levels))
  cov_matrices    <- vector("list", length(group_levels))
  eigen_results   <- vector("list", length(group_levels))
  names(data_by_groups) <- names(cov_matrices) <- names(eigen_results) <- group_levels
  
  for (gr in group_levels) {
    
    idx <- which(df$groups == gr)
    m <- as.matrix(df[idx, c("score_X", "score_Y"), drop = FALSE])
    data_by_groups[[gr]] <- m
    
    if (nrow(m) < min_n) {
      cov_matrices[[gr]]  <- NA
      eigen_results[[gr]] <- NA
    } else {
      cov_m <- cov(m)
      cov_matrices[[gr]] <- cov_m
      
      e <- eigen(cov_m)
      v1 <- e$vectors[, 1]
      v1 <- v1 / sqrt(sum(v1^2))
      
      eigen_results[[gr]] <- list(
        eigenvector1 = v1,
        eigenvalues  = e$values
      )
    }
  }
  
  n_gr <- length(group_levels)
  angle_matrix <- matrix(NA_real_, n_gr, n_gr, dimnames = list(group_levels, group_levels))
  pval_matrix  <- matrix(NA_real_, n_gr, n_gr, dimnames = list(group_levels, group_levels))
  diag(angle_matrix) <- 0
  
  for (i in seq_len(n_gr - 1)) {
    
    gr_i <- group_levels[i]
    if (is.na(cov_matrices[[gr_i]])[1]) next
    
    for (j in (i + 1):n_gr) {
      
      gr_j <- group_levels[j]
      if (is.na(cov_matrices[[gr_j]])[1]) next
      
      cov1 <- cov_matrices[[gr_i]]
      cov2 <- cov_matrices[[gr_j]]
      
      e1 <- eigen(cov1)
      e2 <- eigen(cov2)
      
      v1 <- e1$vectors[, 1]
      v2 <- e2$vectors[, 1]
      v1 <- v1 / sqrt(sum(v1^2))
      v2 <- v2 / sqrt(sum(v2^2))
      
      cs <- sum(v1 * v2)
      cs <- max(min(cs, 1), -1)
      angle_obs <- acos(abs(cs)) * (180 / pi)
      
      angle_matrix[i, j] <- angle_obs
      
      data1 <- data_by_groups[[gr_i]]
      data2 <- data_by_groups[[gr_j]]
      
      pooled <- rbind(data1, data2)
      n1 <- nrow(data1)
      n_total <- nrow(pooled)
      
      angles_perm <- numeric(n_perm)
      
      for (pp in seq_len(n_perm)) {
        idxp <- sample.int(n_total)
        
        g1 <- pooled[idxp[1:n1], , drop = FALSE]
        g2 <- pooled[idxp[(n1 + 1):n_total], , drop = FALSE]
        
        eg1 <- eigen(cov(g1))
        eg2 <- eigen(cov(g2))
        
        vg1 <- eg1$vectors[, 1]
        vg2 <- eg2$vectors[, 1]
        
        vg1 <- vg1 / sqrt(sum(vg1^2))
        vg2 <- vg2 / sqrt(sum(vg2^2))
        
        cs_p <- sum(vg1 * vg2)
        cs_p <- max(min(cs_p, 1), -1)
        
        angles_perm[pp] <- acos(abs(cs_p)) * (180 / pi)
      }
      
      pval_matrix[i, j] <- (1 + sum(angles_perm >= angle_obs)) / (1 + n_perm)
    }
  }
  
  angle_matrix[lower.tri(angle_matrix)] <- t(angle_matrix)[lower.tri(angle_matrix)]
  pval_matrix[lower.tri(pval_matrix)]   <- t(pval_matrix)[lower.tri(pval_matrix)]
  
  out <- list(
    pls = pls_res,
    pls_diagnostics = pls_diagnostics,
    scores = df,
    data_by_groups = data_by_groups,
    cov_matrices = cov_matrices,
    eigen_results = eigen_results,
    angle_matrix = angle_matrix,
    pval_matrix  = pval_matrix,
    params = list(
      pooled_within = pooled_within,
      pls_axis = pls_axis,
      n_perm = n_perm,
      min_n = min_n
    ),
    call = match.call()
  )
  
  class(out) <- "Integration_Patterns"
  out
}
