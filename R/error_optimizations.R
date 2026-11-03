#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Population Error Optimization Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Known graph (constrained maximum likelihood) via GLASSO
# Updated 03.11.2026
known_graph <- function(R, network, tol = 1e-08, max.iter = 1000)
{

  # Identify constrained (non-edge) cells in the upper triangle
  zeros <- which(network == 0 & upper.tri(network), arr.ind = TRUE)

  # Obtain GLASSO output
  output <- suppressWarnings(
    glasso::glasso(
      s = R, rho = 0, zero = zeros, thr = tol,
      maxit = max.iter, trace = FALSE, penalize.diagonal = FALSE
    )
  )

  # Set network
  network <- -cov2cor(output$wi); diag(network) <- 0

  # Return results
  return(
    list(
      network = network,
      W = output$w, Theta = output$wi,
      iterations = output$niter,
      converged = output$del < tol
    )
  )

}

#' @noRd
# General error vector function
# Updated 03.11.2026
error_function <- function(error_vector, R, R_inverse, log_det_R, free_cells, nodes, df, n = 8.3e09)
{

  # Set up error matrix (only free cells)
  error_matrix <- matrix(0, nodes, nodes)
  error_matrix[free_cells] <- error_vector

  # Population correlation matrix with error
  R_population <- R + (error_matrix + t(error_matrix))

  # Check positive definiteness
  population_chol <- try(chol(R_population), silent = TRUE)
  if(is(population_chol, "try-error")){
    return(list(valid = FALSE))
  }

  # Inverse and log-determinant of population matrix
  R_population_inverse <- chol2inv(population_chol)
  log_det_population <- 2 * sum(log(diag(population_chol)))

  # Minimum discrepancy (ML) and chi-square
  #   F = log|R| + tr(Sigma_pop R^-1) - log|Sigma_pop| - p
  discrepancy <- log_det_R + sum(R_population * R_inverse) - log_det_population - nodes
  chi_square <- n * discrepancy

  # Baseline chi-square (independence model): -log|Sigma_pop|
  baseline_chi_square <- -n * log_det_population
  denominator <- max(baseline_chi_square - nodes * (nodes - 1) / 2, 0)

  # Guard against non-finite results
  if(!is.finite(chi_square) || !is.finite(denominator) || denominator == 0){
    return(list(valid = FALSE))
  }

  # Compute CFI
  CFI <- 1 - (max(chi_square - df, 0) / denominator)

  # Return fit
  return(
    list(
      valid = TRUE, R_population = R_population,
      R_population_inverse = R_population_inverse,
      chi_square = chi_square, denominator = denominator,
      CFI = CFI
    )
  )

}

#' @noRd
# General error cost function
# Updated 03.11.2026
error_cost <- function(error_vector, target_CFI, R, R_inverse, log_det_R, free_cells, nodes, df, n = 8.3e09)
{

  # Obtain fit
  fit <- error_function(error_vector, R, R_inverse, log_det_R, free_cells, nodes, df, n)

  # Check for invalid results
  if(!fit$valid){
    return(1000)
  }

  # Smooth cost
  return((fit$CFI - target_CFI)^2)

}

#' @noRd
# General error gradient function
# Updated 03.11.2026
error_gradient <- function(error_vector, target_CFI, R, R_inverse, log_det_R, free_cells, nodes, df, n = 8.3e09)
{

  # Obtain fit
  fit <- error_function(error_vector, R, R_inverse, log_det_R, free_cells, nodes, df, n)

  # Cost is flat when invalid or when CFI is capped at 1 (chi-square <= df)
  if(!fit$valid || fit$chi_square <= df){
    return(numeric(length(error_vector)))
  }

  # Parameter e enters Sigma_pop at (i,j) and (j,i), so
  #   d F / de = 2 * (R^-1 - Sigma_pop^-1)_ij
  #   d baseline / de = n * -2 * (Sigma_pop^-1)_ij
  # CFI = 1 - (chi_square - df) / denominator
  CFI_gradient <- -(n / fit$denominator) * 2 * (R_inverse - fit$R_population_inverse) -
                   (fit$chi_square - df) / fit$denominator^2 * n * 2 * fit$R_population_inverse

  # d cost / de = 2 * (CFI - target) * d CFI / de
  return(2 * (fit$CFI - target_CFI) * CFI_gradient[free_cells])

}