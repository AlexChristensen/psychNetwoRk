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
# Reference (model-implied) matrix with inverse and log-determinant
# Updated 04.11.2026
reference_model <- function(R)
{

  # Cholesky decomposition
  R_chol <- chol(R)

  # Return reference
  return(
    list(
      R = R, inverse = chol2inv(R_chol),
      log_det = 2 * sum(log(diag(R_chol)))
    )
  )

}

#' @noRd
# Fit indices of a model-implied matrix to a population matrix
# Updated 04.11.2026
fit_indices <- function(R_population, reference, df, n = population_size)
{

  # Obtain number of nodes
  nodes <- ncol(R_population)

  # Check positive definiteness
  population_chol <- try(chol(R_population), silent = TRUE)
  if(is(population_chol, "try-error")){
    return(list(valid = FALSE))
  }

  # Inverse and log-determinant of population matrix
  R_population_inverse <- chol2inv(population_chol)
  log_det_population <- 2 * sum(log(diag(population_chol)))

  # Minimum discrepancy (ML) and chi-square
  #   F = log|Sigma_model| + tr(Sigma_pop Sigma_model^-1) - log|Sigma_pop| - p
  discrepancy <- reference$log_det + sum(R_population * reference$inverse) - log_det_population - nodes
  chi_square <- n * discrepancy

  # Baseline chi-square (independence model): -log|Sigma_pop|
  baseline_chi_square <- -n * log_det_population
  baseline_df <- nodes * (nodes - 1) / 2

  # Guard against non-finite results
  if(!is.finite(chi_square) || !is.finite(baseline_chi_square)){
    return(list(valid = FALSE))
  }

  # Chi-square above degrees of freedom
  model_excess <- chi_square - df
  baseline_excess <- baseline_chi_square - baseline_df

  # Check for degenerate baseline
  denominator <- max(baseline_excess, model_excess, 0)
  baseline_ratio <- baseline_chi_square / baseline_df
  if(denominator == 0 || baseline_ratio <= 1){
    return(list(valid = FALSE))
  }

  # Standardized root mean square residual (off-diagonal lower triangle only)
  residual <- R_population - reference$R
  SRMR <- sqrt(mean(residual[lower.tri(residual)]^2))

  # Return fit
  return(
    list(
      valid = TRUE, R_population = R_population,
      R_population_inverse = R_population_inverse,
      discrepancy = discrepancy, chi_square = chi_square,
      baseline_chi_square = baseline_chi_square,
      baseline_df = baseline_df, df = df,
      CFI = 1 - max(model_excess, 0) / denominator,
      TLI = (baseline_ratio - chi_square / df) / (baseline_ratio - 1),
      RMSEA = sqrt(max(model_excess, 0) / (df * n)),
      SRMR = SRMR
    )
  )

}

#' @noRd
# Cudeck and Browne (1992) error function
# Updated 04.11.2026
CB_error <- function(error_vector, reference, free_cells, nodes, df, n = population_size)
{

  # Set up error matrix (only free cells)
  error_matrix <- matrix(0, nodes, nodes)
  error_matrix[free_cells] <- error_vector

  # Population correlation matrix with error
  R_population <- reference$R + (error_matrix + t(error_matrix))

  # Return fit
  return(fit_indices(R_population, reference, df, n))

}

#' @noRd
# Cudeck and Browne (1992) cost function
# Updated 04.11.2026
CB_error_cost <- function(error_vector, metric, target, reference, free_cells, nodes, df, n = population_size)
{

  # Obtain fit
  fit <- CB_error(error_vector, reference, free_cells, nodes, df, n)

  # Check for invalid results
  if(!fit$valid){
    return(1000)
  }

  # Smooth cost
  return((fit[[metric]] - target)^2)

}

#' @noRd
# Cudeck and Browne (1992) gradient function
# Updated 04.11.2026
CB_error_gradient <- function(error_vector, metric, target, reference, free_cells, nodes, df, n = population_size)
{

  # Obtain fit
  fit <- CB_error(error_vector, reference, free_cells, nodes, df, n)

  # Cost is flat when invalid
  if(!fit$valid){
    return(numeric(length(error_vector)))
  }

  # Parameter e enters Sigma_pop at (i,j) and (j,i), so
  #   d chi-square / de = 2 * n * (R^-1 - Sigma_pop^-1)_ij
  #   d baseline / de = -2 * n * (Sigma_pop^-1)_ij
  chi_gradient <- 2 * n * (reference$inverse - fit$R_population_inverse)
  baseline_gradient <- -2 * n * fit$R_population_inverse

  # Chi-square above degrees of freedom
  model_excess <- fit$chi_square - df
  baseline_excess <- fit$baseline_chi_square - fit$baseline_df

  # Gradient of the metric
  metric_gradient <- switch(
    metric,
    "CFI" = {

      # Flat when chi-square is below degrees of freedom (CFI = 1)
      # or when it exceeds the baseline (CFI = 0)
      if(model_excess <= 0 || model_excess >= baseline_excess){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        -(chi_gradient * baseline_excess - model_excess * baseline_gradient) /
          baseline_excess^2
      }

    },
    "TLI" = {

      # B = baseline chi-square / df, M = model chi-square / df
      B <- fit$baseline_chi_square / fit$baseline_df
      M <- fit$chi_square / df
      (-(chi_gradient / df) * (B - 1) + (baseline_gradient / fit$baseline_df) * (M - 1)) /
        (B - 1)^2

    },
    "RMSEA" = {

      # Flat when chi-square is below degrees of freedom (RMSEA = 0)
      if(model_excess <= 0){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        chi_gradient / (2 * df * n * fit$RMSEA)
      }

    },
    "SRMR" = {

      # Flat at zero residual
      if(fit$SRMR == 0){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        (fit$R_population - reference$R) /
          (fit$SRMR * nodes * (nodes - 1) / 2)
      }

    }
  )

  # d cost / de = 2 * (metric - target) * d metric / de
  return(2 * (fit[[metric]] - target) * metric_gradient[free_cells])

}

#' @noRd
# Cudeck and Browne (1992) optimization
# Updated 04.11.2026
CB_optimize <- function(
    metric, target, reference, free_cells, nodes, df,
    n = population_size, tolerance = 0.01
)
{

  # Random direction over free cells
  # (random signs with magnitudes in [0.5, 1] keep the largest error close to
  # the typical error so that more of the feasible range can be used)
  direction <- shuffle_replace(c(-1, 1), length(free_cells)) *
    runif_xoshiro(length(free_cells), min = 0.5, max = 1)

  # Largest scale that keeps the errors within [-1, 1]
  maximum_scale <- 1 / max(abs(direction))

  # Largest scale that keeps the population matrix positive definite
  # (R + scale * E is positive definite when scale * lambda_max(-L^-1 E L^-T) < 1)
  error_matrix <- matrix(0, nodes, nodes)
  error_matrix[free_cells] <- direction
  error_matrix <- error_matrix + t(error_matrix)
  inverse_chol <- backsolve(chol(reference$R), diag(nodes))
  smallest <- min(
    eigen(
      crossprod(inverse_chol, error_matrix %*% inverse_chol),
      symmetric = TRUE, only.values = TRUE
    )$values
  )
  if(smallest < 0){
    maximum_scale <- min(maximum_scale, 0.999 * -1 / smallest)
  }

  # Scale the direction to reach the target metric
  scaled <- solve_scalar(
    FUN = function(scale){
      CB_error(direction * scale, reference, free_cells, nodes, df, n)
    },
    lower = 0, upper = maximum_scale,
    metric = metric, target = target
  )

  # Check whether the target was reached
  if(scaled$valid && abs(scaled[[metric]] - target) <= tolerance){
    return(scaled)
  }

  # Otherwise, optimize error directly (errors are bounded by [-1, 1])
  # starting from the closest valid scaled direction (zero is flat)
  start_scale <- swiftelse(
    scaled$valid, max(scaled$parameter, 0.01 * maximum_scale), 0.01 * maximum_scale
  )

  optimization <- try(
    optim(
      par = direction * start_scale,
      fn = CB_error_cost, gr = CB_error_gradient,
      metric = metric, target = target, reference = reference,
      free_cells = free_cells, nodes = nodes, df = df, n = n,
      method = "L-BFGS-B", lower = -1, upper = 1
    ), silent = TRUE
  )

  # Check for errors
  if(is(optimization, "try-error")){
    return(scaled)
  }

  # Obtain fit
  fit <- CB_error(optimization$par, reference, free_cells, nodes, df, n)

  # Return the solution closest to the target
  if(
    fit$valid && (
      !scaled$valid ||
      abs(fit[[metric]] - target) < abs(scaled[[metric]] - target)
    )
  ){
    return(fit)
  }else{
    return(scaled)
  }

}

#' @noRd
# Fit indices of the known graph (constrained maximum likelihood)
# Updated 04.11.2026
refit_indices <- function(R_population, network, df, n = population_size)
{

  # Obtain known graph
  known <- try(known_graph(R_population, network), silent = TRUE)

  # Check for errors
  if(is(known, "try-error")){
    return(list(valid = FALSE))
  }

  # Obtain fit of the known graph
  fit <- try(
    fit_indices(R_population, reference_model(known$W), df, n),
    silent = TRUE
  )

  # Check for errors
  if(is(fit, "try-error")){
    return(list(valid = FALSE))
  }

  # Return fit
  return(fit)

}

#' @noRd
# Solve for a scalar parameter that reaches the target metric
# Parameter increases with the amount of error (smallest error is preferred)
# Updated 04.11.2026
solve_scalar <- function(FUN, lower, upper, metric, target, grid_size = 20)
{

  # Difference between metric and target (missing when invalid)
  difference <- function(parameter){

    # Obtain fit
    fit <- FUN(parameter)

    # Return difference
    return(swiftelse(fit$valid, fit[[metric]] - target, NA))

  }

  # Evaluate grid of parameters
  grid <- seq(lower, upper, length.out = grid_size)
  values <- unlist(lapply(grid, difference))

  # Check for no valid values
  if(all(is.na(values))){
    return(list(valid = FALSE))
  }

  # Identify exact solutions and sign changes between valid values
  valid <- which(!is.na(values))
  exact <- valid[values[valid] == 0]
  changes <- which(diff(sign(values[valid])) != 0)

  # Choose the solution with the smallest parameter
  if(length(exact) > 0 && (length(changes) == 0 || grid[exact[1]] <= grid[valid[changes[1]]])){

    parameter <- grid[exact[1]]

  }else if(length(changes) > 0){

    # Bracket the smallest root
    bracket <- grid[valid[c(changes[1], changes[1] + 1)]]

    # Find root
    root <- try(
      uniroot(difference, interval = bracket, tol = 1e-10, maxiter = 200),
      silent = TRUE
    )

    # Use root or closest bracket
    parameter <- swiftelse(
      is(root, "try-error"),
      bracket[which.min(abs(values[valid[changes[1] + 0:1]]))],
      root$root
    )

  }else{ # target is not reachable so use closest parameter
    parameter <- grid[valid[which.min(abs(values[valid]))]]
  }

  # Return fit
  return(c(FUN(parameter), list(parameter = parameter)))

}

#' @noRd
# Tucker, Koopman, and Linn (1969) error function
# Updated 04.11.2026
TKL_error <- function(weight, minor, R, network, df, n = population_size)
{

  # Population correlation matrix with error
  # (minor factors contribute `weight` of the variance; diagonal remains one)
  R_population <- (1 - weight) * R + weight * minor

  # Return fit of known graph
  return(refit_indices(R_population, network, df, n))

}

#' @noRd
# Tucker, Koopman, and Linn (1969) optimization
# Updated 04.11.2026
TKL_optimize <- function(
    metric, target, R, network, nodes, df,
    minor_factors, decay, n = population_size
)
{

  # Random loadings on minor factors
  loadings <- matrix(
    rnorm_ziggurat(nodes * minor_factors),
    nrow = nodes, ncol = minor_factors
  )

  # Geometric decay of loadings across minor factors
  loadings <- loadings * rep((1 - decay)^(seq_len(minor_factors) - 1), each = nodes)

  # Scale each variable's loadings to unit length
  loadings <- loadings / sqrt(rowSums(loadings^2))

  # Correlations from minor factors
  minor <- tcrossprod(loadings)

  # Largest weight that keeps the errors within [-1, 1]
  # (error is linear in the weight: weight * (minor - R))
  maximum_weight <- min(
    0.999, 1 / max(abs((minor - R)[upper.tri(R)]))
  )

  # Return solution
  return(
    solve_scalar(
      FUN = function(weight){
        TKL_error(weight, minor, R, network, df, n)
      },
      lower = 0, upper = maximum_weight,
      metric = metric, target = target
    )
  )

}

#' @noRd
# Wu and Browne (2015) error function
# Updated 07.11.2026
WB_error <- function(log_dispersion, draws, R, network, df, n = population_size)
{

  # Obtain precision parameter (m = 1 / v)
  nodes <- ncol(R)
  precision <- exp(-log_dispersion)

  # Bartlett decomposition (fixed random draws, continuous precision)
  bartlett <- draws$normals
  diag(bartlett) <- sqrt(qchisq(draws$uniforms, df = precision - seq_len(nodes) + 1))

  # Wishart draw on the precision matrix, centered on the original precision
  # matrix (E[W] = R^-1), which is an inverted Wishart draw for the covariance
  #   Sigma = W^-1, W ~ Wishart(m, R^-1 / m)
  # Other centerings (e.g., E[Sigma] = R, which uses scale R^-1 / (m - p - 1))
  # differ only by a scalar multiple of W, which is removed when converted to
  # a correlation matrix, so the resulting draws are identical
  factor_inverse <- forwardsolve(draws$cholesky %*% bartlett, diag(nodes))
  R_draw <- cov2cor(crossprod(factor_inverse))

  # Scale error so that the errors are within [-1, 1]
  # (convex combination of positive definite matrices remains positive definite)
  error_matrix <- R_draw - R
  error_matrix <- error_matrix * min(1, 1 / max(abs(error_matrix)))

  # Population correlation matrix with error
  R_population <- R + error_matrix

  # Return fit of known graph
  return(refit_indices(R_population, network, df, n))

}

#' @noRd
# Wu and Browne (2015) optimization
# Updated 04.11.2026
WB_optimize <- function(metric, target, R, network, nodes, df, n = population_size)
{

  # Random draws (held fixed while the dispersion is optimized)
  normals <- matrix(0, nrow = nodes, ncol = nodes)
  normals[lower.tri(normals)] <- rnorm_ziggurat(nodes * (nodes - 1) / 2)
  draws <- list(
    normals = normals, uniforms = runif_xoshiro(nodes),
    cholesky = t(chol(chol2inv(chol(R))))
  )

  # Return solution
  # (dispersion v = 1 / m, with m > nodes - 1; v is approximately RMSEA^2)
  return(
    solve_scalar(
      FUN = function(log_dispersion){
        WB_error(log_dispersion, draws, R, network, df, n)
      },
      lower = log(1e-06), upper = log(1 / nodes),
      metric = metric, target = target
    )
  )

}

#' @noRd
# Hu and Bentler (1998) population correlation matrix
# Updated 04.11.2026
HB_population <- function(size, omissions)
{

  # Add omitted edges (partial correlations of `size`) to the precision matrix
  precision <- omissions$precision - size * omissions$shift

  # Check positive definiteness
  precision_chol <- try(chol(precision), silent = TRUE)
  if(is(precision_chol, "try-error")){
    return(NULL)
  }

  # Return population correlation matrix
  return(cov2cor(chol2inv(precision_chol)))

}

#' @noRd
# Hu and Bentler (1998) error function
# Updated 04.11.2026
HB_error <- function(size, omissions, network, df, n = population_size)
{

  # Population correlation matrix with error
  R_population <- HB_population(size, omissions)

  # Check for invalid results
  if(is.null(R_population)){
    return(list(valid = FALSE))
  }

  # Return fit of known graph
  return(refit_indices(R_population, network, df, n))

}

#' @noRd
# Hu and Bentler (1998) optimization
# Updated 04.11.2026
HB_optimize <- function(
    metric, target, R, network, nodes, df,
    free_cells, omitted, n = population_size
)
{

  # Randomly select omitted edges (at least one)
  omitted_cells <- shuffle(
    free_cells, size = max(1, round(omitted * length(free_cells)))
  )

  # Random signs for omitted edges
  omitted_edges <- matrix(0, nrow = nodes, ncol = nodes)
  omitted_edges[omitted_cells] <- swiftelse(
    runif_xoshiro(length(omitted_cells)) < 0.50, -1, 1
  )
  omitted_edges <- omitted_edges + t(omitted_edges)

  # Precision matrix of the original population and the shift in
  # the precision matrix that adds `size` to the partial correlations
  precision <- chol2inv(chol(R))
  omissions <- list(
    precision = precision,
    shift = omitted_edges * sqrt(tcrossprod(diag(precision)))
  )

  # Check whether size is allowed (positive definite and errors within [-1, 1])
  allowed <- function(size){

    # Obtain population correlation matrix
    R_population <- HB_population(size, omissions)

    # Return whether allowed
    return(!is.null(R_population) && max(abs(R_population - R)) <= 1)

  }

  # Find largest size (partial correlation of omitted edges)
  maximum_size <- 1
  if(!allowed(maximum_size)){

    # Bisection (smallest size is always allowed)
    lower <- 0; upper <- maximum_size
    for(i in seq_len(50)){

      # Update bounds
      middle <- (lower + upper) / 2
      if(allowed(middle)){
        lower <- middle
      }else{
        upper <- middle
      }

    }

    # Stay away from the boundary
    maximum_size <- lower * 0.999

  }

  # Return solution
  return(
    solve_scalar(
      FUN = function(size){
        HB_error(size, omissions, network, df, n)
      },
      lower = 0, upper = maximum_size,
      metric = metric, target = target
    )
  )

}
