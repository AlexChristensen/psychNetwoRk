#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Population Error Optimization Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# The true network (known graph) is fit to the population correlation matrix
# with error by constrained maximum likelihood (GLASSO with zero constraints
# and no penalty). The constrained maximum likelihood problem of a Gaussian
# graphical model is convex with a unique solution (Dempster, 1972). GLASSO
# is always started cold: warm starts from the solution of a different
# population matrix can stall GLASSO's inner coordinate descent (which is not
# bounded by `maxit`)
#
# Fit indices are population values computed from the maximum likelihood
# discrepancy F (no sample size):
#   RMSEA = sqrt(F / df)
#   CFI = 1 - F / F0, where F0 = -log|Sigma| (independence model)
#   TLI = 1 - (F / df) / (F0 / df0), where df0 = p (p - 1) / 2

#' @noRd
# Known graph (constrained maximum likelihood) via GLASSO
# Warm starts use a previous solution (`start = list(W, Theta)`)
# Updated 20.11.2026
known_graph <- function(R, network, tol = 1e-08, max.iter = 1000, start = NULL)
{

  # Identify constrained (non-edge) cells in the upper triangle
  zeros <- which(network == 0 & upper.tri(network), arr.ind = TRUE)

  # Obtain GLASSO output
  output <- suppressWarnings(
    if(is.null(start)){
      glasso::glasso(
        s = R, rho = 0, zero = zeros, thr = tol,
        maxit = max.iter, trace = FALSE, penalize.diagonal = FALSE
      )
    }else{
      glasso::glasso(
        s = R, rho = 0, zero = zeros, thr = tol,
        maxit = max.iter, trace = FALSE, penalize.diagonal = FALSE,
        start = "warm", w.init = start$W, wi.init = start$Theta
      )
    }
  )

  # Set network
  network <- -cov2cor(output$wi); diag(network) <- 0

  # Return results
  # (GLASSO stops when its own convergence criterion is met, so it converged
  # when it stopped before the maximum number of iterations)
  return(
    list(
      network = network,
      W = output$w, Theta = output$wi,
      iterations = output$niter,
      converged = output$niter < max.iter
    )
  )

}

#' @noRd
# Reference (model-implied) matrix with inverse and log-determinant
# Updated 20.11.2026
reference_model <- function(R, R_chol = chol(R))
{

  # Return reference
  return(
    list(
      R = R, inverse = chol2inv(R_chol),
      log_det = 2 * sum(log(diag(R_chol)))
    )
  )

}

#' @noRd
# Fit indices from the maximum likelihood discrepancy
#   F = log|Sigma_model| + tr(Sigma_pop Sigma_model^-1) - log|Sigma_pop| - p
# and the log-determinant of the population matrix (independence baseline
# F0 = -log|Sigma_pop|). Population (infinite sample size) indices:
#   RMSEA = sqrt(F / df)
#   CFI = 1 - F / F0
#   TLI = 1 - (F / df) / (F0 / df0), where df0 = p (p - 1) / 2
# (SRMR is from the residuals or, if supplied, `SRMR` directly)
# Updated 25.11.2026
fit_from_discrepancy <- function(
    discrepancy, log_det_population, residual = NULL, nodes, df,
    R_population = NULL, SRMR = NULL
)
{

  # Discrepancy of the model (zero within rounding error at no error)
  discrepancy <- max(discrepancy, 0)

  # Discrepancy of the independence model
  baseline <- -log_det_population
  baseline_df <- nodes * (nodes - 1) / 2

  # Guard against non-finite results and a degenerate baseline
  if(!is.finite(discrepancy) || !is.finite(baseline) || baseline <= 0){
    return(list(valid = FALSE))
  }

  # Standardized root mean square residual (lower triangle including the
  # diagonal, p * (p + 1) / 2 unique elements; same as {lavaan} for correlations)
  if(is.null(SRMR)){
    SRMR <- sqrt(mean(residual[lower.tri(residual, diag = TRUE)]^2))
  }

  # Return fit
  return(
    list(
      valid = TRUE, R_population = R_population,
      discrepancy = discrepancy, baseline = baseline,
      baseline_df = baseline_df, df = df,
      CFI = 1 - discrepancy / baseline,
      TLI = 1 - (discrepancy / df) / (baseline / baseline_df),
      RMSEA = sqrt(discrepancy / df), SRMR = SRMR
    )
  )

}

#' @noRd
# Fit indices of a model-implied matrix to a population matrix
# Updated 25.11.2026
fit_indices <- function(
    R_population, reference, df,
    population_chol = try(chol(R_population), silent = TRUE)
)
{

  # Check positive definiteness
  if(is(population_chol, "try-error")){
    return(list(valid = FALSE))
  }

  # Obtain number of nodes and log-determinant of population matrix
  nodes <- ncol(R_population)
  log_det_population <- 2 * sum(log(diag(population_chol)))

  # Return fit
  return(
    fit_from_discrepancy(
      discrepancy = reference$log_det + sum(R_population * reference$inverse) -
        log_det_population - nodes,
      log_det_population = log_det_population,
      residual = R_population - reference$R,
      nodes = nodes, df = df, R_population = R_population
    )
  )

}

#' @noRd
# Fit indices of the known graph (constrained maximum likelihood)
# Updated 25.11.2026
refit_indices <- function(R_population, network, df, start = NULL)
{

  # Check positive definiteness
  population_chol <- try(chol(R_population), silent = TRUE)
  if(is(population_chol, "try-error")){
    return(list(valid = FALSE))
  }

  # Obtain known graph
  known <- try(known_graph(R_population, network, start = start), silent = TRUE)

  # Check for errors and non-convergence
  # (a non-converged solution overstates misfit)
  if(is(known, "try-error") || !isTRUE(known$converged)){
    return(list(valid = FALSE))
  }

  # Obtain fit of the known graph
  reference <- try(reference_model(known$W), silent = TRUE)
  if(is(reference, "try-error")){
    return(list(valid = FALSE))
  }
  fit <- fit_indices(R_population, reference, df, population_chol = population_chol)

  # Add known graph
  if(fit$valid){
    fit$model <- known
  }

  # Return fit
  return(fit)

}

#' @noRd
# Refit function of the known graph
# Returns a function of the scalar parameter
# Updated 25.11.2026
scalar_refit <- function(population, network, df)
{

  # Return function
  return(
    function(parameter){

      # Population correlation matrix with error
      R_population <- population(parameter)
      if(is.null(R_population)){
        return(list(valid = FALSE))
      }

      # Return fit of known graph
      return(refit_indices(R_population, network, df))

    }
  )

}

#' @noRd
# Solve for a scalar parameter that reaches the target metric
# The parameter increases with the amount of error and `lower` corresponds to
# (essentially) no error, where the metric is at its best value (CFI/TLI near
# one, RMSEA/SRMR near zero). The smallest error that reaches the target is
# found by searching upward from `lower` (geometric or linear steps) to the
# first change of sign, then refining with Brent's method. Evaluations are
# memoized
# Updated 20.11.2026
solve_scalar <- function(
    FUN, lower, upper, metric, target, steps = 10, tolerance = 1e-10,
    spacing = c("geometric", "linear")
)
{

  # Set spacing of the search
  spacing <- match.arg(spacing)

  # Sign of the difference at no error
  # (CFI/TLI are above and RMSEA/SRMR are below their targets)
  zero_sign <- swiftelse(metric %in% c("CFI", "TLI"), 1, -1)

  # Memoized evaluations
  parameters <- numeric(0)
  fits <- list()
  evaluate <- function(parameter){

    # Check for previous evaluation
    index <- which(parameters == parameter)
    if(length(index) > 0){
      return(fits[[index[1]]])
    }

    # Evaluate
    fit <- FUN(parameter)
    parameters <<- c(parameters, parameter)
    fits[[length(fits) + 1]] <<- fit

    # Return fit
    return(fit)

  }

  # Difference between metric and target (missing when invalid)
  difference <- function(parameter){
    fit <- evaluate(parameter)
    return(swiftelse(isTRUE(fit$valid), fit[[metric]] - target, NA))
  }

  # Search upward from the smallest error: geometric steps
  # (lower + (upper - lower) / 2^k) or equally spaced steps
  grid <- switch(
    spacing,
    "geometric" = lower + (upper - lower) / 2^(rev(seq_len(steps)) - 1),
    "linear" = seq(lower, upper, length.out = steps + 1)[-1]
  )
  previous <- lower; previous_value <- zero_sign * abs(target)
  bracket <- NULL; first_invalid <- NULL
  for(point in grid){

    # Obtain difference
    value <- difference(point)

    # Skip invalid points (bracket only between valid points)
    if(is.na(value)){
      if(is.null(first_invalid)){first_invalid <- point}
      next
    }
    first_invalid <- NULL

    # Check for exact solution
    if(value == 0){
      return(c(evaluate(point), list(parameter = point)))
    }

    # Check for change of sign
    if(sign(value) != zero_sign){
      bracket <- c(previous, point)
      bracket_values <- c(previous_value, value)
      break
    }

    # Update previous valid point
    previous <- point; previous_value <- value

  }

  # The search ended in invalid points (e.g., errors outside of -1 to 1):
  # refine toward the boundary of valid errors looking for the target
  if(is.null(bracket) && !is.null(first_invalid)){

    # Interval between the last valid and the first invalid point
    a <- previous; b <- first_invalid
    for(i in seq_len(30)){

      # Middle point
      middle <- (a + b) / 2
      value <- difference(middle)

      # Update interval
      if(is.na(value)){
        b <- middle
      }else if(value == 0){
        return(c(evaluate(middle), list(parameter = middle)))
      }else if(sign(value) != zero_sign){
        bracket <- c(a, middle)
        bracket_values <- c(previous_value, value)
        break
      }else{
        a <- middle; previous_value <- value
      }

    }

  }

  # Target is not reachable so use the closest evaluated parameter
  if(is.null(bracket)){

    # Valid evaluations
    valid <- vapply(fits, function(x){isTRUE(x$valid)}, logical(1))
    if(!any(valid)){
      return(list(valid = FALSE))
    }

    # Closest
    distances <- vapply(
      fits[valid], function(x){abs(x[[metric]] - target)}, numeric(1)
    )
    parameter <- parameters[valid][which.min(distances)]
    return(c(evaluate(parameter), list(parameter = parameter)))

  }

  # Refine with Brent's method (invalid points inside the bracket are
  # treated as beyond the target so the bracket shrinks toward smaller error)
  root <- try(
    uniroot(
      function(parameter){
        value <- difference(parameter)
        swiftelse(is.na(value), bracket_values[2], value)
      },
      interval = bracket, f.lower = bracket_values[1], f.upper = bracket_values[2],
      tol = tolerance * max(1, abs(upper - lower)), maxiter = 200
    ), silent = TRUE
  )

  # Obtain parameter (closest valid bracket end if root finding failed)
  if(is(root, "try-error")){
    parameter <- bracket[2]
  }else{
    parameter <- root$root
    if(is.na(difference(parameter))){
      parameter <- bracket[2]
    }
  }

  # Return fit
  return(c(evaluate(parameter), list(parameter = parameter)))

}

#' @noRd
# Solve with refits of the known graph
# Updated 25.11.2026
solve_refit <- function(
    population, network, df, lower, upper, metric, target,
    spacing = "geometric"
)
{
  return(
    solve_scalar(
      scalar_refit(population, network, df), lower, upper,
      metric, target, spacing = spacing
    )
  )
}

# Limits on the population correlation matrix ----

#' @noRd
# Check the condition number of a population correlation matrix
# (GLASSO becomes very slow on nearly singular matrices because its inner
# coordinate descent is not bounded by `maxit`)
# Updated 25.11.2026
within_condition <- function(R_population, max_condition)
{

  # No limit
  if(!is.finite(max_condition)){
    return(TRUE)
  }

  # Check condition number
  eigenvalues <- matrix_eigenvalues(R_population)
  return(max(eigenvalues) / min(eigenvalues) <= max_condition)

}

#' @noRd
# Largest parameter (up to one) for which `population` returns a population
# correlation matrix (it returns `NULL` when not allowed). The smallest
# parameter is always allowed
# Updated 25.11.2026
largest_allowed <- function(population)
{

  # Check upper bound
  if(!is.null(population(1))){
    return(1)
  }

  # Bisection
  lower <- 0; upper <- 1
  for(i in seq_len(50)){

    # Update bounds
    middle <- (lower + upper) / 2
    if(!is.null(population(middle))){
      lower <- middle
    }else{
      upper <- middle
    }

  }

  # Stay away from the boundary
  return(lower * 0.999)

}

# Cudeck and Browne (1992) ----

#' @noRd
# Cudeck and Browne (1992) population correlation matrix
# Updated 20.11.2026
CB_population <- function(error_vector, reference, free_cells, nodes)
{

  # Set up error matrix (only free cells)
  error_matrix <- matrix(0, nodes, nodes)
  error_matrix[free_cells] <- error_vector

  # Population correlation matrix with error
  return(reference$R + (error_matrix + t(error_matrix)))

}

#' @noRd
# Cudeck and Browne (1992) error function
# (the known graph is the minimizer so the reference is the original matrix)
# Updated 25.11.2026
CB_error <- function(error_vector, reference, free_cells, nodes, df)
{
  return(
    fit_indices(CB_population(error_vector, reference, free_cells, nodes), reference, df)
  )
}

#' @noRd
# Cudeck and Browne (1992) closed-form fit along a direction
# For R_population = R + s D with eigenvalues l of U^-T D U^-1 (R = U'U):
#   log|R_population| = log|R| + sum(log(1 + s l))
#   F(s) = s sum(l) - sum(log(1 + s l))
#   SRMR(s) = s * SRMR(D)
# (sum(l) = tr(R^-1 D) = 0 because the precision matrix is zero on the cells
# where D is non-zero)
# Updated 25.11.2026
CB_direction_fit <- function(scale, eigenvalues, unit_SRMR, reference, nodes, df)
{

  # Positive definiteness
  shifted <- 1 + scale * eigenvalues
  if(any(shifted <= 0)){
    return(list(valid = FALSE))
  }

  # Return fit
  return(
    fit_from_discrepancy(
      discrepancy = scale * sum(eigenvalues) - sum(log(shifted)),
      log_det_population = reference$log_det + sum(log(shifted)),
      nodes = nodes, df = df, SRMR = scale * unit_SRMR
    )
  )

}

#' @noRd
# Cudeck and Browne (1992) cost function
# Updated 25.11.2026
CB_error_cost <- function(error_vector, metric, target, reference, free_cells, nodes, df)
{

  # Obtain fit
  fit <- CB_error(error_vector, reference, free_cells, nodes, df)

  # Check for invalid results
  if(!fit$valid){
    return(1000)
  }

  # Smooth cost
  return((fit[[metric]] - target)^2)

}

#' @noRd
# Cudeck and Browne (1992) gradient function
# Updated 25.11.2026
CB_error_gradient <- function(error_vector, metric, target, reference, free_cells, nodes, df)
{

  # Obtain fit
  fit <- CB_error(error_vector, reference, free_cells, nodes, df)

  # Cost is flat when invalid
  if(!fit$valid){
    return(numeric(length(error_vector)))
  }

  # Inverse of the population matrix
  R_population_inverse <- chol2inv(chol(fit$R_population))

  # Parameter e enters Sigma_pop at (i,j) and (j,i), so
  #   d F / de = 2 * (R^-1 - Sigma_pop^-1)_ij
  #   d F0 / de = -2 * (Sigma_pop^-1)_ij
  discrepancy_gradient <- 2 * (reference$inverse - R_population_inverse)
  baseline_gradient <- -2 * R_population_inverse

  # Gradient of the metric
  metric_gradient <- switch(
    metric,
    "CFI" = {

      # CFI = 1 - F / F0
      -(discrepancy_gradient * fit$baseline - fit$discrepancy * baseline_gradient) /
        fit$baseline^2

    },
    "TLI" = {

      # TLI = 1 - (F df0) / (F0 df)
      -(fit$baseline_df / df) *
        (discrepancy_gradient * fit$baseline - fit$discrepancy * baseline_gradient) /
        fit$baseline^2

    },
    "RMSEA" = {

      # Flat at no error (RMSEA = 0)
      if(fit$discrepancy <= 0){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        discrepancy_gradient / (2 * df * fit$RMSEA)
      }

    },
    "SRMR" = {

      # Flat at zero residual
      if(fit$SRMR == 0){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        (fit$R_population - reference$R) /
          (fit$SRMR * nodes * (nodes + 1) / 2)
      }

    }
  )

  # d cost / de = 2 * (metric - target) * d metric / de
  return(2 * (fit[[metric]] - target) * metric_gradient[free_cells])

}

#' @noRd
# Cudeck and Browne (1992) optimization
# Updated 25.11.2026
CB_optimize <- function(
    metric, target, reference, free_cells, nodes, df,
    tolerance = 0.01, inverse_chol = backsolve(chol(reference$R), diag(nodes))
)
{

  # Random direction over free cells
  # (random signs with magnitudes in [0.5, 1] keep the largest error close to
  # the typical error so that more of the feasible range can be used)
  direction <- shuffle_replace(c(-1, 1), length(free_cells)) *
    runif_xoshiro(length(free_cells), min = 0.5, max = 1)

  # Error matrix of the direction
  error_matrix <- CB_population(direction, reference, free_cells, nodes) - reference$R

  # Eigenvalues of U^-T D U^-1 (closed-form fit along the direction)
  eigenvalues <- eigen(
    crossprod(inverse_chol, error_matrix %*% inverse_chol),
    symmetric = TRUE, only.values = TRUE
  )$values
  unit_SRMR <- sqrt(mean(error_matrix[lower.tri(error_matrix, diag = TRUE)]^2))

  # Largest scale that keeps the errors within [-1, 1]
  maximum_scale <- 1 / max(abs(direction))

  # Largest scale that keeps the population matrix positive definite
  # (R + scale * D is positive definite when 1 + scale * lambda_min > 0)
  smallest <- min(eigenvalues)
  if(smallest < 0){
    maximum_scale <- min(maximum_scale, 0.999 * -1 / smallest)
  }

  # Scale the direction to reach the target metric (closed form)
  scaled <- solve_scalar(
    FUN = function(scale){
      CB_direction_fit(scale, eigenvalues, unit_SRMR, reference, nodes, df)
    },
    lower = 0, upper = maximum_scale,
    metric = metric, target = target, tolerance = 1e-14
  )

  # Exact fit at the solution
  if(isTRUE(scaled$valid)){
    scaled <- c(
      CB_error(direction * scaled$parameter, reference, free_cells, nodes, df),
      list(parameter = scaled$parameter)
    )
  }

  # Check whether the target was reached
  if(isTRUE(scaled$valid) && abs(scaled[[metric]] - target) <= tolerance){
    return(scaled)
  }

  # Otherwise, optimize error directly (errors are bounded by [-1, 1])
  # starting from the closest valid scaled direction (zero is flat)
  start_scale <- swiftelse(
    isTRUE(scaled$valid), max(scaled$parameter, 0.01 * maximum_scale), 0.01 * maximum_scale
  )

  optimization <- try(
    optim(
      par = direction * start_scale,
      fn = CB_error_cost, gr = CB_error_gradient,
      metric = metric, target = target, reference = reference,
      free_cells = free_cells, nodes = nodes, df = df,
      method = "L-BFGS-B", lower = -1, upper = 1
    ), silent = TRUE
  )

  # Check for errors
  if(is(optimization, "try-error")){
    return(swiftelse(isTRUE(scaled$valid), list(scaled), list(list(valid = FALSE)))[[1]])
  }

  # Obtain fit
  fit <- CB_error(optimization$par, reference, free_cells, nodes, df)

  # Return the solution closest to the target
  if(
    fit$valid && (
      !isTRUE(scaled$valid) ||
      abs(fit[[metric]] - target) < abs(scaled[[metric]] - target)
    )
  ){
    return(fit)
  }else if(isTRUE(scaled$valid)){
    return(scaled)
  }else{
    return(list(valid = FALSE))
  }

}

# Tucker, Koopman, and Linn (1969) ----

#' @noRd
# Tucker, Koopman, and Linn (1969) population correlation matrix
# The true network is rescaled and minor factors are added:
#   S R S + T M T, where S^2 = 1 - share, T^2 = share, and M has a unit diagonal
# Rescaling nodes keeps the network in its model class (the zero pattern of the
# precision matrix and all partial correlations are scale invariant) so the
# minor factors are the only error that is added. Subtracting a diagonal from
# the covariance matrix instead would add conditional dependencies among
# non-neighbors. The diagonal is exactly one (each node's variance is
# `1 - share` from the network and `share` from the minor factors)
# Updated 25.11.2026
TKL_population <- function(variance, minor, R, residual_variances)
{

  # Share of each node's variance that is from the minor factors
  share <- variance * residual_variances

  # Return population correlation matrix
  return(
    tcrossprod(sqrt(1 - share)) * R + tcrossprod(sqrt(share)) * minor
  )

}

#' @noRd
# Tucker, Koopman, and Linn (1969) optimization
# Random minor factors with geometrically decaying loadings contribute a
# proportion `variance` of each node's residual variance, 1 / [R^-1]_ii =
# 1 - R^2_i (a network has no common and unique variance so the variance that
# is not explained by the other nodes is the analog of the uniqueness); the
# network part of each node is rescaled to the remaining variance. This is the
# residual network form of generalized network psychometrics: latent factors
# (here, the minor factors) plus a network. The proportion `variance` (between
# 0 and 1) is solved for the target metric
# Updated 25.11.2026
TKL_optimize <- function(
    metric, target, R, network, nodes, df,
    minor_factors, decay, reference = reference_model(R)
)
{

  # Random loadings on minor factors
  minor_factors <- round(minor_factors)
  loadings <- matrix(
    rnorm_ziggurat(nodes * minor_factors),
    nrow = nodes, ncol = minor_factors
  )

  # Geometric decay of loadings across minor factors
  loadings <- loadings * rep((1 - decay)^(seq_len(minor_factors) - 1), each = nodes)

  # Scale each node's loadings to unit length
  loadings <- loadings / sqrt(rowSums(loadings^2))

  # Correlations from minor factors (unit diagonal)
  minor <- tcrossprod(loadings)

  # Residual variances (variance not explained by the other nodes)
  residual_variances <- 1 / diag(reference$inverse)

  # Return solution (the population correlation matrix is positive definite
  # because it is the sum of a positive definite and positive semi-definite matrix)
  return(
    solve_refit(
      function(variance){TKL_population(variance, minor, R, residual_variances)},
      network, df, lower = 0, upper = 1,
      metric = metric, target = target
    )
  )

}

# Wu and Browne (2015) ----

#' @noRd
# Wu and Browne (2015) population correlation matrix
# Updated 20.11.2026
WB_population <- function(log_dispersion, draws, R)
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
  factor_inverse <- try(
    forwardsolve(draws$cholesky %*% bartlett, diag(nodes)), silent = TRUE
  )
  if(is(factor_inverse, "try-error") || any(!is.finite(factor_inverse))){
    return(NULL)
  }
  R_draw <- cov2cor(crossprod(factor_inverse))

  # Check that the errors are within [-1, 1]
  # (the draw is kept as is so that it remains a pure inverted Wishart draw;
  # draws outside the bounds are treated as invalid)
  if(max(abs(R_draw - R)) > 1){
    return(NULL)
  }

  # Return population correlation matrix
  return(R_draw)

}

#' @noRd
# Wu and Browne (2015) optimization
# Updated 25.11.2026
WB_optimize <- function(metric, target, R, network, nodes, df, reference = reference_model(R))
{

  # Random draws (held fixed while the dispersion is optimized)
  normals <- matrix(0, nrow = nodes, ncol = nodes)
  normals[lower.tri(normals)] <- rnorm_ziggurat(nodes * (nodes - 1) / 2)
  draws <- list(
    normals = normals, uniforms = runif_xoshiro(nodes),
    cholesky = t(chol(reference$inverse))
  )

  # Return solution
  # (dispersion v = 1 / m, with m > nodes - 1; v is approximately RMSEA^2)
  return(
    solve_refit(
      function(log_dispersion){WB_population(log_dispersion, draws, R)},
      network, df, lower = log(1e-06), upper = log(1 / nodes),
      metric = metric, target = target, spacing = "linear"
    )
  )

}

# Hu and Bentler (1998) ----

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
# Hu and Bentler (1998) optimization
# (underparameterized misspecification: omitted edges). Omitted edges have
# random signs and partial correlations that are uniformly distributed between
# 0.5 and 1 times a common size that is solved for the target metric
# Updated 25.11.2026
HB_optimize <- function(
    metric, target, R, network, nodes, df,
    free_cells, omitted, reference = reference_model(R),
    max_condition = Inf
)
{

  # Randomly select omitted edges
  omitted_cells <- shuffle(free_cells, size = omitted)

  # Random relative magnitudes (0.5 to 1) with random signs
  # (uniform so that no omitted edge is negligible and none dominates)
  omitted_edges <- matrix(0, nrow = nodes, ncol = nodes)
  omitted_edges[omitted_cells] <-
    runif_xoshiro(length(omitted_cells), min = 0.5, max = 1) *
    shuffle_replace(c(-1, 1), length(omitted_cells))
  omitted_edges <- omitted_edges + t(omitted_edges)

  # Precision matrix of the original population and the shift in
  # the precision matrix that adds `size` to the partial correlations
  precision <- reference$inverse
  omissions <- list(
    precision = precision,
    shift = omitted_edges * sqrt(tcrossprod(diag(precision)))
  )

  # Population correlation matrix (positive definite and errors within [-1, 1]);
  # when `condition` is `TRUE`, the condition number must also be at most
  # `max_condition` (GLASSO becomes very slow on nearly singular matrices
  # because its inner coordinate descent is not bounded by `maxit`)
  population <- function(size, condition = TRUE){

    # Obtain population correlation matrix
    R_population <- HB_population(size, omissions)

    # Check bounds
    if(is.null(R_population) || max(abs(R_population - R)) > 1){
      return(NULL)
    }

    # Check condition number
    if(condition && !within_condition(R_population, max_condition)){
      return(NULL)
    }

    # Return population correlation matrix
    return(R_population)

  }

  # Find largest size (partial correlation of omitted edges) that is allowed
  # (with and without the condition number limit)
  maximum_size <- largest_allowed(function(size){population(size, TRUE)})
  capped <- is.finite(max_condition) &&
    maximum_size < largest_allowed(function(size){population(size, FALSE)})

  # Solve
  solution <- solve_refit(
    population, network, df, lower = 0, upper = maximum_size,
    metric = metric, target = target
  )

  # Flag whether the condition number limit shortened the search and add the
  # omitted edges (partial correlations in the population)
  if(isTRUE(solution$valid)){
    solution$capped <- capped
    solution$omitted <- solution$parameter * omitted_edges
  }

  # Return solution
  return(solution)

}
