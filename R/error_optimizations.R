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
# and the log-determinant of the population matrix (independence baseline)
# (SRMR is from the residuals or, if supplied, `SRMR` directly)
# Updated 20.11.2026
fit_from_discrepancy <- function(
    discrepancy, log_det_population, residual = NULL, nodes, df,
    n = population_size, R_population = NULL, SRMR = NULL
)
{

  # Chi-square
  chi_square <- n * discrepancy

  # Baseline chi-square (independence model): -log|Sigma_pop|
  baseline_chi_square <- -n * log_det_population
  baseline_df <- nodes * (nodes - 1) / 2

  # Guard against non-finite results
  if(!is.finite(chi_square) || !is.finite(baseline_chi_square)){
    return(list(valid = FALSE))
  }

  # Check for degenerate baseline
  denominator <- max(baseline_chi_square - baseline_df, chi_square - df, 0)
  if(denominator == 0 || baseline_chi_square / baseline_df <= 1){
    return(list(valid = FALSE))
  }

  # CFI, TLI, and RMSEA (shared with `network_fit`)
  indices <- sem_indices(chi_square, df, baseline_chi_square, baseline_df, n)

  # Standardized root mean square residual (lower triangle including the
  # diagonal, p * (p + 1) / 2 unique elements; same as {lavaan} for correlations)
  if(is.null(SRMR)){
    SRMR <- sqrt(mean(residual[lower.tri(residual, diag = TRUE)]^2))
  }

  # Return fit
  return(
    list(
      valid = TRUE, R_population = R_population,
      discrepancy = discrepancy, chi_square = chi_square,
      baseline_chi_square = baseline_chi_square,
      baseline_df = baseline_df, df = df,
      CFI = indices[["CFI"]], TLI = indices[["TLI"]],
      RMSEA = indices[["RMSEA"]], SRMR = SRMR
    )
  )

}

#' @noRd
# Fit indices of a model-implied matrix to a population matrix
# Updated 20.11.2026
fit_indices <- function(
    R_population, reference, df, n = population_size,
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
      nodes = nodes, df = df, n = n, R_population = R_population
    )
  )

}

#' @noRd
# Fit indices of the known graph (constrained maximum likelihood)
# Updated 20.11.2026
refit_indices <- function(R_population, network, df, n = population_size, start = NULL)
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
  fit <- fit_indices(R_population, reference, df, n, population_chol = population_chol)

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
# Updated 20.11.2026
scalar_refit <- function(population, network, df, n = population_size)
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
      return(refit_indices(R_population, network, df, n))

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
# Updated 20.11.2026
solve_refit <- function(
    population, network, df, lower, upper, metric, target,
    n = population_size, spacing = "geometric"
)
{
  return(
    solve_scalar(
      scalar_refit(population, network, df, n), lower, upper,
      metric, target, spacing = spacing
    )
  )
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
# Updated 20.11.2026
CB_error <- function(error_vector, reference, free_cells, nodes, df, n = population_size)
{
  return(
    fit_indices(CB_population(error_vector, reference, free_cells, nodes), reference, df, n)
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
# Updated 20.11.2026
CB_direction_fit <- function(scale, eigenvalues, unit_SRMR, reference, nodes, df, n = population_size)
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
      nodes = nodes, df = df, n = n, SRMR = scale * unit_SRMR
    )
  )

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
# Updated 20.11.2026
CB_error_gradient <- function(error_vector, metric, target, reference, free_cells, nodes, df, n = population_size)
{

  # Obtain fit
  fit <- CB_error(error_vector, reference, free_cells, nodes, df, n)

  # Cost is flat when invalid
  if(!fit$valid){
    return(numeric(length(error_vector)))
  }

  # Inverse of the population matrix
  R_population_inverse <- chol2inv(chol(fit$R_population))

  # Parameter e enters Sigma_pop at (i,j) and (j,i), so
  #   d chi-square / de = 2 * n * (R^-1 - Sigma_pop^-1)_ij
  #   d baseline / de = -2 * n * (Sigma_pop^-1)_ij
  chi_gradient <- 2 * n * (reference$inverse - R_population_inverse)
  baseline_gradient <- -2 * n * R_population_inverse

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
          (fit$SRMR * nodes * (nodes + 1) / 2)
      }

    }
  )

  # d cost / de = 2 * (metric - target) * d metric / de
  return(2 * (fit[[metric]] - target) * metric_gradient[free_cells])

}

#' @noRd
# Cudeck and Browne (1992) optimization
# Updated 20.11.2026
CB_optimize <- function(
    metric, target, reference, free_cells, nodes, df,
    n = population_size, tolerance = 0.01, inverse_chol = backsolve(chol(reference$R), diag(nodes))
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
      CB_direction_fit(scale, eigenvalues, unit_SRMR, reference, nodes, df, n)
    },
    lower = 0, upper = maximum_scale,
    metric = metric, target = target, tolerance = 1e-14
  )

  # Exact fit at the solution
  if(isTRUE(scaled$valid)){
    scaled <- c(
      CB_error(direction * scaled$parameter, reference, free_cells, nodes, df, n),
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
      free_cells = free_cells, nodes = nodes, df = df, n = n,
      method = "L-BFGS-B", lower = -1, upper = 1
    ), silent = TRUE
  )

  # Check for errors
  if(is(optimization, "try-error")){
    return(swiftelse(isTRUE(scaled$valid), list(scaled), list(list(valid = FALSE)))[[1]])
  }

  # Obtain fit
  fit <- CB_error(optimization$par, reference, free_cells, nodes, df, n)

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
# Tucker, Koopman, and Linn (1969) optimization
# Updated 20.11.2026
TKL_optimize <- function(
    metric, target, R, network, nodes, df,
    minor_factors, decay, n = population_size
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
  # (minor factors contribute `weight` of the variance; diagonal remains one)
  return(
    solve_refit(
      function(weight){(1 - weight) * R + weight * minor},
      network, df, lower = 0, upper = maximum_weight,
      metric = metric, target = target, n = n
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
# Updated 20.11.2026
WB_optimize <- function(metric, target, R, network, nodes, df, n = population_size, reference = reference_model(R))
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
      metric = metric, target = target, n = n, spacing = "linear"
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
# Updated 20.11.2026
HB_optimize <- function(
    metric, target, R, network, nodes, df,
    free_cells, omitted, n = population_size, reference = reference_model(R),
    max_condition = Inf
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
    if(condition && is.finite(max_condition)){
      eigenvalues <- matrix_eigenvalues(R_population)
      if(max(eigenvalues) / min(eigenvalues) > max_condition){
        return(NULL)
      }
    }

    # Return population correlation matrix
    return(R_population)

  }

  # Largest size (partial correlation of omitted edges) that is allowed
  # (the smallest size is always allowed)
  largest_size <- function(condition){

    # Check upper bound
    if(!is.null(population(1, condition))){
      return(1)
    }

    # Bisection
    lower <- 0; upper <- 1
    for(i in seq_len(50)){

      # Update bounds
      middle <- (lower + upper) / 2
      if(!is.null(population(middle, condition))){
        lower <- middle
      }else{
        upper <- middle
      }

    }

    # Stay away from the boundary
    return(lower * 0.999)

  }

  # Find largest size (with and without the condition number limit)
  maximum_size <- largest_size(TRUE)
  capped <- is.finite(max_condition) && maximum_size < largest_size(FALSE)

  # Solve
  solution <- solve_refit(
    population, network, df, lower = 0, upper = maximum_size,
    metric = metric, target = target, n = n
  )

  # Flag whether the condition number limit shortened the search
  if(isTRUE(solution$valid)){
    solution$capped <- capped
  }

  # Return solution
  return(solution)

}
