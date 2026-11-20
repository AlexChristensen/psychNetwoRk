#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared Simulation Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Bounds for Weibull distribution
# Updated 20.11.2026
check_bounds <- function(shape, scale)
{
  # The bounds were derived from the 222 empirical datasets
  # from Huth et al. (2025), representing most "standard"
  # conditions that occur in empirical data (round up and down)
  # These ranges will be updated as more empirical data
  # are aggregated
  any((shape < 0.70) | (shape > 1.70) | (scale < 0.03) | (scale > 0.19))
}

#' @noRd
# Bootstrapped Weibull parameters within the empirical bounds
# The predicted parameters plus each empirical residual pair; only pairs
# within the empirical bounds are kept (sampling from these is the same as
# rejection sampling of bootstrapped parameters, without the risk of an
# infinite loop when no residual brings the parameters within bounds)
# Updated 20.11.2026
feasible_weibull <- function(nodes, snr)
{

  # Get Weibull model
  weibull_weights <- get(data("weibull_weights", package = "psychNetwoRk", envir = environment()))

  # Predicted parameters
  predicted <- weibull_parameters(nodes = nodes, snr = snr, bootstrap = FALSE)

  # Bootstrapped parameters
  shape <- predicted[["shape"]] + weibull_weights$shape$residuals
  scale <- predicted[["scale"]] + weibull_weights$scale$residuals

  # Keep parameters within bounds
  keep <- !((shape < 0.70) | (shape > 1.70) | (scale < 0.03) | (scale > 0.19))

  # Return feasible parameters
  return(cbind(shape = unname(shape[keep]), scale = unname(scale[keep])))

}

#' @noRd
# Generate from Weibull distribution
# Updated 08.03.2026
weibull_xoshiro <- function(n, shape, scale)
{
  return(scale * (-log(runif_xoshiro(n)))^(1 / shape))
}

#' @noRd
# Maximum likelihood estimates of Weibull parameters
# (profile likelihood over the shape: scale = mean(x^shape)^(1 / shape))
# Updated 20.11.2026
weibull_mle <- function(x)
{

  # Positive values
  x <- x[x > 0]
  n <- length(x)
  log_x <- log(x)

  # Profile log-likelihood of the shape
  profile <- function(shape){
    n * log(shape) - n * log(mean(x^shape)) + (shape - 1) * sum(log_x) - n
  }

  # Maximize over the shape
  shape <- optimize(profile, interval = c(0.05, 20), maximum = TRUE, tol = 1e-08)$maximum

  # Return parameters
  return(c(shape = shape, scale = mean(x^shape)^(1 / shape)))

}

#' @noRd
# Generate edge values
# Updated 20.11.2026
generate_edges <- function(nonzero, p, snr, max_attempts = 1000)
{

  # Calculate possible edges
  possible_edges <- p * (p - 1) / 2

  # Weibull parameters within the empirical bounds
  feasible <- feasible_weibull(p, snr)
  if(nrow(feasible) == 0){
    stop("Weibull parameters are outside of the empirical bounds", call. = FALSE)
  }

  # Generate until no edges are greater than the maximum weight
  for(attempt in seq_len(max_attempts)){

    # Draw Weibull parameters
    params <- feasible[shuffle(seq_len(nrow(feasible)), 1),]

    # Generate edge weights
    weights <- weibull_xoshiro(possible_edges, shape = params[["shape"]], scale = params[["scale"]])

    # Check for weights greater than the maximum weight
    # (no weights greater than max weight = 0.8452469)
    if(!any(weights >= 0.85)){
      break
    }

    # Check for maximum attempts
    if(attempt == max_attempts){
      stop("Could not generate edge weights below the maximum weight", call. = FALSE)
    }

  }

  # Create weighted sampling of weights
  weighted_weights <- weights^2

  # Set weights
  weights <- sort(
    weighted_shuffle(
      x = weights, size = nonzero, prob = weighted_weights / sum(weighted_weights)
    ), decreasing = TRUE
  )

  # Attach attributes to weights
  attr(weights, "params") <- params

  # Return weights
  return(weights)

}

#' @noRd
# Minimally condition the network
# Updated 11.03.2026
condition_network <- function(network, target_condition)
{

  # Obtain pre-inversion matrix
  Omega <- -network
  diag(Omega) <- 1

  # Compute one plus minimum eigenvalue
  min_eigen <- 1 + abs(min(matrix_eigenvalues(Omega)))

  # Condition network
  opt <- uniroot(
    f = function(lambda){

      # Create precision
      diag(Omega) <- min_eigen + lambda

      # Compute the precision matrix
      K <- try(solve(Omega), silent = TRUE)

      # Catch bad matrices
      if(is(K, "try-error")){
        return(target_condition)
      }else{
        return(fast_kappa(cov2cor(K)) - target_condition)
      }

    }, interval = c(1e-05, 0.30) # max ~23% shrinkage (Peeters et al., 2020)
  )

  # Set diagonal on precision
  diag(Omega) <- min_eigen + opt$root

  # Return results
  return(list(R = cov2cor(solve(Omega)), lambda = opt$root))

}

#' @noRd
# Population correlations, conditioning, and validity checks ----
# Shared by `simulate_sbm` and `simulate_smallworld`: obtains the population
# correlation matrix of the network; when it is not positive definite, the
# network is minimally conditioned (ridge) and the Weibull parameters are
# re-estimated by maximum likelihood from the conditioned edge weights and
# re-checked against the empirical bounds. The condition number (with a
# slack of 1 over `target_condition`) and the maximum correlation are checked
# Updated 20.11.2026
condition_and_validate <- function(
    network, params, target_condition, max_correlation, lower_triangle
)
{

  # Get population values
  R <- silent_call(pcor2cor(network))

  # Set lambda
  lambda <- NA

  # Eigenvalues (positive definiteness and condition number)
  eigenvalues <- swiftelse(anyNA(R), NA, matrix_eigenvalues(R))

  # Check for positive definite
  if(anyNA(eigenvalues) || any(eigenvalues <= .Machine$double.eps)){

    # Try to condition network
    output <- try(condition_network(network, target_condition), silent = TRUE)

    # Check for whether correlation matrix is conditioned
    if(is(output, "try-error")){
      stop("Ill-conditioned", call. = FALSE)
    }

    # Store outputs
    lambda <- output$lambda

    # Update network and correlations
    P <- cor2pcor(output$R)
    P[network == 0] <- 0
    dimnames(P) <- dimnames(network)
    network <- P
    R <- silent_call(pcor2cor(network))

    # Check for positive definite again
    eigenvalues <- swiftelse(anyNA(R), NA, matrix_eigenvalues(R))
    if(anyNA(eigenvalues) || any(eigenvalues <= .Machine$double.eps)){
      stop("Not positive definite", call. = FALSE)
    }

    # Re-estimate Weibull parameters from the conditioned edge weights
    params <- weibull_mle(abs(network[lower_triangle][network[lower_triangle] != 0]))
    if(check_bounds(params[["shape"]], params[["scale"]])){
      stop("Weibull parameters outside of empirical bounds after conditioning", call. = FALSE)
    }

  }

  # Check for condition
  condition <- max(abs(eigenvalues)) / min(abs(eigenvalues))
  if(condition > (target_condition + 1)){
    stop("Condition number exceeds target", call. = FALSE)
  }

  # Check for maximum correlation
  if(max(abs(R[lower_triangle])) > max_correlation){
    stop("Max correlation exceeds target", call. = FALSE)
  }

  # Add names
  dimnames(R) <- dimnames(network)

  # Return results
  return(
    list(
      R = R, network = network, params = params,
      lambda = lambda, condition = condition
    )
  )

}

#' @noRd
# Repeated attempts with rejections ----
# `attempt()` returns a result or signals an error whose message is the
# rejection reason; stops after `max_iterations` with a table of reasons
# Updated 20.11.2026
simulation_loop <- function(attempt, max_iterations, structure)
{

  # Collect rejections
  rejections <- character(length = max_iterations)

  # Loop through attempts
  for(iteration in seq_len(max_iterations)){

    # Attempt
    output <- tryCatch(
      attempt(),
      error = function(e){structure(conditionMessage(e), class = "rejection")}
    )

    # Check for rejection
    if(is(output, "rejection")){
      rejections[iteration] <- unclass(output)
      next
    }

    # Return result
    return(list(output = output, iterations = iteration, rejections = rejections))

  }

  # Collect rejections
  rejection_table <- fast_table(rejections)
  rejection_table <- rejection_table[names(rejection_table) != ""]

  # Return error
  stop(
    paste0(
      "Reached maximum iterations: Could not find weights that converge for the ",
      structure, " structure. Resulting rejections were due to:\n\n",
      paste0(names(rejection_table), " = ", rejection_table, collapse = "\n")
    ), call. = FALSE
  )

}

#' @noRd
# Default proportion of negative edges ----
# Empirical distribution of true negative partial correlations across 194
# psychometric networks (mean = 0.34, SD = 0.086), clamped to the empirical
# minimum (0.083) and `maximum`
# Updated 20.11.2026
default_negative_proportion <- function(maximum)
{
  return(min(max(0.34 + rnorm_ziggurat(1) * 0.086, 0.083), maximum))
}

#' @noRd
# Input checks shared by `simulate_sbm` and `simulate_smallworld` ----
# Updated 20.11.2026
check_simulation_inputs <- function(
    nodes, snr, negative_proportion, maximum_negative, sample_size,
    skew, skew_range, target_condition, max_correlation, max_iterations,
    function_name
)
{

  # Error for 'snr'
  typeof_error(snr, "numeric", function_name)
  length_error(snr, 1, function_name)
  range_error(snr, c(0.01, Inf), function_name)

  # Send warning if beyond bounds
  if((snr < 0.64) | (snr > 1.72)){
    warning(
      paste0(
        "The signal-to-noise ratio specified (", round(snr, 2), ") is beyond ",
        "the bounds of cataloged empirical data (0.64-1.72). "
      ), call. = FALSE
    )
  }

  # Check that Weibull parameters within the empirical bounds can be generated
  if(nrow(feasible_weibull(nodes, snr)) == 0){
    stop(
      paste0(
        "The edge weight distribution predicted for ", nodes, " nodes and a ",
        "signal-to-noise ratio of ", round(snr, 2), " is outside of the bounds of ",
        "the empirical data (Weibull shape between 0.70 and 1.70 and scale between ",
        "0.03 and 0.19). Use a signal-to-noise ratio closer to the empirical range ",
        "(0.64-1.72)"
      ), call. = FALSE
    )
  }

  # Errors for 'negative_proportion'
  typeof_error(negative_proportion, "numeric", function_name)
  length_error(negative_proportion, 1, function_name)
  range_error(negative_proportion, c(0, maximum_negative), function_name)

  # Errors for 'sample_size'
  typeof_error(sample_size, "numeric", function_name)
  length_error(sample_size, 1, function_name)
  range_error(sample_size, c(1, Inf), function_name)

  # Errors for 'skew'
  typeof_error(skew, "numeric", function_name)
  length_error(skew, c(1, nodes), function_name)
  range_error(skew, c(-2, 2), function_name)

  # Errors for 'skew_range'
  if(!is.null(skew_range)){
    typeof_error(skew_range, "numeric", function_name)
    length_error(skew_range, 2, function_name)
    range_error(skew_range, c(-2, 2), function_name)
  }

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric", function_name)
  length_error(target_condition, 1, function_name)
  range_error(target_condition, c(1, 100), function_name)

  # Errors for 'max_correlation'
  typeof_error(max_correlation, "numeric", function_name)
  length_error(max_correlation, 1, function_name)
  range_error(max_correlation, c(0, 1), function_name)

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric", function_name)
  length_error(max_iterations, 1, function_name)
  range_error(max_iterations, c(1, Inf), function_name)

}

#' @noRd
# Check for whole numbers ----
# Updated 20.11.2026
whole_number_error <- function(input, function_name)
{
  if(any(input != round(input))){
    stop(
      paste0(
        "Input into '", deparse(substitute(input)), "' must be whole numbers",
        " (in `", function_name, "`)"
      ), call. = FALSE
    )
  }
}

#' @noRd
# Simulate data
# Skew values are rounded to the nearest 0.05; with `skew_range`, one skew
# per variable is drawn from the values in the range (in increments of 0.05)
# Updated 20.11.2026
simulate_data <- function(n, p, R, skew, skew_range)
{

  # Set skew
  if(!is.null(skew_range)){

    # Possible skews in the range (increments of 0.05)
    possible_skews <- round(
      seq(round(min(skew_range) * 20) / 20, round(max(skew_range) * 20) / 20, by = 0.05), 2
    )

    # Draw one skew per variable
    skew <- shuffle_replace(possible_skews, p)

  }else{

    # Round to the nearest 0.05
    skew <- round(round(skew * 20) / 20, 2)

    # One skew per variable
    if(length(skew) == 1){
      skew <- rep(skew, p)
    }

  }

  # Generate data
  data <- matrix(rnorm_ziggurat(n * p), nrow = n, ncol = p) %*% chol(R)

  # Loop through columns
  for(i in seq_len(p)){
    data[,i] <- skew_continuous(skew = skew[i], data = data[,i])
  }

  # Add column names to data
  names(skew) <- colnames(data) <- node_names(p)

  # Return data and skew
  return(list(data = data, skew = skew))

}
