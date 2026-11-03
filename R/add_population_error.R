#' Adds Population Error to \code{\link[psychNetwoRk]{simulate_sbm}} and
#' \code{\link[psychNetwoRk]{simulate_smallworld}} Data
#'
#' Adds population (model) error to simulated Gaussian graphical model data.
#' The population correlation matrix is perturbed so that the true network
#' (known graph) no longer fits perfectly but instead fits at a specified
#' comparative fit index (CFI). Data are then re-generated from the perturbed
#' population correlation matrix. See examples to get started
#'
#' @param ggm_object Data object from \code{\link[psychNetwoRk]{simulate_sbm}} or
#' \code{\link[psychNetwoRk]{simulate_smallworld}}
#'
#' @param misfit Character or numeric (length = 1).
#' Magnitude of error to add expressed as the CFI
#' of the true network fit to the population correlation matrix with error.
#' Defaults to \code{"close"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"close"} --- Slight deviations from the original population
#' correlation matrix (CFI = 0.99)
#'
#' \item \code{"acceptable"} --- Moderate deviations from the original population
#' correlation matrix (CFI = 0.95)
#'
#' }
#'
#' A single number between 0 and 1 (exclusive) can be used to directly
#' specify the target CFI
#'
#' @param tolerance Numeric (length = 1).
#' Maximum absolute difference between the target CFI and the CFI
#' of the population correlation matrix with error. Ensures that
#' the appropriate amount of population error was added.
#' Defaults to \code{0.001}
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of optimization attempts (random starting values)
#' to reach the target CFI within the specified \code{tolerance}.
#' Defaults to \code{10}
#'
#' @return Returns a list of the same class(es) and structure as the
#' \code{ggm_object} input (see \code{\link[psychNetwoRk]{simulate_sbm}} and
#' \code{\link[psychNetwoRk]{simulate_smallworld}}) with the class \code{"ggm_pe"}
#' appended. The elements below mirror the input; additional elements specific to
#' population error follow:
#'
#' \item{data}{Simulated data from the population correlation matrix with error
#' (same dimensions, column names, and skew as the input data)}
#'
#' \item{parameters}{Parameters from the \code{ggm_object} input (unchanged)}
#'
#' \item{population}{
#' A list containing the population parameters. \code{R} is the population
#' correlation matrix with error added; all other elements are carried over from
#' the \code{ggm_object} input (e.g., \code{Omega}, the true network, and
#' \code{membership} for \code{\link[psychNetwoRk]{simulate_sbm}})
#' }
#'
#' \item{convergence}{Convergence information from the \code{ggm_object} input (unchanged)}
#'
#' \item{original_correlation}{Original population correlation matrix \emph{before}
#' population error was added}
#'
#' \item{population_error}{
#' A list containing the parameters used to generate population error:
#'
#' \itemize{
#'
#' \item \code{misfit} --- Specified misfit value
#'
#' \item \code{target_CFI} --- Target CFI
#'
#' \item \code{CFI} --- CFI achieved
#'
#' \item \code{delta} --- Minimum of the (maximum likelihood) discrepancy
#' function corresponding to the misfit value
#'
#' \item \code{error} --- Matrix of the error that was added to the
#' original population correlation matrix
#'
#' \item \code{R_difference} --- Largest absolute difference between the
#' original population correlation matrix and the known graph (constrained
#' maximum likelihood) fit to the population correlation matrix with error.
#' Values near zero indicate that the error is orthogonal to the true network
#'
#' \item \code{iterations} --- Number of optimization attempts used
#'
#' }
#'
#' }
#'
#' \item{original_results}{Original \code{ggm_object} input into function}
#'
#' @examples
#' # Generate SBM data
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- 0.90
#'
#' sbm <- simulate_sbm(
#'   nodes = 6, blocks = 3,
#'   density_matrix = dm,
#'   sample_size = 1000
#' )
#'
#' # Add small population error
#' sbm_close <- add_population_error(
#'   ggm_object = sbm, misfit = "close"
#' )
#'
#' \donttest{
#' # Add moderate population error
#' sbm_acceptable <- add_population_error(
#'   ggm_object = sbm, misfit = "acceptable"
#' )
#'
#' # Generate small-world data and add population error
#' sw <- simulate_smallworld(
#'   nodes = 20, density = 0.20, rewire = 0.10,
#'   sample_size = 1000
#' )
#'
#' sw_close <- add_population_error(ggm_object = sw)}
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#'
#' @references
#' Cudeck, R., & Browne, M.W. (1992).
#' Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value.
#' \emph{Psychometrika}, \emph{57}, 357--369.
#'
#' @export
#'
# Add population error to simulated data ----
# Updated 03.11.2026
add_population_error <- function(
    ggm_object, misfit = c("close", "acceptable"),
    tolerance = 0.01, max_iterations = 10
)
{

  # Check for missing misfit
  if(missing(misfit)){
    misfit <- "close"
  }

  # Check inputs
  target_CFI <- add_population_error_errors(
    ggm_object, misfit, tolerance, max_iterations
  )

  # Obtain population parameters
  R <- ggm_object$population$R
  network <- ggm_object$population$Omega
  nodes <- ncol(R)

  # Identify free cells (non-edges in the upper triangle)
  # Error is placed on cells that the true network constrains to zero
  free_cells <- which(network == 0 & upper.tri(network))

  # Check for free cells
  if(length(free_cells) == 0){
    stop(
      paste(
        "The network in `ggm_object` is fully connected so there are no",
        "constrained (zero) edges to add population error to.",
        "\n\nUse a sparser network to add population error"
      )
    )
  }

  # Degrees of freedom of the known graph
  df <- length(free_cells)

  # Inverse and log-determinant of original population correlation matrix
  R_chol <- chol(R)
  R_inverse <- chol2inv(R_chol)
  log_det_R <- 2 * sum(log(diag(R_chol)))

  # Initialize best solution
  best <- NULL
  best_difference <- Inf

  # Loop through attempts
  for(iteration in seq_len(max_iterations)){

    # Random starting values (zero is a flat region of the cost)
    start <- rnorm_ziggurat(df) * 0.005

    # Optimize error
    optimization <- try(
      optim(
        par = start, fn = error_cost, gr = error_gradient,
        target_CFI = target_CFI, R = R, R_inverse = R_inverse,
        log_det_R = log_det_R, free_cells = free_cells,
        nodes = nodes, df = df, method = "L-BFGS-B",
        lower = -0.25, upper = 0.25
      ), silent = TRUE
    )

    # Check for errors
    if(is(optimization, "try-error")){
      next
    }

    # Obtain fit
    fit <- error_function(
      optimization$par, R, R_inverse, log_det_R,
      free_cells, nodes, df
    )

    # Check for invalid fit
    if(!fit$valid){
      next
    }

    # Obtain difference from target
    difference <- abs(fit$CFI - target_CFI)

    # Check for improvement
    if(difference < best_difference){

      # Store results
      best <- list(
        fit = fit, error_vector = optimization$par,
        iterations = iteration
      )
      best_difference <- difference

    }

    # Check for convergence
    if(best_difference <= tolerance){
      break
    }

  }

  # Check for no solution
  if(is.null(best)){
    stop(
      paste(
        "Population error could not be generated with a positive",
        "definite matrix. \n\n",
        "Try increasing the number of iterations using the",
        "`max_iterations` argument."
      )
    )
  }

  # Check for convergence
  if(best_difference > tolerance){

    warning(
      paste0(
        "Convergence counter has exceeded its limit. ",
        "There were issues converging to the target CFI = ", target_CFI,
        " within `tolerance` = ", tolerance, ". \n\n",
        "Using the solution with the closest CFI = ", round(best$fit$CFI, 4)
      ), call. = FALSE
    )

  }

  # Obtain population correlation matrix with error
  R_population <- best$fit$R_population
  dimnames(R_population) <- dimnames(R)

  # Refit the known graph to verify that the error is orthogonal to the network
  known <- known_graph(R_population, network)
  R_difference <- max(abs(known$W - R))

  # Generate data
  data_output <- simulate_data(
    n = nrow(ggm_object$data), p = nodes, R = R_population,
    skew = ggm_object$parameters$skew, skew_range = NULL
  )

  # Update population
  population <- ggm_object$population
  population$R <- R_population

  # Populate results
  # (same elements and order as `simulate_sbm` and `simulate_smallworld`,
  # followed by the elements specific to population error)
  results <- list(
    data = data_output$data,
    parameters = ggm_object$parameters,
    population = population,
    convergence = ggm_object$convergence,
    original_correlation = R,
    population_error = list(
      misfit = misfit,
      target_CFI = target_CFI,
      CFI = best$fit$CFI,
      delta = best$fit$chi_square / formals(error_function)$n,
      error = R_population - R,
      R_difference = R_difference,
      iterations = best$iterations
    ),
    original_results = ggm_object
  )

  # Add class
  class(results) <- c(class(ggm_object), "ggm_pe")

  # Return results
  return(results)

}

# Input checking ----
#' @noRd
# Updated 03.11.2026
add_population_error_errors <- function(ggm_object, misfit, tolerance, max_iterations)
{

  # Check for appropriate class
  if(!is(ggm_object, "ggm_simulate")){

    # Produce error
    stop(
      paste(
        "`ggm_object` input is not class \"ggm_simulate\" from the",
        "`simulate_sbm` or `simulate_smallworld` functions.",
        "\n\nInput class(es) of current `ggm_object`:",
        paste0("\"", class(ggm_object), "\"", collapse = ", "),
        "\n\nUse `simulate_sbm` or `simulate_smallworld` to generate your",
        "data to input into this function"
      )
    )

  }

  # Check that population error has not yet been added
  if(is(ggm_object, "ggm_pe")){

    # Produce error
    stop(
      paste(
        "`ggm_object` input is class \"ggm_pe\" from the `add_population_error` function.",
        "Population error has already been added.",
        "\n\nUse the `original_results` element of this object or generate new data with",
        "`simulate_sbm` or `simulate_smallworld`"
      )
    )

  }

  # Check for appropriate misfit
  length_error(misfit, 1)

  # Obtain target CFI
  if(is.character(misfit)){

    # Standardize and check options
    misfit <- tolower(misfit)

    # Obtain target CFI
    target_CFI <- switch(
      misfit,
      "close" = 0.99,
      "acceptable" = 0.95,
      stop(
        paste0(
          "`misfit` = \"", misfit, "\" is not available.",
          "\n\nUse `misfit = \"close\"`, `misfit = \"acceptable\"`, or a number between 0 and 1"
        )
      )
    )

  }else{

    # Check for numeric
    typeof_error(misfit, c("integer", "double"))

    # Check for range
    if(!is.finite(misfit) || misfit <= 0 || misfit >= 1){
      stop(
        paste0(
          "`misfit` = ", misfit, " is not a valid CFI.",
          "\n\nNumeric `misfit` must be between 0 and 1 (exclusive)"
        )
      )
    }

    # Set target CFI
    target_CFI <- misfit

  }

  # Check tolerance and maximum iterations
  length_error(tolerance, 1)
  typeof_error(tolerance, c("integer", "double"))
  length_error(max_iterations, 1)
  typeof_error(max_iterations, c("integer", "double"))

  # Return target CFI
  return(target_CFI)

}

# Debugging ----
# dm <- matrix(0.20, nrow = 3, ncol = 3); diag(dm) <- 0.90
# ggm_object <- simulate_sbm(
#   nodes = 6, blocks = 3, density_matrix = dm, sample_size = 1000
# )
# misfit = "close"; tolerance = 0.001; max_iterations = 10
