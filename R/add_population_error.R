#' Adds Population Error to \code{\link[psychNetwoRk]{simulate_sbm}} and
#' \code{\link[psychNetwoRk]{simulate_smallworld}} Data
#'
#' Adds population (model) error to simulated Gaussian graphical model data.
#' The population correlation matrix is perturbed so that the true network
#' (known graph) no longer fits perfectly but instead fits at a specified
#' value of a fit index (CFI, TLI, RMSEA, or SRMR). Data are then
#' re-generated from the perturbed population correlation matrix.
#' See examples to get started
#'
#' @details
#' The error is generated with one of four methods (see \code{method}) and
#' optimized so that the \code{metric} of the true network, fit to the
#' population correlation matrix with error using constrained maximum
#' likelihood (i.e., the known graph), equals \code{fit_value}. The fit indices
#' are computed at the population level using the maximum likelihood
#' discrepancy function between the population correlation matrix with error
#' and the correlation matrix implied by the true network. The largest
#' absolute error added to any cell of the original population correlation
#' matrix is limited by \code{max_residual}.
#'
#' Each optimization attempt uses new random starting values or random draws
#' of the error (up to \code{max_iterations}) and the attempt closest to
#' \code{fit_value} is returned. If \code{fit_value} cannot be reached within
#' \code{tolerance} (e.g., because \code{max_residual} is too small), a warning
#' is given and the closest solution is used
#'
#' @param ggm_object Data object from \code{\link[psychNetwoRk]{simulate_sbm}} or
#' \code{\link[psychNetwoRk]{simulate_smallworld}}
#'
#' @param method Character (length = 1).
#' Method used to generate the population error.
#' Defaults to \code{"cb"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"cb"} --- Cudeck and Browne (1992). Error is placed only on
#' the cells that the true network constrains to zero, so the true network
#' remains the best-fitting model of the population correlation matrix with
#' error (the error is orthogonal to the true network)
#'
#' \item \code{"tkl"} --- Tucker, Koopman, and Linn (1969). Error is
#' the contribution of many random minor factors whose loadings decay
#' geometrically. The proportion of variance attributed to the minor factors
#' is optimized
#'
#' \item \code{"wb"} --- Wu and Browne (2015). The population correlation
#' matrix with error is a random draw from an inverted Wishart distribution
#' centered on the original population correlation matrix. The dispersion
#' of the distribution is optimized
#'
#' \item \code{"hb"} --- Hu and Bentler (1998). Underparameterized
#' misspecification: a random subset of the edges that are zero in the true
#' network are non-zero in the population (omitted edges). The size of the
#' omitted edges is optimized
#'
#' }
#'
#' @param metric Character (length = 1).
#' Fit index used as the objective criterion. The true network is fit to the
#' population correlation matrix with error (constrained maximum likelihood)
#' and the error is optimized so that \code{metric} equals \code{fit_value}.
#' Defaults to \code{"CFI"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"CFI"} --- Comparative fit index
#'
#' \item \code{"TLI"} --- Tucker-Lewis index
#'
#' \item \code{"RMSEA"} --- Root mean square error of approximation
#'
#' \item \code{"SRMR"} --- Standardized root mean square residual
#'
#' }
#'
#' @param fit_value Numeric (length = 1).
#' Target value of \code{metric} for the true network fit to the population
#' correlation matrix with error. Must be between 0 and 1 (exclusive).
#' Defaults to the conventional cutoff for the \code{metric}:
#' \code{0.95} for \code{"CFI"} and \code{"TLI"}, \code{0.06} for
#' \code{"RMSEA"}, and \code{0.08} for \code{"SRMR"}
#'
#' @param max_residual Numeric (length = 1).
#' Maximum absolute error that can be added to any cell of the original
#' population correlation matrix. If the target \code{fit_value} cannot be reached
#' without exceeding this value, the closest solution is used (with a warning).
#' Defaults to \code{0.10}
#'
#' @param tolerance Numeric (length = 1).
#' Maximum absolute difference between the target \code{fit_value} and the
#' \code{metric} of the population correlation matrix with error. Ensures that
#' the appropriate amount of population error was added.
#' Defaults to \code{0.001}
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of optimization attempts (random starting values or random
#' draws of the error) to reach the target \code{fit_value} within the specified
#' \code{tolerance}.
#' Defaults to \code{10}
#'
#' @param control List.
#' Method-specific options. Options that do not apply to \code{method} are
#' ignored. Available options:
#'
#' \itemize{
#'
#' \item \code{minor_factors} --- Number of minor factors for
#' \code{method = "tkl"}. Defaults to \code{max(50, 2 * nodes)}
#'
#' \item \code{decay} --- Geometric decay (\eqn{\epsilon}) of the minor factor
#' loadings for \code{method = "tkl"}. Must be between 0 and 1.
#' Defaults to \code{0.02}
#'
#' \item \code{omitted} --- Proportion of the zero edges in the true network
#' that are omitted edges for \code{method = "hb"}. At least one edge is always
#' omitted. Defaults to \code{0.10}
#'
#' }
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
#' \item \code{method} --- Method used to generate the population error
#'
#' \item \code{metric} --- Fit index used as the objective criterion
#'
#' \item \code{fit_value} --- Target value of \code{metric}
#'
#' \item \code{achieved} --- Value of \code{metric} achieved by the true network
#' fit to the population correlation matrix with error
#'
#' \item \code{fit} --- Named vector of all fit indices (\code{CFI},
#' \code{TLI}, \code{RMSEA}, and \code{SRMR}) achieved
#'
#' \item \code{delta} --- Minimum of the (maximum likelihood) discrepancy
#' function corresponding to the achieved fit
#'
#' \item \code{error} --- Matrix of the error that was added to the
#' original population correlation matrix
#'
#' \item \code{max_error} --- Largest absolute error that was added (never
#' larger than \code{max_residual})
#'
#' \item \code{max_residual} --- Maximum absolute error allowed
#'
#' \item \code{R_difference} --- Largest absolute difference between the
#' original population correlation matrix and the known graph (constrained
#' maximum likelihood) fit to the population correlation matrix with error.
#' Values near zero indicate that the error is orthogonal to the true network
#' (expected for \code{method = "cb"} only)
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
#' # Add population error (CFI = 0.95 using Cudeck and Browne)
#' sbm_cfi <- add_population_error(ggm_object = sbm)
#'
#' \donttest{
#' # Use Wu and Browne to target TLI = 0.90 with a larger maximum error
#' sbm_wb <- add_population_error(
#'   ggm_object = sbm, method = "wb",
#'   metric = "TLI", fit_value = 0.90, max_residual = 0.30
#' )
#'
#' # Target RMSEA = 0.04 with a different maximum error
#' sbm_rmsea <- add_population_error(
#'   ggm_object = sbm, metric = "RMSEA",
#'   fit_value = 0.04, max_residual = 0.20
#' )
#'
#' # Generate small-world data and add population error
#' sw <- simulate_smallworld(
#'   nodes = 20, density = 0.20, rewire = 0.10,
#'   sample_size = 1000
#' )
#'
#' sw_srmr <- add_population_error(
#'   ggm_object = sw, metric = "SRMR", fit_value = 0.05
#' )}
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{Cudeck and Browne (\code{"cb"})} \cr
#' Cudeck, R., & Browne, M.W. (1992).
#' Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value.
#' \emph{Psychometrika}, \emph{57}, 357--369.
#'
#' \strong{Tucker, Koopman, and Linn (\code{"tkl"})} \cr
#' Tucker, L.R., Koopman, R.F., & Linn, R.L. (1969).
#' Evaluation of factor analytic research procedures by means of simulated correlation matrices.
#' \emph{Psychometrika}, \emph{34}, 421--459.
#'
#' \strong{Wu and Browne (\code{"wb"})} \cr
#' Wu, H., & Browne, M.W. (2015).
#' Quantifying adventitious error in a covariance structure as a random effect.
#' \emph{Psychometrika}, \emph{80}, 571--600.
#'
#' \strong{Hu and Bentler (\code{"hb"}) and fit indices} \cr
#' Hu, L., & Bentler, P.M. (1998).
#' Fit indices in covariance structure modeling: Sensitivity to underparameterized model misspecification.
#' \emph{Psychological Methods}, \emph{3}, 424--453.
#'
#' @export
#'
# Add population error to simulated data ----
# Updated 04.11.2026
add_population_error <- function(
    ggm_object,
    method = c("cb", "tkl", "wb", "hb"),
    metric = c("CFI", "TLI", "RMSEA", "SRMR"),
    fit_value = NULL, max_residual = 0.10,
    tolerance = 0.001, max_iterations = 10,
    control = list()
)
{

  # Set defaults
  method <- set_default(method, "cb", add_population_error)
  metric <- toupper(set_default(metric, "cfi", add_population_error))

  # Check inputs
  inputs <- add_population_error_errors(
    ggm_object, method, metric, fit_value,
    max_residual, tolerance, max_iterations, control
  )
  target <- inputs$fit_value
  control <- inputs$control

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

  # Original population correlation matrix as the reference
  reference <- reference_model(R)

  # Initialize best solution
  best <- NULL
  best_difference <- Inf

  # Loop through attempts
  for(iteration in seq_len(max_iterations)){

    # Optimize error
    solution <- switch(
      method,
      "cb" = CB_optimize(
        metric, target, reference, free_cells,
        nodes, df, max_residual
      ),
      "tkl" = TKL_optimize(
        metric, target, R, network, nodes, df, max_residual,
        control$minor_factors, control$decay
      ),
      "wb" = WB_optimize(
        metric, target, R, network, nodes, df, max_residual
      ),
      "hb" = HB_optimize(
        metric, target, R, network, nodes, df, max_residual,
        free_cells, control$omitted
      )
    )

    # Check for invalid fit
    if(!solution$valid){
      next
    }

    # Obtain difference from target
    difference <- abs(solution[[metric]] - target)

    # Check for improvement
    if(difference < best_difference){

      # Store results
      best <- list(fit = solution, iterations = iteration)
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

  # Obtain population correlation matrix with error
  R_population <- best$fit$R_population
  dimnames(R_population) <- dimnames(R)

  # Obtain largest absolute error
  max_error <- max(abs(R_population - R))

  # Check for convergence
  if(best_difference > tolerance){

    # Check whether the maximum error was the limiting factor
    limit_text <- swiftelse(
      max_error >= max_residual * (1 - 1e-06),
      paste0(
        " The largest absolute error reached `max_residual` = ",
        max_residual, ". Try increasing `max_residual`."
      ),
      ""
    )

    warning(
      paste0(
        "Convergence counter has exceeded its limit. ",
        "There were issues converging to the target ", metric, " = ", target,
        " within `tolerance` = ", tolerance, ". \n\n",
        "Using the solution with the closest ", metric, " = ",
        round(best$fit[[metric]], 4), ".", limit_text
      ), call. = FALSE
    )

  }

  # Refit the known graph to obtain the achieved fit
  # (also verifies whether the error is orthogonal to the network)
  known <- known_graph(R_population, network)
  R_difference <- max(abs(known$W - R))
  fit <- refit_indices(R_population, network, df)

  # Check for invalid refit
  if(!fit$valid){
    fit <- best$fit
  }

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
      method = method,
      metric = metric,
      fit_value = target,
      achieved = fit[[metric]],
      fit = unlist(fit[c("CFI", "TLI", "RMSEA", "SRMR")]),
      delta = fit$discrepancy,
      error = R_population - R,
      max_error = max_error,
      max_residual = max_residual,
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
# Updated 04.11.2026
add_population_error_errors <- function(
    ggm_object, method, metric, fit_value,
    max_residual, tolerance, max_iterations, control
)
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

  # Check for available methods
  available_methods <- c("cb", "tkl", "wb", "hb")
  if(!method %in% available_methods){
    stop(
      paste0(
        "`method` = \"", method, "\" is not available.",
        "\n\nUse ", paste0("`method = \"", available_methods, "\"`", collapse = ", ")
      )
    )
  }

  # Check for metric
  if(!metric %in% c("CFI", "TLI", "RMSEA", "SRMR")){
    stop(
      paste0(
        "`metric` = \"", metric, "\" is not available.",
        "\n\nUse `metric = \"CFI\"`, `metric = \"TLI\"`, ",
        "`metric = \"RMSEA\"`, or `metric = \"SRMR\"`"
      )
    )
  }

  # Set default fit_value (conventional cutoffs)
  if(is.null(fit_value)){
    fit_value <- switch(
      metric,
      "CFI" = 0.95,
      "TLI" = 0.95,
      "RMSEA" = 0.06,
      "SRMR" = 0.08
    )
  }

  # Check for appropriate fit_value
  typeof_error(fit_value, c("integer", "double"), "add_population_error")
  length_error(fit_value, 1, "add_population_error")

  # Check for range
  if(!is.finite(fit_value) || fit_value <= 0 || fit_value >= 1){
    stop(
      paste0(
        "`fit_value` = ", fit_value, " is not a valid ", metric, ".",
        "\n\nNumeric `fit_value` must be between 0 and 1 (exclusive)"
      )
    )
  }

  # Check maximum residual
  typeof_error(max_residual, c("integer", "double"), "add_population_error")
  length_error(max_residual, 1, "add_population_error")
  if(!is.finite(max_residual) || max_residual <= 0){
    stop(
      paste0(
        "`max_residual` = ", max_residual, " is not valid.",
        "\n\nNumeric `max_residual` must be greater than 0"
      )
    )
  }

  # Check tolerance and maximum iterations
  typeof_error(tolerance, c("integer", "double"), "add_population_error")
  length_error(tolerance, 1, "add_population_error")
  typeof_error(max_iterations, c("integer", "double"), "add_population_error")
  length_error(max_iterations, 1, "add_population_error")
  range_error(max_iterations, c(1, Inf), "add_population_error")

  # Check control
  if(!is.list(control)){
    stop("`control` must be a list. See `?add_population_error` for options")
  }

  # Set control defaults
  control_defaults <- list(
    minor_factors = max(50, 2 * ncol(ggm_object$population$R)),
    decay = 0.02, omitted = 0.10
  )

  # Check for unknown options
  unknown <- setdiff(names(control), names(control_defaults))
  if(length(unknown) > 0 || (length(control) > 0 && is.null(names(control)))){
    stop(
      paste0(
        "`control` option(s) not available: ",
        paste0("\"", unknown, "\"", collapse = ", "),
        "\n\nUse ", paste0("\"", names(control_defaults), "\"", collapse = ", ")
      )
    )
  }

  # Update defaults
  control_defaults[names(control)] <- control
  control <- control_defaults

  # Check control options
  for(option in names(control)){
    typeof_error(control[[option]], c("integer", "double"), "add_population_error")
    length_error(control[[option]], 1, "add_population_error")
  }
  range_error(control$minor_factors, c(1, Inf), "add_population_error")
  range_error(control$decay, c(0, 1), "add_population_error")
  range_error(control$omitted, c(0, 1), "add_population_error")

  # Return inputs
  return(list(fit_value = fit_value, control = control))

}

# Debugging ----
# dm <- matrix(0.20, nrow = 3, ncol = 3); diag(dm) <- 0.90
# ggm_object <- simulate_sbm(
#   nodes = 6, blocks = 3, density_matrix = dm, sample_size = 1000
# )
# method = "cb"; metric = "CFI"; fit_value = NULL; max_residual = 0.10
# tolerance = 0.001; max_iterations = 10; control = list()
