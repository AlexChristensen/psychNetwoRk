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
#' likelihood (i.e., the known graph), equals \code{fit_value}. The error added
#' to any cell of the original population correlation matrix is limited to the
#' range -1 to 1 and the population correlation matrix with error is always
#' positive definite.
#'
#' The fit indices are population values (infinite sample size) computed from
#' the maximum likelihood discrepancy \eqn{F} between the population
#' correlation matrix with error, \eqn{\Sigma}, and the correlation matrix
#' implied by the true network fit to it (\eqn{df} are the degrees of freedom
#' of the true network, the number of edges that are zero):
#' \eqn{RMSEA = \sqrt{F / df}} and \eqn{CFI = 1 - F / F_0}, where
#' \eqn{F_0 = -\log|\Sigma|} is the discrepancy of the independence model
#' (\eqn{TLI = 1 - (F / df) / (F_0 / df_0)} with
#' \eqn{df_0 = p (p - 1) / 2} for \eqn{p} nodes). These are the values of
#' \code{\link{network_fit}} for the known graph as the sample size approaches
#' infinity
#'
#' Each optimization attempt uses new random draws of the error. Attempts
#' continue (up to \code{max_iterations}) until \code{fit_value} is matched
#' (to within \eqn{10^{-6}}) and the attempt closest to \code{fit_value} is
#' returned. If \code{fit_value} cannot be reached within
#' \code{tolerance} (e.g., because the target requires more error than can be
#' added while remaining positive definite), a warning is given and the closest
#' solution is used
#'
#' @param pn_object Data object from \code{\link[psychNetwoRk]{simulate_sbm}} or
#' \code{\link[psychNetwoRk]{simulate_smallworld}}
#'
#' @param method Character (length = 1).
#' Method used to generate the population error.
#' There is no default; \code{method} must be specified.
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
#' geometrically. The minor factors take a proportion of each node's residual
#' variance (the variance of the node that is not explained by the other
#' nodes, \eqn{1 - R^2}, the analog of a uniqueness in a network) and this
#' proportion is optimized. The network is rescaled to the remaining variance
#' of each node (the partial correlations of the true network are unchanged),
#' so the population is the true network plus latent minor factors (a residual
#' network model)
#'
#' \item \code{"wb"} --- Wu and Browne (2015). The population correlation
#' matrix with error is the inverse of a Wishart draw for the precision matrix
#' (i.e., an inverted Wishart draw for the covariance matrix), where the
#' Wishart is centered on the original (sparse) precision matrix and the result
#' is standardized to a correlation matrix. Alternative scalings of the
#' centering (e.g., centering the covariance rather than the precision matrix)
#' give identical results after standardization. The dispersion
#' of the distribution is optimized. Draws are not rescaled; draws with an
#' error outside of -1 to 1 are treated as infeasible
#'
#' \item \code{"hb"} --- Hu and Bentler (1998). Underparameterized
#' misspecification: a random subset of the edges that are zero in the true
#' network are non-zero in the population (omitted edges). The omitted edges
#' have random signs and partial correlations that vary uniformly between one
#' half and one times a common size, and the common size is optimized
#'
#' }
#'
#' @param metric Character (length = 1).
#' Fit index used as the objective criterion. The true network is fit to the
#' population correlation matrix with error (constrained maximum likelihood)
#' and the error is optimized so that \code{metric} equals \code{fit_value}.
#' Defaults to the fit index that is consistent with \code{method}:
#' \code{"CFI"} for \code{"cb"} and \code{"tkl"}, \code{"RMSEA"} for
#' \code{"wb"} (the dispersion is approximately RMSEA\eqn{^2}), and
#' \code{"SRMR"} for \code{"hb"} (error is defined through the size of
#' omitted edges).
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
#' (computed over the lower triangle of the correlation matrices including
#' the diagonal, following \pkg{lavaan})
#'
#' }
#'
#' @param fit_value Numeric (length = 1).
#' Target value of \code{metric} for the true network fit to the population
#' correlation matrix with error. Must be between 0 and 1 (exclusive).
#' Defaults to the conventional cutoff for the \code{metric}:
#' \code{0.95} for \code{"CFI"} and \code{"TLI"}, \code{0.06} for
#' \code{"RMSEA"}, and \code{0.05} for \code{"SRMR"}
#'
#' @param tolerance Numeric (length = 1).
#' Maximum absolute difference between the target \code{fit_value} and the
#' \code{metric} of the population correlation matrix with error that is
#' accepted without a warning (must be positive). Ensures that the appropriate
#' amount of population error was added.
#' Defaults to \code{0.01}
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of optimization attempts (random starting values or random
#' draws of the error) to reach the target \code{fit_value} within the specified
#' \code{tolerance}.
#' Defaults to \code{10}
#'
#' @param seed Numeric (length = 1).
#' Seed for reproducible results. The same \code{seed} gives the same results
#' (and does not change the random number stream of the R session).
#' Defaults to \code{NULL} (random seed)
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
#' The loading of minor factor \eqn{k} is multiplied by
#' \eqn{(1 - \epsilon)^{k - 1}}, so the variance it contributes is proportional
#' to \eqn{(1 - \epsilon)^{2(k - 1)}}. \code{decay} therefore controls how
#' concentrated the error is across the minor factors:
#'
#' \itemize{
#'
#' \item \emph{Smaller values} (e.g., \code{0.02}) --- Minor factors decay
#' slowly, so many of them contribute about equally. The error approaches
#' diffuse, unstructured noise (the minor-factor correlation matrix
#' approaches an identity matrix): many small error correlations. This
#' resembles very many unrelated, tiny influences and
#' yields little coherent structure in the error
#'
#' \item \emph{Larger values} (e.g., \code{0.30}, the default) --- Minor factors
#' decay quickly, so a handful of factors dominate the error. The error is
#' more structured: groups of nodes share small unmodeled common causes
#' (e.g., method or wording effects, omitted variables), which is more like
#' the misspecification expected in psychometric data. The same amount of
#' misfit is carried by fewer, larger error correlations
#'
#' \item \emph{Very large values} (e.g., \code{> 0.60}) --- Essentially one or two
#' minor factors carry all of the error, which produces a few large
#' error correlations that may be unrealistic as general population
#' error and can make the target \code{fit_value} harder to reach
#' within the error bounds (\eqn{[-1, 1]})
#'
#' }
#'
#' The \code{minor_factors} option matters little once \code{decay} is large
#' because later factors contribute negligible variance.
#' Defaults to \code{0.30}
#'
#' \item \code{omitted} --- Number of edges that are zero in the true network
#' that are omitted edges for \code{method = "hb"} (at least one). If not
#' specified, the number is limited to the number of zero edges. Defaults to
#' \code{2}
#'
#' \item \code{max_condition} --- Largest condition number of the population
#' correlation matrix with error for \code{method = "hb"}. The estimation of
#' the true network (GLASSO) becomes very slow for nearly singular matrices,
#' which occur at the largest errors the omitted edges can add. In
#' simulations across a wide range of networks, a limit of \code{1000}
#' keeps nearly all of the reachable misfit (CFI, TLI, and SRMR reach their
#' worst fit well before the limit; RMSEA loses about 2\% on average) and
#' avoids fits that take many seconds. If \code{fit_value} requires more error
#' than the limit allows, a warning is given and the closest solution is used.
#' Use \code{Inf} for no limit. Defaults to \code{1000}
#'
#' }
#'
#' @return Returns a list of the same class(es) and structure as the
#' \code{pn_object} input (see \code{\link[psychNetwoRk]{simulate_sbm}} and
#' \code{\link[psychNetwoRk]{simulate_smallworld}}) with the class \code{"pn_pe"}
#' appended. The elements below mirror the input; additional elements specific to
#' population error follow:
#'
#' \item{data}{Simulated data from the population correlation matrix with error
#' (same dimensions, column names, and skew as the input data)}
#'
#' \item{parameters}{Parameters from the \code{pn_object} input (unchanged)}
#'
#' \item{population}{
#' A list containing the population parameters. Elements carried over from
#' the \code{pn_object} input (e.g., \code{Omega}, the true network, and
#' \code{membership} for \code{\link[psychNetwoRk]{simulate_sbm}}) are
#' unchanged except for \code{R}. Added or updated elements:
#'
#' \itemize{
#'
#' \item \code{R} --- Population correlation matrix with error added
#'
#' \item \code{Omega_pseudo} --- Pseudo-true network (constrained maximum
#' likelihood fit of the true network to the population correlation matrix
#' with error)
#'
#' \item \code{Omega_error} --- Network with error (partial correlations of the
#' population correlation matrix with error; saturated in general)
#'
#' \item \code{discrepancy} --- Difference between the pseudo-true and true
#' networks (\code{Omega_pseudo - Omega})
#'
#' }
#' }
#'
#' \item{convergence}{Convergence information from the \code{pn_object} input (unchanged)}
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
#' \code{TLI}, \code{RMSEA}, and \code{SRMR}) achieved (population values for
#' the known graph; see Details)
#'
#' \item \code{delta} --- Minimum of the (maximum likelihood) discrepancy
#' function corresponding to the achieved fit
#'
#' \item \code{error} --- Matrix of the error that was added to the
#' original population correlation matrix
#'
#' \item \code{max_error} --- Largest absolute error that was added
#'
#' \item \code{R_difference} --- Largest absolute difference between the
#' original population correlation matrix and the known graph (constrained
#' maximum likelihood) fit to the population correlation matrix with error.
#' Values near zero indicate that the error is orthogonal to the true network
#' (expected for \code{method = "cb"} only)
#'
#' \item \code{iterations} --- Number of optimization attempts used
#'
#' \item \code{omitted} --- Only for \code{method = "hb"}. Matrix of the omitted
#' edges that were added to the population (partial correlations of the
#' population correlation matrix with error; true network values are zero).
#' The magnitude of each omitted edge is between one half and one times the
#' common size
#'
#' }
#'
#' }
#'
#' \item{original_results}{Original \code{pn_object} input into function}
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
#' sbm_cfi <- add_population_error(pn_object = sbm, method = "cb")
#' \donttest{
#' # Use Wu and Browne to target TLI = 0.90
#' sbm_wb <- add_population_error(
#'   pn_object = sbm, method = "wb",
#'   metric = "TLI", fit_value = 0.90
#' )
#'
#' # Target RMSEA = 0.04
#' sbm_rmsea <- add_population_error(
#'   pn_object = sbm, method = "cb", metric = "RMSEA",
#'   fit_value = 0.04
#' )
#'
#' # Generate small-world data and add population error
#' sw <- simulate_smallworld(
#'   nodes = 20, density = 0.20, rewire = 0.10,
#'   sample_size = 1000
#' )
#'
#' sw_srmr <- add_population_error(
#'   pn_object = sw, method = "hb", fit_value = 0.05
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
# Updated 25.11.2026
add_population_error <- function(
    pn_object,
    method = c("cb", "tkl", "wb", "hb"),
    metric = c("cfi", "tli", "rmsea", "srmr"),
    fit_value = NULL,
    tolerance = 0.01, max_iterations = 10,
    control = list(), seed = NULL
)
{

  # Set seed (reproducible results)
  seed_token <- seed_start(seed)
  on.exit(seed_end(seed_token), add = TRUE)

  # Obtain options specified by the user
  specified <- names(control)

  # Check inputs
  # (`missing()` is evaluated here: it does not carry over to a function
  # that receives an argument with a default)
  inputs <- add_population_error_errors(
    pn_object, method, metric, fit_value,
    tolerance, max_iterations, control,
    method_missing = missing(method)
  )
  method <- inputs$method
  metric <- inputs$metric
  target <- inputs$fit_value
  control <- inputs$control

  # Obtain population parameters
  R <- pn_object$population$R
  network <- pn_object$population$Omega
  nodes <- ncol(R)

  # Identify free cells (non-edges in the upper triangle)
  # Error is placed on cells that the true network constrains to zero
  free_cells <- which(network == 0 & upper.tri(network))

  # Check for free cells
  if(length(free_cells) == 0){
    stop(
      paste(
        "The network in `pn_object` is fully connected so there are no",
        "constrained (zero) edges to add population error to.",
        "\n\nUse a sparser network to add population error"
      )
    )
  }

  # Degrees of freedom of the known graph
  df <- length(free_cells)

  # Check number of omitted edges (cannot exceed the zero edges)
  if(method == "hb" && control$omitted > length(free_cells)){

    # Specified by the user
    if("omitted" %in% specified){
      stop(
        paste0(
          "`control$omitted` = ", control$omitted, " is more than the number ",
          "of edges that are zero in the true network (", length(free_cells), ")",
          "\n\nUse `control = list(omitted = ", length(free_cells), ")` or fewer"
        ), call. = FALSE
      )
    }

    # Default
    control$omitted <- length(free_cells)

  }

  # Original population correlation matrix as the reference
  R_chol <- chol(R)
  reference <- reference_model(R, R_chol)
  inverse_chol <- backsolve(R_chol, diag(nodes))

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
        nodes, df, tolerance = tolerance, inverse_chol = inverse_chol
      ),
      "tkl" = TKL_optimize(
        metric, target, R, network, nodes, df,
        control$minor_factors, control$decay, reference = reference
      ),
      "wb" = WB_optimize(
        metric, target, R, network, nodes, df, reference = reference
      ),
      "hb" = HB_optimize(
        metric, target, R, network, nodes, df,
        free_cells, control$omitted, reference = reference,
        max_condition = control$max_condition
      )
    )

    # Check for invalid fit
    if(!isTRUE(solution$valid)){
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

    # Check for an exact solution (target is matched)
    # (attempts continue past `tolerance` until the target is matched so that
    # the achieved fit is as close to `fit_value` as the error allows)
    if(best_difference <= 1e-06){
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

    warning(
      paste0(
        "Convergence counter has exceeded its limit. ",
        "There were issues converging to the target ", metric, " = ", target,
        " within `tolerance` = ", tolerance, ". \n\n",
        "Using the solution with the closest ", metric, " = ",
        round(best$fit[[metric]], 4), ".",
        swiftelse(
          isTRUE(best$fit$capped),
          paste0(
            "\n\nThe search range was shortened by `control$max_condition` = ",
            control$max_condition, " (largest condition number of the ",
            "population correlation matrix with error), so a larger error may ",
            "reach the target. Increase it (or use `Inf`) to allow larger errors; ",
            "estimation is slower for larger values"
                    ), ""
        )
      ), call. = FALSE
    )

  }

  # Known graph fit to the population correlation matrix with error
  # (from the optimization; for Cudeck and Browne, the known graph is fit once
  # here and verifies that the error is orthogonal to the network)
  known <- best$fit$model
  if(is.null(known)){
    known <- known_graph(R_population, network)
  }
  if(!known$converged){
    warning(
      paste0(
        "The pseudo-true network (`Omega_pseudo`) did not converge. ",
        "Interpret `Omega_pseudo` with caution"
      ), call. = FALSE
    )
  }
  R_difference <- max(abs(known$W - R))

  # Fit indices and discrepancy (population values; the same as `network_fit`
  # of the known graph as the sample size approaches infinity)
  fit_names <- c("CFI", "TLI", "RMSEA", "SRMR")
  fit_vector <- unlist(best$fit[fit_names])
  delta <- best$fit$discrepancy

  # Generate data
  data_output <- simulate_data(
    n = nrow(pn_object$data), p = nodes, R = R_population,
    skew = pn_object$parameters$skew, skew_range = NULL
  )

  # Pseudo-true network (constrained ML fit of the true graph)
  pseudo_true <- known$network
  dimnames(pseudo_true) <- dimnames(R)

  # Network with error (partial correlations of the population correlation
  # matrix with error; saturated in general)
  Omega_error <- cor2pcor(R_population)
  dimnames(Omega_error) <- dimnames(R)

  # Update population
  population <- pn_object$population
  population$R <- R_population
  population$Omega_pseudo <- pseudo_true
  population$Omega_error <- Omega_error
  population$discrepancy <- pseudo_true - network
  dimnames(population$discrepancy) <- dimnames(R)

  # Populate results
  # (same elements and order as `simulate_sbm` and `simulate_smallworld`,
  # followed by the elements specific to population error)
  results <- list(
    data = data_output$data,
    parameters = pn_object$parameters,
    population = population,
    convergence = pn_object$convergence,
    population_error = list(
      method = method,
      metric = metric,
      fit_value = target,
      achieved = unname(fit_vector[metric]),
      fit = fit_vector,
      delta = delta,
      error = R_population - R,
      max_error = max_error,
      R_difference = R_difference,
      iterations = best$iterations
    ),
    original_results = pn_object
  )

  # Add omitted edges (Hu and Bentler)
  if(method == "hb"){
    omitted <- best$fit$omitted
    dimnames(omitted) <- dimnames(R)
    results$population_error$omitted <- omitted
  }

  # Add class
  class(results) <- c(class(pn_object), "pn_pe")

  # Return results
  return(results)

}

# Input checking ----
#' @noRd
# Updated 04.11.2026
add_population_error_errors <- function(
    pn_object, method, metric, fit_value,
    tolerance, max_iterations, control, method_missing = FALSE
)
{

  # Check for appropriate class
  if(!is(pn_object, "pn_simulate")){

    # Produce error
    stop(
      paste(
        "`pn_object` input is not class \"pn_simulate\" from the",
        "`simulate_sbm` or `simulate_smallworld` functions.",
        "\n\nInput class(es) of current `pn_object`:",
        paste0("\"", class(pn_object), "\"", collapse = ", "),
        "\n\nUse `simulate_sbm` or `simulate_smallworld` to generate your",
        "data to input into this function"
      )
    )

  }

  # Check that population error has not yet been added
  if(is(pn_object, "pn_pe")){

    # Produce error
    stop(
      paste(
        "`pn_object` input is class \"pn_pe\" from the `add_population_error` function.",
        "Population error has already been added.",
        "\n\nUse the `original_results` element of this object or generate new data with",
        "`simulate_sbm` or `simulate_smallworld`"
      )
    )

  }

  # Check for method (no default)
  if(method_missing || is.null(method)){
    stop(
      paste0(
        "`method` must be specified. ",
        "\n\nUse `method = \"cb\"`, `method = \"tkl\"`, ",
        "`method = \"wb\"`, or `method = \"hb\"`. ",
        "See `?add_population_error` for details"
      ), call. = FALSE
    )
  }
  typeof_error(method, "character", "add_population_error")
  length_error(method, 1, "add_population_error")
  method <- tolower(method)

  # Check for available methods
  available_methods <- c("cb", "tkl", "wb", "hb")
  if(!method %in% available_methods){
    stop(
      paste0(
        "`method` = \"", method, "\" is not available.",
        "\n\nUse ", paste0("`method = \"", available_methods, "\"`", collapse = ", ")
      ), call. = FALSE
    )
  }

  # Set default metric (consistent with the method) and check for metric
  default_metric <- switch(
    method,
    "cb" = "cfi", "tkl" = "cfi", "wb" = "rmsea", "hb" = "srmr"
  )
  metric <- toupper(set_default(metric, default_metric, add_population_error))

  # Set default fit_value (conventional cutoffs)
  if(is.null(fit_value)){
    fit_value <- switch(
      metric,
      "CFI" = 0.95,
      "TLI" = 0.95,
      "RMSEA" = 0.06,
      "SRMR" = 0.05
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

  # Check tolerance and maximum iterations
  typeof_error(tolerance, c("integer", "double"), "add_population_error")
  length_error(tolerance, 1, "add_population_error")
  if(!is.finite(tolerance) || tolerance <= 0){
    stop("`tolerance` must be a positive number", call. = FALSE)
  }
  typeof_error(max_iterations, c("integer", "double"), "add_population_error")
  length_error(max_iterations, 1, "add_population_error")
  range_error(max_iterations, c(1, Inf), "add_population_error")

  # Check control
  if(!is.list(control)){
    stop("`control` must be a list. See `?add_population_error` for options")
  }

  # Set control defaults
  control_defaults <- list(
    minor_factors = max(50, 2 * ncol(pn_object$population$R)),
    decay = 0.30, omitted = 2, max_condition = 1000
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
  range_error(control$omitted, c(1, Inf), "add_population_error")
  if(control$omitted %% 1 != 0){
    stop(
      paste0("`control$omitted` = ", control$omitted, " must be a whole number"),
      call. = FALSE
    )
  }
  range_error(control$max_condition, c(1, Inf), "add_population_error")

  # Return inputs
  return(
    list(
      method = method, metric = metric,
      fit_value = fit_value, control = control
    )
  )

}

# Debugging ----
# dm <- matrix(0.20, nrow = 3, ncol = 3); diag(dm) <- 0.90
# pn_object <- simulate_sbm(
#   nodes = 6, blocks = 3, density_matrix = dm, sample_size = 1000
# )
# method = "cb"; metric = "CFI"; fit_value = NULL
# tolerance = 0.001; max_iterations = 10; control = list()
