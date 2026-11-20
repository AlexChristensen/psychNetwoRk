#' Predict Weibull Parameters for Edge Weight Distributions
#'
#' @description
#' Predicts the shape and scale parameters of a Weibull distribution that
#' characterizes the absolute partial correlation edge weights of a psychometric
#' network, given its number of nodes and signal-to-noise ratio. Parameter estimates are
#' derived from a Seemingly Unrelated Regression (SUR) model fitted to empirical
#' network data from Huth et al. (2025), where absolute partial correlations were
#' found to follow a Weibull distribution more consistently than Beta, Gamma, or
#' log-normal alternatives.
#'
#' @param nodes Numeric (length = 1).
#' The number of nodes (variables) in the network (at least 2).
#' Values outside of 8 to 54, the range of the empirical networks used to fit
#' the underlying SUR model, extrapolate the model and give a warning.
#'
#' @param snr Numeric (length = 1).
#' Signal-to-noise ratio of partial correlations (\eqn{\bar{|w|} / \mathrm{SD}(|w|)});
#' must be positive.
#' Values less than 1 indicate wider range of partial correlations (\eqn{w}) whereas
#' values greater than 1 indicate narrower range.
#' Defaults to \code{1} where the mean of the partial correlations (\eqn{\bar{|w|}})
#' is equal to the standard deviation (\eqn{\mathrm{SD}(|w|)})
#'
#' @param seed Numeric (length = 1).
#' Seed for reproducible results when \code{bootstrap = TRUE}. The same
#' \code{seed} gives the same results (and does not change the random number
#' stream of the R session). Defaults to \code{NULL} (random seed)
#'
#' @param bootstrap Boolean (length = 1). If \code{TRUE}, a randomly sampled residual from
#' the SUR model fit is added to each predicted parameter, introducing
#' empirically grounded variability suitable for use in simulation or
#' bootstrapping contexts. Defaults to \code{FALSE}.
#'
#' @details
#' The shape and scale parameters are predicted from two network descriptors
#' used in both SUR equations:
#'
#' \describe{
#'   \item{\code{snr}}{The signal-to-noise ratio of the absolute partial
#'     correlations, \eqn{\bar{|w|} / \mathrm{SD}(|w|)}.}
#'   \item{\code{rlp}}{The reciprocal log of the number of nodes,
#'     \eqn{1 / \log(p)}, capturing the diminishing marginal effect of network
#'     size on edge weight distributions.}
#' }
#'
#' Sample size is not a predictor: the predicted parameters describe the
#' population edge weight distribution, independent of sampling precision.
#'
#' These predictors enter two SUR equations whose coefficients are stored in
#' the internal \code{weibull_weights} dataset. SUR was used to account for the
#' correlated residuals between the shape and scale equations across networks
#' (residual correlation = 0.330). Model fit was strong: the shape equation
#' achieved \eqn{R^2 = 0.887} (RMSE = 0.048) and the scale equation achieved
#' \eqn{R^2 = 0.885} (RMSE = 0.014). Shape residuals were normally distributed
#' (Shapiro-Wilk W = 0.990, p = 0.173). Scale residuals showed a modest
#' departure from normality (Shapiro-Wilk W = 0.985, p = 0.034), consistent
#' with slight right skew in the scale outcome and test sensitivity at n = 194
#' rather than a substantive violation. Heteroskedasticity was detected in both
#' the shape equation (Breusch-Pagan p < 0.001) and the scale equation
#' (Breusch-Pagan p < 0.001); robust standard errors (HC3) were used for
#' inference. Multicollinearity among predictors in the scale equation was
#' negligible (VIF \eqn{\leq} 1.60); the shape equation contains only two
#' predictors with no multicollinearity concern.
#'
#' \code{nodes} influences predicted shape and scale via \code{rlp}, and
#' \code{snr} influences both parameters directly.
#'
#' Empirically, shape values ranged from approximately 0.72 to 1.63
#' (M = 1.07, SD = 0.14) and scale values from approximately 0.03 to 0.19
#' (M = 0.10, SD = 0.03) across the 194 networks used to fit the model. Shape
#' values near 1 indicate approximately exponential edge weight distributions;
#' values above 1 indicate a rising hazard (mode-bearing distribution).
#'
#' When \code{bootstrap = TRUE}, residuals are drawn via \code{shuffle()} --
#' a random sampling without replacement -- from the empirical SUR residuals,
#' preserving the observed marginal residual distribution.
#'
#' @return A named numeric vector of length 2:
#' \describe{
#'   \item{\code{shape}}{The predicted Weibull shape parameter (\eqn{k > 0}).}
#'   \item{\code{scale}}{The predicted Weibull scale parameter (\eqn{\lambda > 0}).}
#' }
#'
#' @examples
#' # Predict parameters for a 10-node network
#' weibull_parameters(nodes = 10)
#'
#' # With bootstrapped residuals for use in simulation
#' weibull_parameters(nodes = 10, bootstrap = TRUE)
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., &
#' Marsman, M. (2025). Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#' \doi{10.1038/s41562-025-02314-2}
#'
#' @export
#'
# Predict Weibull parameters ----
# Updated 18.11.2026
weibull_parameters <- function(nodes, snr = 1, bootstrap = FALSE, seed = NULL)
{

  # Set seed (reproducible results)
  seed_token <- seed_start(seed)
  on.exit(seed_end(seed_token), add = TRUE)

  # Check inputs
  typeof_error(nodes, "numeric", "weibull_parameters")
  length_error(nodes, 1, "weibull_parameters")
  range_error(nodes, c(2, Inf), "weibull_parameters")
  typeof_error(snr, "numeric", "weibull_parameters")
  length_error(snr, 1, "weibull_parameters")
  range_error(snr, c(.Machine$double.eps, Inf), "weibull_parameters")
  typeof_error(bootstrap, "logical", "weibull_parameters")
  length_error(bootstrap, 1, "weibull_parameters")

  # Warn outside of the range of the empirical networks
  if(nodes < 8 || nodes > 54){
    warning(
      paste0(
        "`nodes` (", nodes, ") is outside of the range of the empirical networks ",
        "used to fit the model (8-54); parameters are extrapolated"
      ), call. = FALSE
    )
  }

  # Get Weibull model
  weibull_weights <- get(data("weibull_weights", package = "psychNetwoRk", envir = environment()))

  # Compute descriptive parameters
  scale_parameters <- shape_parameters <- c(snr = snr, rlp = 1 / log(nodes))

  # Compute Weibull parameters (coefficients matched by name)
  shape_coefficients <- weibull_weights$shape$coefficients
  scale_coefficients <- weibull_weights$scale$coefficients
  shape <- unname(
    shape_coefficients[1] + sum(shape_coefficients[-1] * shape_parameters[sub("^shape_", "", names(shape_coefficients)[-1])])
  )
  scale <- unname(
    scale_coefficients[1] + sum(scale_coefficients[-1] * scale_parameters[sub("^scale_", "", names(scale_coefficients)[-1])])
  )

  # Bootstrap residuals
  if(bootstrap){

    # Draw index
    index <- shuffle(seq_along(weibull_weights$shape$residuals), 1)

    # Set parameters with residuals
    shape <- shape + weibull_weights$shape$residuals[index]
    scale <- scale + weibull_weights$scale$residuals[index]

  }

  # Return parameters
  return(c(shape = shape, scale = scale))

}
