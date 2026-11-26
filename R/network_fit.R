#' @title Traditional Fit Metrics for Networks
#'
#' @description Computes several traditional fit metrics for networks including
#'
#' \itemize{
#'
#' \item chi-square (\eqn{\chi^2})
#'
#' \item root mean square error of approximation (RMSEA) with confidence intervals
#'
#' \item confirmatory fit index (CFI)
#'
#' \item Tucker-Lewis index (TLI)
#'
#' \item standardized root mean residual (SRMR)
#'
#' \item log-likelihood
#'
#' \item Akaike's information criterion (AIC)
#'
#' \item Bayesian information criterion (BIC)
#'
#' }
#'
#' @param network Matrix or data frame.
#' A p by p square network matrix
#'
#' @param n Numeric (length = 1).
#' Sample size
#'
#' @param S Matrix or data frame.
#' A p by p square zero-order correlation matrix corresponding
#' with the input \code{network} (must be positive definite)
#'
#' @param ci Numeric (length = 1).
#' Confidence interval for RMSEA (between 0 and 1), or \code{NULL} to skip
#' the confidence interval.
#' Defaults to \code{0.95}
#'
#' @return A named numeric vector (class \code{"pn_fit"}, with \code{print}
#' and \code{summary} methods; indexing such as \code{fit["CFI"]} returns plain
#' numbers) of traditional and likelihood-based fit
#' indices. The vector always contains the following elements:
#'
#' \describe{
#'
#' \item{\code{chisq}}{Chi-square statistic
#' (\eqn{\chi^2 = n \cdot F_{ML}}), where \eqn{F_{ML}} is the maximum
#' likelihood discrepancy function between the model-implied and empirical
#' correlation matrices}
#'
#' \item{\code{df}}{Degrees of freedom: total number of unique off-diagonal
#' correlations minus the number of non-zero edges in \code{network}}
#'
#' \item{\code{chisq.p.value}}{p-value for the chi-square test of exact fit
#' (H0: model-implied covariance equals the population covariance)}
#'
#' \item{\code{RMSEA}}{Root mean square error of approximation (0 when
#' \code{df = 0}; the p-values are \code{NA} and CFI and TLI are 1 in that
#' case, following \pkg{lavaan}). Values \eqn{\leq} 0.05 indicate close
#' fit; values \eqn{\leq} 0.08 indicate acceptable fit}
#'
#' \item{\code{RMSEA.XX.lower}, \code{RMSEA.XX.upper}}{Lower and upper
#' bounds of the \code{ci}-level confidence interval for RMSEA, where
#' \code{XX} is the percentage (e.g., \code{RMSEA.95.lower} and
#' \code{RMSEA.95.upper} for a 95\% CI and \code{RMSEA.97.5.lower} for a
#' 97.5\% CI)}
#'
#' \item{\code{RMSEA.p.value}}{p-value for the one-sided test of close fit
#' (H0: RMSEA \eqn{\leq} 0.05)}
#'
#' \item{\code{CFI}}{Comparative fit index, comparing the target model to
#' an independence (null) baseline. Values \eqn{\geq} 0.95 indicate
#' acceptable fit}
#'
#' \item{\code{TLI}}{Tucker-Lewis index (non-normed fit index). Values
#' \eqn{\geq} 0.95 indicate acceptable fit; can fall outside \eqn{[0, 1]}
#' for severely misspecified models}
#'
#' \item{\code{SRMR}}{Standardized root mean residual: the root mean squared
#' difference between the lower triangle elements (including the diagonal) of
#' the model-implied and observed correlation matrices, following
#' \pkg{lavaan}. Values \eqn{\leq} 0.08 indicate acceptable fit}
#'
#' \item{\code{logLik}}{Gaussian log-likelihood of the model-implied
#' correlation matrix, assuming zero mean structure (means are not
#' estimated)}
#'
#' \item{\code{AIC}}{Akaike's information criterion:
#' \eqn{-2 \cdot \log L + 2 \cdot E}, where \eqn{E} is the number of
#' non-zero edges in \code{network}}
#'
#' \item{\code{BIC}}{Bayesian information criterion:
#' \eqn{-2 \cdot \log L + E \cdot \log(n)}}
#'
#' }
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Simulate small-world network data
#' sw <- simulate_smallworld(
#'   nodes = 10, density = 0.30, rewire = 0.20,
#'   sample_size = 500
#' )
#'
#' # Obtain correlation matrix
#' S <- cor(sw$data)
#'
#' # Obtain fit of the true network
#' # (expects continuous variables!)
#' network_fit(network = sw$population$Omega, n = nrow(sw$data), S = S)
#'
#' @references
#' Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017).
#' Generalized network psychometrics: Combining network and latent variable models.
#' \emph{Psychometrika}, \emph{82}(4), 904–927.
#'
#' @export
#'
# Compute network fit statistics ----
# Updated 20.11.2026
network_fit <- function(network, n, S, ci = 0.95)
{

  # Check inputs
  inputs <- network_fit_errors(network, n, S, ci)
  network <- inputs$network; S <- inputs$S

  # Obtain lower triangle indices
  lower_triangle <- lower.tri(network)
  srmr_triangle <- lower.tri(network, diag = TRUE) # SRMR includes the diagonal

  # Obtain model-implied zero-order correlation matrix
  R <- pcor2cor(network)

  # Compute number of parameters
  p <- dim(network)[2]
  zero_parameters <- p * (p - 1) / 2
  model_parameters <- sum(network[lower_triangle] != 0)
  df <- zero_parameters - model_parameters

  # Decompositions
  R_chol <- chol(R)
  R_inverse <- chol2inv(R_chol)
  log_det_R <- 2 * sum(log(diag(R_chol)))
  log_det_S <- 2 * sum(log(diag(chol(S))))
  trace_S_R <- sum(S * R_inverse)

  # Baseline (independence model with unit variances)
  baseline_chi_square <- n * (sum(diag(S)) - log_det_S - p)

  # Model discrepancy and chi-square
  chi_square <- n * (log_det_R + trace_S_R - log_det_S - p)
  nDF <- n * df
  rmsea_null <- nDF * 0.0025 # 0.05^2

  # log-likelihood
  # (assumes no mean structure or that all means are equal to zero)
  loglik <- -(n / 2) * (p * log(2 * pi) + log_det_R + trace_S_R)

  # Obtain RMSEA confidence intervals
  if(!is.null(ci)) rmsea_cis <- rmsea_ci(chi_square, df, n, nDF, ci) else rmsea_cis <- NULL

  # Fit indices
  indices <- sem_indices(chi_square, df, baseline_chi_square, zero_parameters, n)

  # p-values (undefined when degrees of freedom are zero)
  if(df > 0){
    chisq_p <- pchisq(chi_square, df = df, lower.tail = FALSE)
    rmsea_p <- suppressWarnings(
      pchisq(chi_square, df = df, ncp = rmsea_null, lower.tail = FALSE)
    )
  }else{
    chisq_p <- NA
    rmsea_p <- NA
  }

  # Get fit indices
  fit_indices <- c(
    # Traditional fit measures
    chisq = chi_square, df = df, chisq.p.value = chisq_p,
    RMSEA = indices[["RMSEA"]],
    rmsea_cis,
    RMSEA.p.value = rmsea_p,
    CFI = indices[["CFI"]],
    TLI = indices[["TLI"]],
    SRMR = sqrt(mean((R[srmr_triangle] - S[srmr_triangle])^2)),
    # Gaussian log-likelihood measures
    logLik = loglik,
    AIC = -2 * loglik + 2 * model_parameters,
    BIC = -2 * loglik + model_parameters * log(n)
  )

  # Rename confidence intervals (e.g., 95 or 97.5)
  names(fit_indices)[names(fit_indices) %in% c("lower", "upper")] <-
    paste("RMSEA", format(ci * 100, trim = TRUE), c("lower", "upper"), sep = ".")

  # Add class
  class(fit_indices) <- "pn_fit"

  # Return fit indices
  return(fit_indices)

}

#' @noRd
# Errors ----
# Updated 20.11.2026
network_fit_errors <- function(network, n, S, ci)
{

  # Check network
  object_error(network, c("matrix", "data.frame", "tibble"), "network_fit")
  network <- as.matrix(network)
  typeof_error(network, "numeric", "network_fit")
  if(!is_symmetric(network)){
    stop("`network` must be a square, symmetric matrix", call. = FALSE)
  }

  # Check correlation matrix
  object_error(S, c("matrix", "data.frame", "tibble"), "network_fit")
  S <- as.matrix(S)
  typeof_error(S, "numeric", "network_fit")
  if(!is_symmetric(S) || !all(dim(S) == dim(network))){
    stop("`S` must be a square, symmetric matrix with the same dimensions as `network`", call. = FALSE)
  }
  if(is(try(chol(S), silent = TRUE), "try-error")){
    stop("`S` must be positive definite", call. = FALSE)
  }

  # Check sample size
  typeof_error(n, "numeric", "network_fit")
  length_error(n, 1, "network_fit")
  range_error(n, c(1, Inf), "network_fit")

  # Check confidence interval
  if(!is.null(ci)){
    typeof_error(ci, "numeric", "network_fit")
    length_error(ci, 1, "network_fit")
    if(ci <= 0 || ci >= 1){
      stop("`ci` must be between 0 and 1 (or NULL)", call. = FALSE)
    }
  }

  # Return inputs
  return(list(network = network, S = S))

}

#' @noRd
# CFI, TLI, and RMSEA from chi-squares ----
# Follows {lavaan}: RMSEA = 0 and TLI = 1 when the degrees of freedom are
# zero; CFI = 1 when its denominator is zero (CFI denominator follows {lavaan})
# Updated 20.11.2026
sem_indices <- function(chi_square, df, baseline_chi_square, baseline_df, n)
{

  # Chi-square above degrees of freedom
  model_excess <- max(chi_square - df, 0)
  denominator <- max(chi_square - df, baseline_chi_square - baseline_df, 0)

  # Indices
  if(df > 0){
    baseline_ratio <- baseline_chi_square / baseline_df
    RMSEA <- sqrt(model_excess / (n * df))
    TLI <- (baseline_ratio - (chi_square / df)) / (baseline_ratio - 1)
  }else{
    RMSEA <- 0
    TLI <- 1
  }

  # Return indices
  return(
    c(
      CFI = swiftelse(denominator == 0, 1, 1 - model_excess / denominator),
      TLI = TLI, RMSEA = RMSEA
    )
  )

}

#' @noRd
# RMSEA Confidence Intervals ----
# Follows {lavaan} version 0.6.19
# Updated 01.11.2024
rmsea_ci <- function(chi_square, df, n, nDF, ci)
{

  # Set up CI
  lower_ci <- 1 - (1 - ci) / 2
  upper_ci <- 1 - lower_ci

  # Internal function for finding RMSEA confidence intervals
  # (same as {lavann} version 0.6.19)
  find_lambda <- function(lambda, ci){
    suppressWarnings(pchisq(chi_square, df = df, ncp = lambda)) - ci
  }

  # Find lower bound
  if(df < 1 || find_lambda(0, lower_ci) < 0){
    rmsea_lower <- 0
  }else{

    # Try uniroot
    lambda <- try(
      uniroot(f = find_lambda, lower = 0, upper = chi_square, ci = lower_ci)$root,
      silent = TRUE
    )

    # Determine if there was an error
    rmsea_lower <- swiftelse(is(lambda, "try-error"), NA, sqrt(lambda / nDF))

  }

  # Find upper bound
  N_RMSEA <- max(n, chi_square * 4)
  if(df < 1 || find_lambda(N_RMSEA, upper_ci) > 0 || find_lambda(0, upper_ci) < 0){
    rmsea_upper <- 0
  }else{

    # Try uniroot
    lambda <- try(
      uniroot(f = find_lambda, lower = 0, upper = N_RMSEA, ci = upper_ci)$root,
      silent = TRUE
    )

    # Determine if there was an error
    rmsea_upper <- swiftelse(is(lambda, "try-error"), NA, sqrt(lambda / nDF))

  }

  # Return confidence interval
  return(c(lower = rmsea_lower, upper = rmsea_upper))

}