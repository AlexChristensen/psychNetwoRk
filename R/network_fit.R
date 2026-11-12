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
#' with the input \code{network}
#'
#' @param ci Numeric (length = 1).
#' Confidence interval for RMSEA.
#' Defaults to \code{0.95}
#'
#' @return A named numeric vector of traditional and likelihood-based fit
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
#' \code{XX} is the integer percentage (e.g., \code{RMSEA.95.lower} and
#' \code{RMSEA.95.upper} for a 95\% CI)}
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
#' difference between the off-diagonal elements of the model-implied and
#' observed correlation matrices (the diagonal is excluded, so this will be
#' larger than the \pkg{lavaan} value, which includes the diagonal).
#' Values \eqn{\leq} 0.08 indicate acceptable fit}
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
# Updated 11.11.2026
network_fit <- function(network, n, S, ci = 0.95)
{

  # Obtain lower triangle indices
  lower_triangle <- lower.tri(network)

  # Obtain model-implied zero-order correlation matrix
  R <- pcor2cor(network)

  # Compute number of parameters
  p <- dim(network)[2]
  zero_parameters <- p * (p - 1) / 2
  model_parameters <- sum(network[lower_triangle] != 0)

  # Compute baseline
  baseline <- diag(1, nrow = p, ncol = p)
  baseline_ML <- log(det(baseline)) + sum(diag(S %*% solve(baseline))) - log(det(S)) - p
  baseline_chi_square <- n * baseline_ML
  baseline_tli <- baseline_chi_square / zero_parameters

  # Compute traditional SEM measures
  loglik_ML <- log(det(R)) + sum(diag(S %*% solve(R))) - log(det(S)) - p
  chi_square <- n * loglik_ML
  df <- zero_parameters - model_parameters
  chi_max <- max(chi_square - df, 0)
  nDF <- n * df
  rmsea_null <- nDF * 0.0025 # 0.05^2

  # log-likelihood
  loglik <- -(n / 2) * (p * log(2 * pi) + log(det(R)) + sum(diag(S %*% solve(R))))
  # Assumes no mean structure (or that all means are equal to zero)

  # Obtain RMSEA confidence intervals
  if(!is.null(ci)) rmsea_cis <- rmsea_ci(chi_square, df, n, nDF, ci) else rmsea_cis <- NULL

  # Compute indices that are undefined when degrees of freedom are zero
  # (follows {lavaan}: RMSEA = 0, CFI = 1, TLI = 1, p-values are NA)
  if(df > 0){

    chisq_p <- 1 - pchisq(chi_square, df = df)
    rmsea <- sqrt(chi_max / nDF)
    rmsea_p <- 1 - pchisq(chi_square, df = df, ncp = rmsea_null)
    tli <- (baseline_tli - (chi_square / df)) / (baseline_tli - 1)

  }else{

    chisq_p <- NA
    rmsea <- 0
    rmsea_p <- NA
    tli <- 1

  }

  # CFI (denominator follows {lavaan})
  cfi_denominator <- max(chi_square - df, baseline_chi_square - zero_parameters, 0)
  cfi <- swiftelse(cfi_denominator == 0, 1, 1 - (chi_max / cfi_denominator))

  # Get fit indices
  fit_indices <- c(
    # Traditional fit measures
    chisq = chi_square, df = df, chisq.p.value = chisq_p,
    RMSEA = rmsea,
    rmsea_cis,
    RMSEA.p.value = rmsea_p,
    CFI = cfi,
    TLI = tli,
    SRMR = sqrt(mean((R[lower_triangle] - S[lower_triangle])^2)),
    # Gaussian log-likelihood measures
    logLik = loglik,
    AIC = -2 * loglik + 2 * model_parameters,
    BIC = -2 * loglik + model_parameters * log(n)
  )

  # Rename confidence intervals
  names(fit_indices)[names(fit_indices) %in% c("lower", "upper")] <-
    paste("RMSEA", format_integer(ci * 100, 1), c("lower", "upper"), sep = ".")

  # Return log-likelihood
  return(fit_indices)

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
    pchisq(chi_square, df = df, ncp = lambda) - ci
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