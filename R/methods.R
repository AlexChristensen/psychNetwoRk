#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### psychNetwoRk S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Conventions (shared with {latentFactoR}): a title, a rule of dashes,
# left-justified labels (`%-24s`), values with three decimals, and
# `print` and `summary` both return the input invisibly. `summary` adds
# details (matrices, fit indices) to the `print` overview

# Updated 20.11.2026

# pn_simulate ----

# Print `pn_simulate`
# Handles output from both `simulate_sbm` and `simulate_smallworld`,
# which share this class (and `add_population_error`, which adds `pn_pe`)
# Updated 20.11.2026
#' @export
print.pn_simulate <- function(x, ...)
{

  # Print compact overview
  pn_simulate_describe(x, detailed = FALSE)

  # Return input invisibly
  invisible(x)

}

# Summary `pn_simulate`
# Updated 20.11.2026
#' @export
summary.pn_simulate <- function(object, ...)
{

  # Print detailed overview
  pn_simulate_describe(object, detailed = TRUE)

  # Return input invisibly
  invisible(object)

}

# pn_fit ----

# Print `pn_fit` (`network_fit`)
# Updated 20.11.2026
#' @export
print.pn_fit <- function(x, ...)
{

  # Header
  pn_header("Network Fit")

  # Chi-square
  cat(
    sprintf(
      "%-24s%.3f (df = %d, p = %s)\n", "Chi-square:", x[["chisq"]], as.integer(x[["df"]]),
      pn_format_p(x[["chisq.p.value"]])
    )
  )

  # RMSEA with confidence interval
  ci_names <- grep("^RMSEA\\..*\\.(lower|upper)$", names(x), value = TRUE)
  rmsea_text <- sprintf("%.3f", x[["RMSEA"]])
  if(length(ci_names) == 2){
    rmsea_text <- sprintf(
      "%s [%.3f, %.3f] (%s%% CI)", rmsea_text, x[[ci_names[1]]], x[[ci_names[2]]],
      sub("^RMSEA\\.(.*)\\.lower$", "\\1", ci_names[1])
    )
  }
  cat(sprintf("%-24s%s\n", "RMSEA:", rmsea_text))

  # Incremental and residual fit
  cat(sprintf("%-24s%.3f\n", "CFI:", x[["CFI"]]))
  cat(sprintf("%-24s%.3f\n", "TLI:", x[["TLI"]]))
  cat(sprintf("%-24s%.3f\n", "SRMR:", x[["SRMR"]]))

  # Return input invisibly
  invisible(x)

}

# Summary `pn_fit`
# Updated 20.11.2026
#' @export
summary.pn_fit <- function(object, ...)
{

  # Overview
  print.pn_fit(object)

  # Additional measures
  cat(sprintf("%-24s%s\n", "RMSEA p-value (close):", pn_format_p(object[["RMSEA.p.value"]])))
  cat("\nLikelihood-based measures:\n")
  cat(sprintf("  %-24s%.3f\n", "Log-likelihood:", object[["logLik"]]))
  cat(sprintf("  %-24s%.3f\n", "AIC:", object[["AIC"]]))
  cat(sprintf("  %-24s%.3f\n", "BIC:", object[["BIC"]]))

  # Return input invisibly
  invisible(object)

}

# Shared helpers ----

#' @noRd
# Print a title with a rule of dashes
# Updated 20.11.2026
pn_header <- function(title)
{
  cat(title, "\n", sep = "")
  cat(strrep("-", 40), "\n", sep = "")
}

#' @noRd
# Section header
# Updated 20.11.2026
pn_section <- function(title)
{
  cat("\n", title, ":\n", sep = "")
}

#' @noRd
# Format p-values
# Updated 20.11.2026
pn_format_p <- function(p)
{
  if(is.na(p)){
    return("NA")
  }else if(p < 0.001){
    return("< .001")
  }else{
    return(sprintf("%.3f", p))
  }
}

#' @noRd
# Formats the range and mean of a numeric vector for print/summary
# (a single value is printed once)
# Updated 20.11.2026
pn_simulate_format_range <- function(values)
{

  # Single value
  if(length(unique(values)) == 1){
    return(sprintf("%.3f (all)", values[1]))
  }

  # Range and mean
  sprintf(
    "%.3f to %.3f (mean = %.3f)",
    min(values), max(values), mean(values)
  )

}

#' @noRd
# Print a matrix (rounded) or, for large matrices, a summary of its values
# Updated 20.11.2026
pn_print_matrix <- function(matrix, title, maximum = 20)
{

  # Large matrices
  if(ncol(matrix) > maximum){
    values <- matrix[lower.tri(matrix)]
    values <- values[values != 0]
    cat(
      sprintf(
        "\n%s (%d nodes; non-zero values): %s\n", title, ncol(matrix),
        swiftelse(length(values) > 0, pn_simulate_format_range(values), "none")
      )
    )
  }else{
    cat("\n", title, ":\n", sep = "")
    print(round(matrix, 3))
  }

}

#' @noRd
# Shared print/summary body for `pn_simulate` objects
# Updated 20.11.2026
pn_simulate_describe <- function(x, detailed = FALSE)
{

  parameters <- x$parameters
  network <- x$population$Omega
  sbm <- !is.null(parameters$blocks)
  dimensions <- dim(x$data)

  # Edge information
  edges <- network[lower.tri(network)]
  edges <- edges[edges != 0]
  possible_edges <- length(network[lower.tri(network)])

  # Header
  pn_header(
    swiftelse(
      isTRUE(sbm), "Simulated Stochastic Block Model GGM Data",
      "Simulated Small-World GGM Data"
    )
  )

  cat(sprintf("%-24s%d\n", "Sample size:", dimensions[1]))
  cat(sprintf("%-24s%d\n", "Nodes:", dimensions[2]))

  if(isTRUE(sbm)){
    cat(sprintf("%-24s%d\n", "Blocks:", as.integer(parameters$blocks)))
    cat(
      sprintf("%-24s", "Nodes per block:"),
      paste(parameters$nodes, collapse = ", "), "\n", sep = ""
    )
  }else{
    cat(sprintf("%-24s%.3f\n", "Rewire:", parameters$rewire))
    cat(sprintf("%-24s%d\n", "Neighbors:", as.integer(parameters$neighbors)))
  }

  cat(sprintf("%-24s%.3f\n", "Density:", length(edges) / possible_edges))
  cat(
    sprintf(
      "%-24s%.3f (requested %.3f%s)\n", "Negative edges:", mean(edges < 0),
      parameters$negative_proportion,
      swiftelse(isTRUE(sbm), " of between-community edges", "")
    )
  )
  cat(sprintf("%-24s%s\n", "Edge weights:", pn_simulate_format_range(abs(edges))))
  cat(sprintf("%-24s%.3f\n", "Omega:", parameters$omega))
  cat(sprintf("%-24s%.3f\n", "Modularity (Q):", parameters$Q))

  # Skew
  if(any(parameters$skew != 0)){
    cat(sprintf("%-24s%s\n", "Skew:", pn_simulate_format_range(parameters$skew)))
  }

  # Detailed: generation diagnostics
  if(isTRUE(detailed)){
    cat(
      sprintf(
        "%-24sshape = %.3f, scale = %.3f\n", "Weibull:",
        parameters$weibull[["shape"]], parameters$weibull[["scale"]]
      )
    )
    cat(sprintf("%-24s%.3f\n", "Condition number:", x$convergence$condition))
    if(!is.na(x$convergence$lambda)){
      cat(sprintf("%-24s%.5f\n", "Ridge (lambda):", x$convergence$lambda))
    }
    cat(sprintf("%-24s%d\n", "Attempts:", as.integer(x$convergence$iterations)))
  }

  # Population error
  if(is(x, "pn_pe")){
    pn_describe_pe(x, detailed)
  }

  # Detailed: network
  if(isTRUE(detailed)){
    pn_print_matrix(network, "Population network (Omega)")
    if(is(x, "pn_pe")){
      pn_print_matrix(
        x$population$Omega_pseudo,
        "Pseudo-true network (Omega_pseudo)"
      )
    }
  }

}

#' @noRd
# Population error (`add_population_error`)
# Updated 20.11.2026
pn_describe_pe <- function(x, detailed = FALSE)
{

  # Obtain population error
  error <- x$population_error

  # Describe
  pn_section("Population error")
  cat(sprintf("  %-24s%s\n", "Method:", toupper(error$method)))
  cat(sprintf("  %-24s%.3f\n", paste0("Target ", error$metric, ":"), error$fit_value))
  cat(sprintf("  %-24s%.3f\n", paste0("Achieved ", error$metric, ":"), error$achieved))
  cat(sprintf("  %-24s%.5f\n", "Delta:", error$delta))
  cat(sprintf("  %-24s%.3f\n", "Maximum error:", error$max_error))

  # Detailed fit
  if(isTRUE(detailed)){
    error_values <- error$error[upper.tri(error$error)]
    error_values <- error_values[error_values != 0]
    cat(sprintf("  %-24s%s\n", "Error (absolute):", pn_simulate_format_range(abs(error_values))))
    cat(
      sprintf(
        "  %-24s%s\n", "Fit indices:",
        paste(names(error$fit), sprintf("%.3f", error$fit), collapse = ", ")
      )
    )
    if(error$method == "cb"){
      cat(sprintf("  %-24s%.2e\n", "Orthogonality (max):", error$R_difference))
    }
    cat(sprintf("  %-24s%d\n", "Optimization attempts:", as.integer(error$iterations)))
  }

}
