#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### psychNetwoRk S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Updated 03.11.2026

# print() Methods ----

# Print `ggm_simulate`
# Handles output from both `simulate_sbm` and `simulate_smallworld`,
# which share this class (and `add_population_error`, which adds `ggm_pe`)
# Updated 03.11.2026
#' @export
print.ggm_simulate <- function(x, ...)
{

  # Print compact overview
  ggm_simulate_describe(x, detailed = FALSE)

  # Return input invisibly
  invisible(x)

}

# summary() Methods ----

# Summary `ggm_simulate`
# Updated 03.11.2026
#' @export
summary.ggm_simulate <- function(object, ...)
{

  # Print detailed overview
  ggm_simulate_describe(object, detailed = TRUE)

  # Return input invisibly
  invisible(object)

}

# Shared helpers ----

#' @noRd
# Formats the range and mean of a numeric vector for `ggm_simulate` print/summary
# Updated 03.11.2026
ggm_simulate_format_range <- function(values)
{
  sprintf(
    "%.3f to %.3f (mean = %.3f)",
    min(values), max(values), mean(values)
  )
}

#' @noRd
# Shared print/summary body for `ggm_simulate` objects
# Updated 03.11.2026
ggm_simulate_describe <- function(x, detailed = FALSE)
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
  cat(
    if(isTRUE(sbm)){
      "Simulated Stochastic Block Model GGM Data\n"
    }else{
      "Simulated Small-World GGM Data\n"
    }
  )
  cat(strrep("-", 40), "\n", sep = "")

  cat(sprintf("%-24s%d\n", "Sample size:", dimensions[1]))
  cat(sprintf("%-24s%d\n", "Nodes:", dimensions[2]))

  if(isTRUE(sbm)){
    cat(sprintf("%-24s%d\n", "Blocks:", parameters$blocks))
    cat(
      sprintf("%-24s", "Nodes per block:"),
      paste(parameters$nodes, collapse = ", "), "\n", sep = ""
    )
  }else{
    cat(sprintf("%-24s%.3f\n", "Rewire:", parameters$rewire))
    cat(sprintf("%-24s%d\n", "Neighbors:", parameters$neighbors))
  }

  cat(sprintf("%-24s%.3f\n", "Density:", length(edges) / possible_edges))
  cat(sprintf("%-24s%.3f\n", "Negative proportion:", mean(edges < 0)))
  cat(sprintf("%-24s%s\n", "Edge weights:", ggm_simulate_format_range(abs(edges))))
  cat(sprintf("%-24s%.3f\n", "Omega:", parameters$omega))
  cat(sprintf("%-24s%.3f\n", "Modularity (Q):", parameters$Q))

  # Skew
  if(any(parameters$skew != 0)){
    cat(sprintf("%-24s%s\n", "Skew:", ggm_simulate_format_range(parameters$skew)))
  }

  # Population error
  if(is(x, "ggm_pe")){

    error <- x$population_error

    cat("\nPopulation error:\n")
    cat(sprintf("  %-24s%.3f\n", "Target CFI:", error$target_CFI))
    cat(sprintf("  %-24s%.3f\n", "CFI:", error$CFI))
    cat(sprintf("  %-24s%.5f\n", "Delta:", error$delta))

    if(isTRUE(detailed)){
      error_values <- error$error[upper.tri(error$error)]
      error_values <- error_values[error_values != 0]
      cat(sprintf("  %-24s%s\n", "Error (absolute):", ggm_simulate_format_range(abs(error_values))))
      cat(sprintf("  %-24s%.2e\n", "Orthogonality (max):", error$R_difference))
      cat(sprintf("  %-24s%d\n", "Optimization attempts:", error$iterations))
    }

  }

  if(isTRUE(detailed)){
    cat("\nPopulation network (Omega):\n")
    print(round(network, 3))
  }

}
