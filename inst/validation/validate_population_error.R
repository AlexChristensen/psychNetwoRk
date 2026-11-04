#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Validation of `add_population_error` Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Run from the package root:
#   Rscript inst/validation/validate_population_error.R
#
# Environment variables:
#   VALIDATE_CORES = number of cores for parallel checks (default 10)
#   VALIDATE_QUICK = 1 for fewer Monte Carlo draws and replications
#   VALIDATE_CSV   = path to save every check (method, check, pass, detail)
#
# Standard (acceptable) cut-offs from Hu and Bentler (1998):
#   CFI = TLI = 0.95, RMSEA = 0.06, SRMR = 0.08
#
# Every method is checked for (1) the common contract (target reached,
# fit indices correct against an independent implementation, valid matrices,
# `max_residual`, data-level behavior, input validation) and (2) the specific
# way in which the method is described to add error (see sections below)
#
# Note: the package random number generators are seeded from the clock and
# process ID, so networks differ between runs. Checks are therefore
# written as properties that hold for any network

# Settings ----
quick <- identical(Sys.getenv("VALIDATE_QUICK"), "1")
cores <- as.integer(Sys.getenv("VALIDATE_CORES", "10"))
draws <- if(quick){100}else{300} # Monte Carlo draws
reps <- if(quick){3}else{5} # replications for random methods
cutoffs <- c(CFI = 0.95, TLI = 0.95, RMSEA = 0.06, SRMR = 0.08)
methods <- c("cb", "tkl", "wb", "hb")
metrics <- names(cutoffs)
tolerance <- 0.001
start_time <- Sys.time()

# Load package ----
suppressMessages(library(parallel))
if(file.exists("DESCRIPTION") && any(grepl("^Package: psychNetwoRk", readLines("DESCRIPTION")))){
  suppressMessages(devtools::load_all(".", quiet = TRUE))
}else{
  suppressMessages(library(psychNetwoRk))
}

# Internal function access (works with `load_all` and installed package)
internal <- function(name)
{
  return(get(name, envir = asNamespace("psychNetwoRk")))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Logging and parallel helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Results log ----
new_log <- function()
{

  # Initialize rows
  rows <- list()

  # Add a check
  add <- function(method, check, pass, detail = ""){
    rows[[length(rows) + 1]] <<- data.frame(
      method = method, check = check,
      pass = isTRUE(as.logical(pass)),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  # Get checks
  get <- function(){
    if(length(rows) == 0){
      return(NULL)
    }
    return(do.call(rbind, rows))
  }

  return(list(add = add, get = get))

}

# Format scientific numbers ----
sci <- function(x)
{
  return(sprintf("%.2e", x))
}

# Run `add_population_error` capturing warnings and errors ----
run_pe <- function(g, ...)
{

  # Initialize warnings
  warns <- character()

  # Run
  out <- withCallingHandlers(
    tryCatch(add_population_error(g, ...), error = function(e) e),
    warning = function(w){
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  return(list(out = out, warns = warns))

}

# Run a list of jobs in parallel (each returns rows from a log) ----
run_jobs <- function(specs)
{

  # Run job with error capture
  run_spec <- function(spec){
    tryCatch(
      do.call(spec$fun, spec$args),
      error = function(e){
        data.frame(
          method = spec$label, check = "job completed", pass = FALSE,
          detail = conditionMessage(e), stringsAsFactors = FALSE
        )
      }
    )
  }

  # Run
  results <- mclapply(specs, run_spec, mc.cores = cores, mc.preschedule = FALSE)

  # Replace failures from the parallel process itself
  results <- lapply(seq_along(results), function(i){
    if(is.data.frame(results[[i]])){
      return(results[[i]])
    }
    data.frame(
      method = specs[[i]]$label, check = "job completed", pass = FALSE,
      detail = "parallel worker failed", stringsAsFactors = FALSE
    )
  })

  return(do.call(rbind, results))

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Independent fit-index implementation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Constrained ML (known graph) refit with tight convergence ----
ind_refit <- function(Rp, net)
{

  zeros <- which(net == 0 & upper.tri(net), arr.ind = TRUE)

  # No constraints: saturated model reproduces the matrix
  if(nrow(zeros) == 0){
    return(Rp)
  }

  W <- suppressWarnings(
    glasso::glasso(
      Rp, rho = 0, zero = zeros, thr = 1e-10, maxit = 10000,
      penalize.diagonal = FALSE
    )$w
  )

  return((W + t(W)) / 2)

}

# Log-determinant ----
ind_logdet <- function(M)
{
  return(as.numeric(determinant(M, logarithm = TRUE)$modulus))
}

# ML discrepancy: F = log|S| + tr(P S^-1) - log|P| - p ----
ind_discrepancy <- function(Sigma, Rp)
{
  return(ind_logdet(Sigma) + sum(diag(Rp %*% solve(Sigma))) - ind_logdet(Rp) - ncol(Rp))
}

# Fit indices (Hu & Bentler, 1998, Table 1) at the population level ----
# `n` = 8.3e09 is the population-level constant used by the package
ind_indices <- function(Rp, Sigma, df, n = 8.3e09)
{

  p <- ncol(Rp)
  Fm <- ind_discrepancy(Sigma, Rp)
  chi <- n * Fm
  chi_b <- -n * ind_logdet(Rp) # independence model (identity)
  df_b <- p * (p - 1) / 2
  residual <- (Rp - Sigma)[lower.tri(Rp, diag = TRUE)]

  return(
    c(
      CFI = 1 - max(chi - df, 0) / max(chi_b - df_b, chi - df, 0),
      TLI = ((chi_b / df_b) - (chi / df)) / ((chi_b / df_b) - 1),
      RMSEA = sqrt(max(chi - df, 0) / (df * n)),
      SRMR = sqrt(mean(residual^2)),
      F = Fm
    )
  )

}

# Is `achieved` on the better-fitting side of `target`? ----
better_side <- function(metric, achieved, target, tol = tolerance)
{
  if(metric %in% c("CFI", "TLI")){
    return(achieved >= target - tol)
  }
  return(achieved <= target + tol)
}

# Target unreachable because the matrix sits on the positive definiteness frontier ----
# (a warning was given, the closest solution is on the better-fitting side of
# the target, and the matrix is numerically singular)
on_frontier <- function(run, metric, target)
{
  if(inherits(run$out, "error") || length(run$warns) == 0){
    return(FALSE)
  }
  eigen_min <- min(eigen(run$out$population$R, symmetric = TRUE, only.values = TRUE)$values)
  return(
    better_side(metric, run$out$population_error$achieved, target) && eigen_min < 1e-3
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Networks to validate ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

cat("Simulating networks...\n")
dm <- matrix(0.20, nrow = 3, ncol = 3); diag(dm) <- 0.90
nets <- list(
  sbm18 = simulate_sbm(nodes = 6, blocks = 3, density_matrix = dm, sample_size = 1000),
  sw20 = simulate_smallworld(nodes = 20, density = 0.20, rewire = 0.10, sample_size = 1000),
  sw10 = simulate_smallworld(nodes = 10, density = 0.40, rewire = 0.20, sample_size = 500)
)

# Degrees of freedom of the known graph ----
get_df <- function(g)
{
  net <- g$population$Omega
  return(sum(net == 0 & upper.tri(net)))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Checks that apply to every method (A) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Independent recomputation of the reported fit ----
independent_checks <- function(out, g, log, label, tag)
{

  net <- g$population$Omega
  df <- get_df(g)
  pe <- out$population_error
  Rp <- out$population$R
  names_fit <- c("CFI", "TLI", "RMSEA", "SRMR")

  # Independent refit and indices
  ind <- ind_indices(Rp, ind_refit(Rp, net), df)

  # Reported fit versus independent
  differences <- abs(pe$fit[names_fit] - ind[names_fit])
  log$add(
    label, paste("fit indices match independent computation", tag),
    all(differences < 1e-5),
    paste("max difference", sci(max(differences)))
  )

  # Reported discrepancy versus independent
  log$add(
    label, paste("delta matches independent discrepancy", tag),
    abs(pe$delta - ind[["F"]]) < 1e-6,
    paste("difference", sci(abs(pe$delta - ind[["F"]])))
  )

  # Reported achieved value is the requested metric
  log$add(
    label, paste("achieved equals reported", pe$metric, tag),
    abs(pe$achieved - pe$fit[[pe$metric]]) < 1e-12,
    ""
  )

  # Population RMSEA identity: RMSEA^2 = F / df - 1 / n
  log$add(
    label, paste("RMSEA^2 = delta / df - 1/n", tag),
    abs(pe$delta / df - 1 / 8.3e09 - pe$fit[["RMSEA"]]^2) < 1e-9,
    ""
  )

}

# Matrix validity ----
validity_checks <- function(out, g, log, label, tag)
{

  R <- g$population$R
  pe <- out$population_error
  Rp <- out$population$R
  E <- pe$error

  log$add(label, paste("symmetric", tag), isSymmetric(unname(Rp), tol = 1e-12), "")
  log$add(
    label, paste("unit diagonal", tag),
    all(abs(diag(Rp) - 1) < 1e-12), ""
  )
  log$add(
    label, paste("positive definite", tag),
    min(eigen(Rp, symmetric = TRUE, only.values = TRUE)$values) > 0, ""
  )
  log$add(
    label, paste("error = new R - original R", tag),
    isTRUE(all.equal(unname(E), unname(Rp - R), tolerance = 1e-12)), ""
  )
  log$add(
    label, paste("original_correlation unchanged", tag),
    isTRUE(all.equal(unname(out$original_correlation), unname(R), tolerance = 1e-12)), ""
  )
  log$add(
    label, paste("error diagonal is zero", tag),
    all(abs(diag(E)) < 1e-12), ""
  )
  log$add(
    label, paste("max_error equals max |error|", tag),
    abs(pe$max_error - max(abs(E))) < 1e-12, ""
  )
  log$add(
    label, paste("data dimensions preserved", tag),
    identical(dim(out$data), dim(g$data)), ""
  )

}

# Common job: one method x metric x network ----
common_job <- function(method, metric, net_name)
{

  g <- nets[[net_name]]
  R <- g$population$R
  net <- g$population$Omega
  target <- cutoffs[[metric]]
  tag <- paste0("[", metric, " = ", target, ", ", net_name, "]")
  log <- new_log()

  # Run at the standard cut-off (cap does not interfere)
  run <- run_pe(g, method = method, metric = metric, fit_value = target, max_residual = 1, max_iterations = 50)
  if(inherits(run$out, "error")){
    log$add(method, paste("runs", tag), FALSE, conditionMessage(run$out))
    return(log$get())
  }
  out <- run$out; pe <- out$population_error

  # 1. Target reached without warning (or provably beyond the PD frontier)
  frontier <- on_frontier(run, metric, target)
  log$add(
    method, paste("target reached without warning", tag),
    (abs(pe$achieved - target) <= tolerance && length(run$warns) == 0) || frontier,
    sprintf(
      "achieved %.5f, warnings %d%s", pe$achieved, length(run$warns),
      if(frontier){" (target beyond positive definiteness frontier)"}else{""}
    )
  )

  # 2-4. Independent computation, validity
  independent_checks(out, g, log, method, tag)
  validity_checks(out, g, log, method, tag)

  # 5. Worse than perfect fit
  log$add(
    method, paste("achieved fit is imperfect", tag),
    if(metric %in% c("CFI", "TLI")){pe$achieved < 1}else{pe$achieved > 0},
    ""
  )

  # 6. `max_residual` cap (unreachable target)
  cap <- 0.02
  capped <- run_pe(g, method = method, metric = metric, fit_value = target, max_residual = cap)
  if(inherits(capped$out, "error")){
    log$add(method, paste("max_residual = 0.02 runs", tag), FALSE, conditionMessage(capped$out))
  }else{
    cp <- capped$out$population_error
    reached <- abs(cp$achieved - target) <= tolerance
    log$add(
      method, paste("max_residual = 0.02 respected", tag),
      cp$max_error <= cap + 1e-12,
      sprintf("max error %.5f", cp$max_error)
    )
    log$add(
      method, paste("unreachable target warns and stays on better side", tag),
      reached || (length(capped$warns) > 0 && better_side(metric, cp$achieved, target)),
      sprintf("achieved %.4f, warnings %d", cp$achieved, length(capped$warns))
    )
  }

  return(log$get())

}

# Stricter targets produce more error (on average) ----
strictness_job <- function(method, metric)
{

  g <- nets$sbm18
  log <- new_log()
  strict <- list(
    CFI = c(0.99, 0.97, 0.95), TLI = c(0.99, 0.97, 0.95),
    RMSEA = c(0.02, 0.04, 0.06), SRMR = c(0.03, 0.05, 0.07)
  )[[metric]]

  # Mean error norm across replications
  norms <- vapply(strict, function(target){
    mean(vapply(seq_len(reps * 3), function(i){
      out <- suppressWarnings(
        add_population_error(g, method = method, metric = metric, fit_value = target, max_residual = 1)
      )
      if(method == "hb"){ # size of the omitted edges (partial correlation)
        P <- -cov2cor(solve(out$population$R)); diag(P) <- 0
        max(abs(P - g$population$Omega))
      }else{ # size of the error added to the correlations
        sqrt(sum(out$population_error$error^2))
      }
    }, numeric(1)))
  }, numeric(1))

  log$add(
    method, paste("worse targets give larger error [", metric, "]"),
    all(diff(norms) > 0),
    paste("sizes at", paste(strict, collapse = "/"), ":", paste(round(norms, 3), collapse = " < "))
  )

  return(log$get())

}

# Data-level end-to-end: large sample ----
data_job <- function(method, metric)
{

  log <- new_log()
  target <- cutoffs[[metric]]
  tag <- paste0("[", metric, " = ", target, "]")
  g <- nets_big
  net <- g$population$Omega
  df <- get_df(g)
  N <- nrow(g$data)

  run <- run_pe(g, method = method, metric = metric, fit_value = target, max_residual = 1, max_iterations = 50)
  if(inherits(run$out, "error")){
    log$add(method, paste("large-sample run", tag), FALSE, conditionMessage(run$out))
    return(log$get())
  }
  if(on_frontier(run, metric, target)){ # target not attainable for this network
    log$add(method, paste("large-sample checks skipped (target beyond PD frontier)", tag), TRUE, "")
    return(log$get())
  }
  out <- run$out; Rp <- out$population$R

  # Sample correlation recovers population matrix
  S <- cor(out$data)
  log$add(
    method, paste("sample correlation matches population R (N = 2e5)", tag),
    max(abs(S - Rp)) < 0.012,
    paste("max difference", sci(max(abs(S - Rp))))
  )

  # Sample fit of true network is close to the target (RMSEA, SRMR)
  sample_fit <- ind_indices(S, ind_refit(S, net), df, n = N)
  key <- if(metric %in% c("RMSEA", "SRMR")){metric}else{NULL}
  if(!is.null(key)){
    log$add(
      method, paste("sample", key, "within 0.01 of target", tag),
      abs(sample_fit[[key]] - target) < 0.01,
      sprintf("sample %s = %.4f", key, sample_fit[[key]])
    )
  }else{ # CFI and TLI: sample value close to the population value
    log$add(
      method, paste("sample", metric, "within 0.01 of target", tag),
      abs(sample_fit[[metric]] - out$population_error$achieved) < 0.01 ||
        abs(sample_fit[[metric]] - target) < 0.01,
      sprintf("sample (N-based) %s = %.4f", metric, sample_fit[[metric]])
    )
  }

  return(log$get())

}

# Input validation ----
input_job <- function()
{

  log <- new_log()
  g <- nets$sbm18

  # Expect an error from each call
  expect_error <- function(label, expr_fun){
    result <- tryCatch({expr_fun(); "no error"}, error = function(e) conditionMessage(e))
    log$add("input", label, !identical(result, "no error") && nchar(result) > 10, substr(result, 1, 70))
  }

  expect_error("bad method", function() add_population_error(g, method = "zzz"))
  expect_error("bad metric", function() add_population_error(g, metric = "aic"))
  expect_error("fit_value = 0", function() add_population_error(g, fit_value = 0))
  expect_error("fit_value = 1", function() add_population_error(g, fit_value = 1))
  expect_error("fit_value negative", function() add_population_error(g, fit_value = -0.1))
  expect_error("fit_value length 2", function() add_population_error(g, fit_value = c(0.9, 0.95)))
  expect_error("fit_value character", function() add_population_error(g, fit_value = "close"))
  expect_error("max_residual = 0", function() add_population_error(g, max_residual = 0))
  expect_error("max_residual negative", function() add_population_error(g, max_residual = -1))
  expect_error("max_residual character", function() add_population_error(g, max_residual = "a"))
  expect_error("unknown control option", function() add_population_error(g, control = list(foo = 1)))
  expect_error("control not a list", function() add_population_error(g, control = 5))
  expect_error("decay out of range", function() add_population_error(g, method = "tkl", control = list(decay = 2)))
  expect_error("non-ggm object", function() add_population_error(list(a = 1)))

  # Population error cannot be added twice
  once <- suppressWarnings(add_population_error(g))
  expect_error("ggm_pe object passed again", function() add_population_error(once))

  # Fully connected network
  full <- g
  full$population$Omega[] <- 0.1; diag(full$population$Omega) <- 0
  expect_error("fully connected network", function() add_population_error(full))

  # Case-insensitive options
  ok <- tryCatch({
    o <- suppressWarnings(add_population_error(g, method = "CB", metric = "rmsea"))
    o$population_error$method == "cb" && o$population_error$metric == "RMSEA"
  }, error = function(e) FALSE)
  log$add("input", "method and metric are case-insensitive", ok, "")

  # Defaults
  d <- suppressWarnings(add_population_error(g))
  log$add(
    "input", "defaults: cb, CFI = 0.95, max_residual = 0.10",
    d$population_error$method == "cb" && d$population_error$metric == "CFI" &&
      d$population_error$fit_value == 0.95 && d$population_error$max_residual == 0.10 &&
      d$population_error$max_error <= 0.10 + 1e-12,
    ""
  )
  for(m in metrics){
    dm_out <- suppressWarnings(add_population_error(g, metric = m))
    log$add(
      "input", paste("default fit_value for", m),
      dm_out$population_error$fit_value == cutoffs[[m]], ""
    )
  }

  return(log$get())

}

# Zero error gives perfect fit (cut-off reference) ----
perfect_job <- function()
{

  log <- new_log()
  for(name in names(nets)){
    g <- nets[[name]]
    idx <- ind_indices(g$population$R, ind_refit(g$population$R, g$population$Omega), get_df(g))
    log$add(
      "all", paste("zero error: CFI = 1, TLI >= 1, RMSEA = 0, SRMR = 0 [", name, "]"),
      abs(idx[["CFI"]] - 1) < 1e-8 && idx[["TLI"]] >= 1 - 1e-8 &&
        idx[["RMSEA"]] < 1e-5 && idx[["SRMR"]] < 1e-6,
      paste(round(idx[c("CFI", "TLI", "RMSEA", "SRMR")], 8), collapse = " ")
    )
  }

  return(log$get())

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Method-specific checks (B): concept vs. behavior ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# ---- CB: Cudeck & Browne (1992) ----
# Error only on cells the network constrains to zero; the true network
# (with the original R) remains the minimizer of the discrepancy
cb_checks <- function(out, g, log, label, tag)
{

  R <- g$population$R
  net <- g$population$Omega
  p <- ncol(R)
  df <- get_df(g)
  pe <- out$population_error
  E <- pe$error; Rp <- out$population$R

  # Error is zero on edges and the diagonal; non-zero on constrained cells
  kept <- net != 0; diag(kept) <- TRUE
  constrained <- net == 0 & upper.tri(net)
  log$add(
    label, paste("error is zero on edge cells and diagonal", tag),
    all(abs(E[kept]) < 1e-12),
    paste("max", sci(max(abs(E[kept]))))
  )
  log$add(
    label, paste("error is non-zero on constrained cells", tag),
    any(abs(E[constrained]) > 1e-8),
    paste(sum(abs(E[constrained]) > 1e-8), "of", sum(constrained), "cells")
  )

  # Orthogonality: known graph refit reproduces original R
  refit <- ind_refit(Rp, net)
  log$add(
    label, paste("known graph refit reproduces original R", tag),
    max(abs(refit - R)) < 1e-6 && pe$R_difference < 1e-6,
    paste("max difference", sci(max(abs(refit - R))))
  )

  # Minimizer: no other network-consistent matrix fits Rp better than R
  F0 <- ind_discrepancy(R, Rp)
  K <- solve(R)
  allowed <- net != 0; diag(allowed) <- TRUE
  alternatives <- vapply(seq_len(200), function(i){
    noise <- matrix(rnorm(p * p, sd = 0.02), p, p)
    noise <- (noise + t(noise)) / 2
    noise[!allowed] <- 0
    Kc <- K + noise
    if(min(eigen(Kc, symmetric = TRUE, only.values = TRUE)$values) <= 1e-6){
      return(NA_real_)
    }
    ind_discrepancy(solve(Kc), Rp)
  }, numeric(1))
  log$add(
    label, paste("original R is the minimizer (200 alternatives)", tag),
    sum(!is.na(alternatives)) >= 100 && all(alternatives >= F0 - 1e-9, na.rm = TRUE),
    sprintf("min alternative - F0 = %s", sci(min(alternatives - F0, na.rm = TRUE)))
  )
  log$add(
    label, paste("delta is the minimum discrepancy", tag),
    abs(pe$delta - F0) < 1e-6, ""
  )

  # Closed-form fit (against R) equals refit fit
  closed <- ind_indices(Rp, R, df)
  refit_fit <- ind_indices(Rp, refit, df)
  log$add(
    label, paste("closed-form fit equals refit fit", tag),
    max(abs(closed[1:4] - refit_fit[1:4])) < 1e-5,
    paste("max difference", sci(max(abs(closed[1:4] - refit_fit[1:4]))))
  )

  # Misfit grows monotonically with the amount of error
  s <- seq(0, 1, by = 0.1)
  path <- vapply(s, function(x) ind_indices(R + x * E, R, df), numeric(5))
  log$add(
    label, paste("fit worsens monotonically as error is scaled 0 to 1", tag),
    abs(path["F", 1]) < 1e-12 && all(diff(path["F", ]) > 0) &&
      all(diff(path["CFI", ]) <= 1e-12) && all(diff(path["TLI", ]) <= 1e-12) &&
      all(diff(path["RMSEA", ]) >= -1e-12) && all(diff(path["SRMR", ]) >= -1e-12),
    ""
  )

}

cb_job <- function(net_name)
{
  g <- nets[[net_name]]
  log <- new_log()
  for(metric in metrics){
    tag <- paste0("[", metric, ", ", net_name, "]")
    run <- run_pe(g, method = "cb", metric = metric, fit_value = cutoffs[[metric]], max_residual = 1)
    if(inherits(run$out, "error")){
      log$add("cb", paste("runs", tag), FALSE, conditionMessage(run$out))
      next
    }
    cb_checks(run$out, g, log, "cb", tag)
  }
  return(log$get())
}

# ---- TKL: Tucker, Koopman & Linn (1969) ----
# R_pop = (1 - w) R + w Q, Q = L L' from k minor factors with geometrically
# decaying loadings (rows of unit length)

# Minor-factor structure: R_pop - (1 - w) R = w Q has rank <= k, PSD, diagonal w
tkl_structure <- function(out, g, k)
{

  R <- g$population$R
  Rp <- out$population$R
  p <- ncol(R)

  # Residual rank beyond k for a given weight
  excess <- function(w){
    values <- eigen(Rp - (1 - w) * R, symmetric = TRUE, only.values = TRUE)$values
    sum(sort(values, decreasing = TRUE)[-seq_len(k)]^2)
  }

  # Grid then refine
  grid <- seq(0.0005, 0.9995, length.out = 400)
  start <- grid[which.min(vapply(grid, excess, numeric(1)))]
  best <- optimize(excess, interval = c(max(start - 0.003, 1e-6), min(start + 0.003, 1 - 1e-6)), tol = 1e-14)
  w <- best$minimum
  M <- Rp - (1 - w) * R
  min_eigen <- min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

  return(
    list(
      w = w, relative_residual = sqrt(best$objective) / w,
      min_eigen = min_eigen, diag_error = max(abs(diag(M) - w))
    )
  )

}

tkl_checks <- function(out, g, log, label, tag)
{

  net <- g$population$Omega
  pe <- out$population_error
  E <- pe$error
  edge_cells <- net != 0

  # Error is not confined to constrained cells (unlike CB)
  log$add(
    label, paste("error also lands on true edges", tag),
    any(abs(E[edge_cells]) > 1e-8) && pe$R_difference > 1e-6,
    paste("R_difference", sci(pe$R_difference))
  )

  # Minor-factor structure for small and large numbers of factors
  for(k in c(1, 2, 5)){
    run <- run_pe(g, method = "tkl", metric = "CFI", fit_value = 0.95, max_residual = 1, control = list(minor_factors = k))
    if(inherits(run$out, "error")){
      log$add(label, paste("minor_factors =", k, "runs", tag), FALSE, conditionMessage(run$out))
      next
    }
    s <- tkl_structure(run$out, g, k)
    log$add(
      label, paste0("R_pop - (1-w)R has rank <= ", k, ", is PSD, diagonal w ", tag),
      s$relative_residual < 1e-5 && s$min_eigen > -1e-8 && s$diag_error < 1e-10,
      sprintf("w = %.4f, rel. residual %s, min eig %s", s$w, sci(s$relative_residual), sci(s$min_eigen))
    )
  }

  # Default minor factors (max(50, 2 * nodes)) give a full-rank minor-factor component
  run <- run_pe(g, method = "tkl", metric = "CFI", fit_value = 0.95, max_residual = 1)
  if(!inherits(run$out, "error")){
    s <- tkl_structure(run$out, g, ncol(g$population$R) - 1)
    log$add(
      label, paste("default factors: component has rank >= p - 1", tag),
      s$min_eigen > -1e-8, sprintf("w = %.4f", s$w)
    )
  }

  # Random effect: two calls differ but both reach the target
  a <- run_pe(g, method = "tkl", metric = "RMSEA", fit_value = 0.06, max_residual = 1)
  b <- run_pe(g, method = "tkl", metric = "RMSEA", fit_value = 0.06, max_residual = 1)
  if(!inherits(a$out, "error") && !inherits(b$out, "error")){
    log$add(
      label, paste("repeated calls differ (random minor factors) but reach target", tag),
      max(abs(a$out$population$R - b$out$population$R)) > 1e-6 &&
        abs(a$out$population_error$achieved - 0.06) <= tolerance &&
        abs(b$out$population_error$achieved - 0.06) <= tolerance,
      ""
    )
  }

}

tkl_job <- function(net_name)
{
  g <- nets[[net_name]]
  log <- new_log()
  tag <- paste0("[", net_name, "]")
  run <- run_pe(g, method = "tkl", metric = "CFI", fit_value = 0.95, max_residual = 1)
  if(inherits(run$out, "error")){
    log$add("tkl", paste("runs", tag), FALSE, conditionMessage(run$out))
    return(log$get())
  }
  tkl_checks(run$out, g, log, "tkl", tag)
  return(log$get())
}

# Spectrum concentration depends on decay: larger decay -> fewer effective factors
tkl_decay_job <- function(decay)
{

  g <- nets$sbm18
  R <- g$population$R; net <- g$population$Omega
  p <- ncol(R); df <- get_df(g)
  TKL_optimize <- internal("TKL_optimize")

  shares <- vapply(seq_len(40), function(i){
    sol <- TKL_optimize("CFI", 0.95, R, net, p, df, 1, 50, decay)
    if(!isTRUE(sol$valid) || sol$parameter <= 0){
      return(NA_real_)
    }
    Q <- (sol$R_population - (1 - sol$parameter) * R) / sol$parameter
    eigen(Q, symmetric = TRUE, only.values = TRUE)$values[1] / p
  }, numeric(1))

  data.frame(
    method = "tkl", check = "decay", pass = TRUE,
    detail = paste(decay, mean(shares, na.rm = TRUE), sd(shares, na.rm = TRUE), sum(!is.na(shares))),
    stringsAsFactors = FALSE
  )

}

# ---- WB: Wu & Browne (2015) ----
# Sigma ~ inverted Wishart centered on R, dispersion v = 1/m ~ RMSEA^2
wb_checks <- function(out, g, log, label, tag)
{

  net <- g$population$Omega
  pe <- out$population_error
  E <- pe$error
  edge_cells <- net != 0
  non_edge <- net == 0 & upper.tri(net)

  # Error is spread across all cells (edges and non-edges)
  log$add(
    label, paste("error lands on both true edges and constrained cells", tag),
    any(abs(E[edge_cells]) > 1e-8) && any(abs(E[non_edge]) > 1e-8) && pe$R_difference > 1e-6,
    paste("R_difference", sci(pe$R_difference))
  )

  # Random effect: two calls differ but both reach the target
  a <- run_pe(g, method = "wb", metric = "RMSEA", fit_value = 0.06, max_residual = 1)
  b <- run_pe(g, method = "wb", metric = "RMSEA", fit_value = 0.06, max_residual = 1)
  if(!inherits(a$out, "error") && !inherits(b$out, "error")){
    log$add(
      label, paste("repeated calls differ (random adventitious error) but reach target", tag),
      max(abs(a$out$population$R - b$out$population$R)) > 1e-6 &&
        abs(a$out$population_error$achieved - 0.06) <= tolerance &&
        abs(b$out$population_error$achieved - 0.06) <= tolerance,
      ""
    )
  }

}

wb_job <- function(net_name)
{
  g <- nets[[net_name]]
  log <- new_log()
  tag <- paste0("[", net_name, "]")
  run <- run_pe(g, method = "wb", metric = "CFI", fit_value = 0.95, max_residual = 1)
  if(inherits(run$out, "error")){
    log$add("wb", paste("runs", tag), FALSE, conditionMessage(run$out))
    return(log$get())
  }
  wb_checks(run$out, g, log, "wb", tag)
  return(log$get())
}

# Monte Carlo of the inverted-Wishart draw at fixed dispersion v
wb_monte_carlo <- function(g, v, n_draws)
{

  R <- g$population$R; net <- g$population$Omega
  p <- ncol(R); df <- get_df(g)
  WB_error <- internal("WB_error")
  rnorm_ziggurat <- internal("rnorm_ziggurat")
  runif_xoshiro <- internal("runif_xoshiro")

  draw_one <- function(i){
    normals <- matrix(0, nrow = p, ncol = p)
    normals[lower.tri(normals)] <- rnorm_ziggurat(p * (p - 1) / 2)
    d <- list(
      normals = normals, uniforms = runif_xoshiro(p),
      cholesky = t(chol(chol2inv(chol(R))))
    )
    fit <- WB_error(log(v), d, R, net, df, Inf)
    if(!isTRUE(fit$valid)){
      return(NULL)
    }
    list(error = (fit$R_population - R)[upper.tri(R)], RMSEA = fit$RMSEA)
  }

  results <- mclapply(seq_len(n_draws), draw_one, mc.cores = cores)
  results <- results[!vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))]

  list(
    errors = do.call(rbind, lapply(results, `[[`, "error")),
    RMSEA = vapply(results, `[[`, numeric(1), "RMSEA")
  )

}

wb_monte_carlo_checks <- function(net_name, log)
{

  g <- nets[[net_name]]
  tag <- paste0("[", net_name, "]")
  v_values <- c(0.0004, 0.0036, 0.01)
  variances <- numeric(length(v_values))

  for(i in seq_along(v_values)){

    v <- v_values[i]
    mc <- wb_monte_carlo(g, v, draws)
    errors <- mc$errors

    # Centering on the original R
    z <- colMeans(errors) / (apply(errors, 2, sd) / sqrt(nrow(errors)))
    log$add(
      "wb", paste0("draws centered on original R (v = ", v, ", ", nrow(errors), " draws) ", tag),
      max(abs(z)) < 5,
      sprintf("max |z| = %.2f", max(abs(z)))
    )

    # Paper (Eq. 16): v ~ RMSEA^2
    ratio <- mean(mc$RMSEA) / sqrt(v)
    log$add(
      "wb", paste0("mean RMSEA ~ sqrt(v) (v = ", v, ") ", tag),
      abs(ratio - 1) < 0.10,
      sprintf("mean RMSEA %.4f vs sqrt(v) %.4f (ratio %.3f)", mean(mc$RMSEA), sqrt(v), ratio)
    )

    variances[i] <- mean(errors^2)

  }

  # Error variance scales proportionally with dispersion
  slope <- unname(coef(lm(log(variances) ~ log(v_values)))[2])
  log$add(
    "wb", paste("error variance grows proportionally with v", tag),
    slope > 0.9 && slope < 1.15,
    sprintf("log-log slope = %.3f", slope)
  )

}

# ---- HB: Hu & Bentler (1998) ----
# Underparameterized misspecification: omitted edges in the true network
hb_checks <- function(out, g, log, label, tag, omitted = 0.10)
{

  R <- g$population$R
  net <- g$population$Omega
  df <- get_df(g)
  pe <- out$population_error
  Rp <- out$population$R

  # Partial correlations of a correlation matrix
  pcor <- function(M){
    P <- -cov2cor(solve(M)); diag(P) <- 0
    P
  }

  # Precondition: original R implies the network
  log$add(
    label, paste("original R implies Omega", tag),
    max(abs(pcor(R) - net)) < 1e-6, ""
  )

  # Changes in partial correlations
  D <- pcor(Rp) - net
  changed <- which(abs(D) > 1e-9 & upper.tri(D))
  expected <- max(1, round(omitted * df))
  log$add(
    label, paste("number of omitted edges =", expected, tag),
    length(changed) == expected,
    paste("found", length(changed))
  )
  log$add(
    label, paste("omitted edges are zero edges of the true network", tag),
    length(changed) > 0 && all(net[changed] == 0), ""
  )
  sizes <- abs(D[changed])
  log$add(
    label, paste("omitted edges have equal magnitude", tag),
    length(sizes) > 0 && (max(sizes) - min(sizes)) < 1e-8,
    sprintf("size %.4f", mean(sizes))
  )

  # Data-generating network (with omitted edges) fits perfectly, true network does not
  added <- matrix(0, nrow = nrow(net), ncol = ncol(net))
  added[changed] <- D[changed]
  full <- net + added + t(added)
  log$add(
    label, paste("network including omitted edges fits exactly; true network does not", tag),
    {
      Sigma_full <- ind_refit(Rp, full)
      ind_discrepancy(Sigma_full, Rp) < 1e-8 && pe$delta > 1e-4
    },
    sprintf("delta %.4f", pe$delta)
  )

  # Error propagates beyond omitted cells
  log$add(
    label, paste("error propagates beyond the omitted cells", tag),
    pe$R_difference > 1e-6, paste("R_difference", sci(pe$R_difference))
  )

}

hb_job <- function(net_name)
{
  g <- nets[[net_name]]
  log <- new_log()
  tag <- paste0("[", net_name, "]")
  run <- run_pe(g, method = "hb", metric = "CFI", fit_value = 0.95, max_residual = 1)
  if(inherits(run$out, "error")){
    log$add("hb", paste("runs", tag), FALSE, conditionMessage(run$out))
    return(log$get())
  }
  hb_checks(run$out, g, log, "hb", tag)
  return(log$get())
}

# The `omitted` control: more omitted edges, smaller edges at the same fit
hb_control_job <- function()
{

  g <- nets$sbm18
  net <- g$population$Omega
  df <- get_df(g)
  log <- new_log()
  pcor <- function(M){P <- -cov2cor(solve(M)); diag(P) <- 0; P}
  levels <- c(0.05, 0.10, 0.25)

  summary_rows <- lapply(levels, function(omitted){
    runs <- lapply(seq_len(reps), function(i){
      out <- suppressWarnings(
        add_population_error(
          g, method = "hb", metric = "CFI", fit_value = 0.95,
          max_residual = 1, control = list(omitted = omitted)
        )
      )
      D <- pcor(out$population$R) - net
      changed <- which(abs(D) > 1e-9 & upper.tri(D))
      c(k = length(changed), size = mean(abs(D[changed])), achieved = out$population_error$achieved)
    })
    do.call(rbind, runs)
  })

  expected <- pmax(1, round(levels * df))
  log$add(
    "hb", "control$omitted sets the number of omitted edges",
    all(vapply(seq_along(levels), function(i) all(summary_rows[[i]][, "k"] == expected[i]), logical(1))),
    paste("expected", paste(expected, collapse = "/"))
  )
  mean_sizes <- vapply(summary_rows, function(m) mean(m[, "size"]), numeric(1))
  log$add(
    "hb", "more omitted edges -> smaller per-edge size at the same fit",
    all(diff(mean_sizes) < 0),
    paste("mean sizes", paste(round(mean_sizes, 4), collapse = " > "))
  )

  # Both signs of omitted edges occur
  out <- suppressWarnings(
    add_population_error(g, method = "hb", metric = "CFI", fit_value = 0.95, max_residual = 1, control = list(omitted = 0.5))
  )
  D <- pcor(out$population$R) - net
  changed <- which(abs(D) > 1e-9 & upper.tri(D))
  log$add(
    "hb", "omitted edges have random signs (omitted = 0.5)",
    any(D[changed] > 0) && any(D[changed] < 0),
    paste(sum(D[changed] > 0), "positive,", sum(D[changed] < 0), "negative")
  )

  return(log$get())

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Run the validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

results <- list()
section <- function(title)
{
  cat(sprintf("\n[%s] %s\n", format(round(difftime(Sys.time(), start_time, units = "secs"))), title))
}

# A. Common checks
section("A. Common contract (4 methods x 4 metrics x 3 networks)")
specs <- list()
for(m in methods) for(k in metrics) for(n in names(nets)){
  specs[[length(specs) + 1]] <- list(fun = common_job, args = list(m, k, n), label = m)
}
results$common <- run_jobs(specs)

section("A. Stricter targets give more error")
specs <- list()
for(m in methods) for(k in metrics){
  specs[[length(specs) + 1]] <- list(fun = strictness_job, args = list(m, k), label = m)
}
results$strict <- run_jobs(specs)

section("A. Data-level checks (N = 2e5)")
nets_big <- simulate_sbm(nodes = 6, blocks = 3, density_matrix = dm, sample_size = 2e5)
specs <- list()
for(m in methods) for(k in metrics){
  specs[[length(specs) + 1]] <- list(fun = data_job, args = list(m, k), label = m)
}
results$data <- run_jobs(specs)

section("A. Input validation and cut-off reference")
results$input <- run_jobs(list(
  list(fun = input_job, args = list(), label = "input"),
  list(fun = perfect_job, args = list(), label = "all")
))

# B. Method-specific
section("B1. CB: error on constrained cells only, original R is the minimizer")
results$cb <- run_jobs(lapply(names(nets), function(n) list(fun = cb_job, args = list(n), label = "cb")))

section("B2. TKL: minor-factor structure")
results$tkl <- run_jobs(lapply(names(nets), function(n) list(fun = tkl_job, args = list(n), label = "tkl")))

decay_rows <- run_jobs(list(
  list(fun = tkl_decay_job, args = list(0.02), label = "tkl"),
  list(fun = tkl_decay_job, args = list(0.50), label = "tkl")
))
decay_stats <- lapply(strsplit(decay_rows$detail, " "), as.numeric)
low <- decay_stats[[1]]; high <- decay_stats[[2]]
decay_log <- new_log()
decay_se <- sqrt(low[3]^2 / low[4] + high[3]^2 / high[4])
decay_log$add(
  "tkl", "larger decay concentrates the minor-factor spectrum (first eigenvalue share)",
  high[2] - low[2] > 3 * decay_se,
  sprintf("decay 0.02: %.3f, decay 0.50: %.3f (SE %.3f)", low[2], high[2], decay_se)
)
results$tkl_decay <- decay_log$get()

section("B3. WB: inverted-Wishart adventitious error")
results$wb <- run_jobs(lapply(names(nets), function(n) list(fun = wb_job, args = list(n), label = "wb")))
wb_log <- new_log()
for(n in c("sbm18", "sw20")){
  wb_monte_carlo_checks(n, wb_log)
}
results$wb_mc <- wb_log$get()

section("B4. HB: omitted edges (underparameterized misspecification)")
results$hb <- run_jobs(lapply(names(nets), function(n) list(fun = hb_job, args = list(n), label = "hb")))
results$hb_control <- run_jobs(list(list(fun = hb_control_job, args = list(), label = "hb")))

# C. Self-tests: the checks must fail when applied to the wrong method
section("C. Self-tests (checks must detect wrong behavior)")
g <- nets$sbm18
self <- new_log()
count_failures <- function(fun){
  inner <- new_log()
  fun(inner)
  rows <- inner$get()
  sum(!rows$pass)
}
tkl_out <- suppressWarnings(add_population_error(g, method = "tkl"))
cb_out <- suppressWarnings(add_population_error(g, method = "cb"))
wb_out <- suppressWarnings(add_population_error(g, method = "wb"))
self$add(
  "self-test", "CB checks fail on TKL output",
  count_failures(function(l) cb_checks(tkl_out, g, l, "x", "")) >= 1, ""
)
self$add(
  "self-test", "HB checks fail on CB output",
  count_failures(function(l) hb_checks(cb_out, g, l, "x", "")) >= 1, ""
)
self$add(
  "self-test", "TKL rank-1 structure check fails on WB output",
  {s <- tkl_structure(wb_out, g, 1); s$relative_residual > 1e-3 || s$min_eigen < -1e-6}, ""
)
self$add(
  "self-test", "WB checks fail on CB output",
  count_failures(function(l) wb_checks(cb_out, g, l, "x", "")) >= 1, ""
)
tampered <- cb_out
tampered$population_error$fit[["RMSEA"]] <- tampered$population_error$fit[["RMSEA"]] * 1.1
self$add(
  "self-test", "independent computation detects a tampered fit index",
  count_failures(function(l) independent_checks(tampered, g, l, "x", "")) >= 1, ""
)
results$self <- self$get()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Report ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

all_results <- do.call(rbind, results)
all_results$method <- factor(all_results$method, levels = c(methods, "input", "all", "self-test"))

# Optionally save all results (VALIDATE_CSV = file path)
if(nzchar(Sys.getenv("VALIDATE_CSV"))){
  write.csv(all_results, Sys.getenv("VALIDATE_CSV"), row.names = FALSE)
}

cat("\n\n==================== Summary ====================\n")
summary_table <- do.call(rbind, lapply(split(all_results, all_results$method), function(d){
  data.frame(checks = nrow(d), passed = sum(d$pass), failed = sum(!d$pass))
}))
print(summary_table)

failures <- all_results[!all_results$pass, ]
if(nrow(failures) > 0){
  cat("\n==================== Failures ====================\n")
  for(i in seq_len(nrow(failures))){
    cat(sprintf("[%s] %s\n    %s\n", failures$method[i], failures$check[i], failures$detail[i]))
  }
}

cat(sprintf(
  "\n%d checks, %d passed, %d failed (%s, %d cores%s)\n",
  nrow(all_results), sum(all_results$pass), nrow(failures),
  format(round(difftime(Sys.time(), start_time, units = "mins"), 1)), cores,
  if(quick){", quick mode"}else{""}
))

if(nrow(failures) > 0){
  quit(status = 1)
}
