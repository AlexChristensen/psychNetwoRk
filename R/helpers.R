#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### {psychNetwoRk} Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Numeric vector/matrix output apply ----
# Updated 06.07.2023
nvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = numeric(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Unlist `lapply` ----
# Updated 14.07.2023
ulapply <- function(X, FUN, ..., recursive = TRUE)
{
  return(unlist(lapply(X, FUN, ...), recursive = recursive))
}

#' @noRd
# Reproducible random number generation ----
# Exported functions with a `seed` argument open a seed scope: while it is
# active, every call to the C random number generators (without its own
# seed) uses a deterministic seed derived from the scope's seed and a call
# counter. All random numbers in the package come from the C generators; R's
# random number generator is also set with `set.seed(seed)` so that functions
# from other packages used within (e.g., {igraph}'s random graphs and Leiden
# algorithm) are reproducible as well. Nested calls (e.g., `simulate_sbm`
# calling `smallworldness`) share the outer scope. When the scope ends, R's
# random number generator state is restored, so a seeded call does not change
# the user's random stream
# Updated 20.11.2026
pn_rng <- new.env(parent = emptyenv())

#' @noRd
# Start a seed scope (returns a token for `seed_end`; NULL when no scope
# was started)
# Updated 20.11.2026
seed_start <- function(seed)
{

  # No seed or a scope is already active
  if(is.null(seed) || isTRUE(pn_rng$active)){
    return(NULL)
  }

  # Check seed
  if(
    !is.numeric(seed) || length(seed) != 1 || !is.finite(seed) ||
    seed != round(seed) || abs(seed) > .Machine$integer.max
  ){
    stop("`seed` must be a single whole number (or NULL for a random seed)", call. = FALSE)
  }

  # Store user's random number generator state
  global <- globalenv()
  token <- list(
    exists = exists(".Random.seed", envir = global, inherits = FALSE)
  )
  if(token$exists){
    token$state <- get(".Random.seed", envir = global, inherits = FALSE)
  }

  # Start scope
  pn_rng$active <- TRUE
  pn_rng$base <- abs(seed) + 1
  pn_rng$counter <- 0

  # Set R's random number generator (for functions from other packages)
  set.seed(seed)

  # Return token
  return(token)

}

#' @noRd
# End a seed scope (restores the user's random number generator state)
# Updated 20.11.2026
seed_end <- function(token)
{

  # Check for scope
  if(is.null(token)){
    return(invisible(NULL))
  }

  # End scope
  pn_rng$active <- FALSE

  # Restore user's random number generator state
  global <- globalenv()
  if(isTRUE(token$exists)){
    assign(".Random.seed", token$state, envir = global)
  }else if(exists(".Random.seed", envir = global, inherits = FALSE)){
    rm(".Random.seed", envir = global)
  }

  # Return invisibly
  return(invisible(NULL))

}

#' @noRd
# Seed for a call to the C random number generators
# (0 = random seed from the clock, process ID, and a counter)
# Updated 20.11.2026
call_seed <- function(seed)
{

  # Supplied seed
  if(!is.null(seed)){
    return(as.numeric(seed))
  }

  # No active scope
  if(!isTRUE(pn_rng$active)){
    return(0)
  }

  # Deterministic seed from the scope's seed and the call counter
  # (non-zero and exactly representable as a double)
  pn_rng$counter <- pn_rng$counter + 1
  return(pn_rng$base * 2097152 + pn_rng$counter)

}

#' @noRd
# Generate uniform data ----
# Allows adjustment of range (values are in the open interval (min, max))
# Updated 20.11.2026
runif_xoshiro <- function(n, min = 0, max = 1, seed = NULL)
{

  # Get values
  values <- .Call(
    "r_xoshiro_uniform",
    as.integer(n),
    call_seed(seed),
    PACKAGE = "psychNetwoRk"
  )

  # Check for changes to minimum and maximum
  if(min != 0 || max != 1){ # transform
    values <- min + (max - min) * values
  }

  # Return call from C
  return(values)

}

#' @noRd
# Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 20.11.2026
shuffle <- function(x, size = length(x), seed = NULL)
{

  # Check size
  if(size > length(x)){
    stop("Cannot take a sample larger than the population without replacement", call. = FALSE)
  }

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle",
      as.integer(seq_along(x)),
      call_seed(seed),
      PACKAGE = "psychNetwoRk"
    )][seq_len(size)]
  )

}

#' @noRd
# Weighted Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 20.11.2026
weighted_shuffle <- function(x, size = length(x), prob, seed = NULL)
{

  # Check size
  if(size > length(x)){
    stop("Cannot take a sample larger than the population without replacement", call. = FALSE)
  }

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_weighted_shuffle",
      as.integer(seq_along(x)), as.numeric(prob),
      call_seed(seed),
      PACKAGE = "psychNetwoRk"
    )][seq_len(size)]
  )

}

#' @noRd
# Shuffle (with replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 20.11.2026
shuffle_replace <- function(x, size = length(x), seed = NULL)
{

  # No draws
  if(size == 0){
    return(x[0])
  }

  # Check for values to sample from
  if(length(x) == 0){
    stop("Cannot sample from an empty vector", call. = FALSE)
  }

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle_replace",
      x, as.integer(size), call_seed(seed),
      PACKAGE = "psychNetwoRk"
    )]
  )

}

#' @noRd
# Random normal generation with Ziggurat ----
# https://people.sc.fsu.edu/~jburkardt/cpp_src/ziggurat/ziggurat.html
# Updated 20.11.2026
rnorm_ziggurat <- function(n, seed = NULL)
{

  # Return call from C
  return(
    .Call(
      "r_ziggurat",
      as.integer(n),
      call_seed(seed),
      PACKAGE = "psychNetwoRk"
    )
  )

}

#' @noRd
# Rank with ties broken at random (C random number generator) ----
# Updated 20.11.2026
rank_random <- function(x)
{

  # Order by value, then by a random tie-breaker
  ranks <- integer(length(x))
  ranks[order(x, runif_xoshiro(length(x)))] <- seq_along(x)

  # Return ranks
  return(ranks)

}

#' @noRd
# Node names ("V1", ..., or zero-padded "V01", ...) ----
# Updated 20.11.2026
node_names <- function(nodes)
{
  return(paste0("V", format_integer(seq_len(nodes), digits(nodes) - 1)))
}

#' @noRd
# Ring distances between nodes ----
# min(|i - j|, n - |i - j|)
# Updated 20.11.2026
ring_distance <- function(nodes)
{

  # Node sequence
  node_sequence <- seq_len(nodes)

  # Distances
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))

  # Return ring distances
  return(pmin(distance_matrix, nodes - distance_matrix))

}

#' @noRd
# Sample skewness ----
# Same as `psych::skew` (type = 3): b1 = m3 / m2^(3/2) * ((n - 1) / n)^(3/2)
# Updated 20.11.2026
sample_skewness <- function(x)
{

  # Center data
  n <- length(x)
  centered <- x - mean(x)

  # Moments
  m2 <- sum(centered^2) / n
  m3 <- sum(centered^3) / n

  # Return skewness
  return((m3 / m2^1.5) * ((n - 1) / n)^1.5)

}

#' @noRd
# Adds skew to a single variable ----
# Updated 30.11.2022
skew_single_variable <- function(data, skew_values){

  # Categorize biased data with updated thresholds
  for(i in (length(skew_values) + 1):1){

    # First category
    if(i == 1){
      data[data < skew_values[i]] <- i
    }else if(i == length(skew_values) + 1){ # Last category
      data[data >= skew_values[i-1]] <- i
    }else{ # Middle category
      data[data >= skew_values[i-1] & data < skew_values[i]] <- i
    }

  }

  # Return skewed data
  return(data)

}

# Based on
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Cache of fitted skew transformations ----
# Updated 20.11.2026
skew_cache <- new.env(parent = emptyenv())

#' @noRd
# Fit the skew transformation in the population ----
# Finds the sinh-arcsinh transformation sinh(k * (asinh(x) + m)) whose skew
# on the standard normal distribution (evaluated on a fine grid of normal
# quantiles) equals the target, along with the population mean and standard
# deviation of the transformed variable (for standardization). Results are
# cached by skew value
# Updated 20.11.2026
skew_transformation <- function(absolute, tolerance = 1e-06, max_iterations = 1000)
{

  # Check cache
  key <- format(absolute, digits = 15)
  if(!is.null(skew_cache[[key]])){
    return(skew_cache[[key]])
  }

  # Standard normal distribution (quantile grid)
  arcsinh_reference <- asinh(qnorm(ppoints(20000)))

  # Optimize location for a given tail weight
  objective <- function(kurtosis){
    optimize(
      f = function(x){
        abs(absolute - sample_skewness(sinh(kurtosis * (arcsinh_reference + x))))
      }, interval = c(0, absolute + 2), tol = tolerance
    )
  }

  # Increase tail weight until the target skew is reached
  kurtosis <- 1
  best <- NULL
  for(iteration in seq_len(max_iterations)){

    # Get result
    result <- objective(kurtosis)

    # Store best result
    if(is.null(best) || result$objective < best$objective){
      best <- list(kurtosis = kurtosis, location = result$minimum, objective = result$objective)
    }

    # Check for tolerance
    if(result$objective < tolerance){
      break
    }

    # Update kurtosis
    kurtosis <- kurtosis + result$objective * 0.25 # learning rate

  }

  # Population moments of the transformed variable
  transformed <- sinh(best$kurtosis * (arcsinh_reference + best$location))
  best$mean <- mean(transformed)
  best$sd <- sqrt(mean((transformed - best$mean)^2))

  # Store and return
  skew_cache[[key]] <- best
  return(best)

}

#' @noRd
# Generates skewed data for continuous data ----
# A sinh-arcsinh transformation whose skew on the standard normal
# distribution equals `skew` (the skew is a population value; sample skews
# vary around it). Negative skew transforms the reflected data and reflects
# back, -f(-x), so the transformation is increasing and correlations keep
# their sign. Skewed variables are standardized with the population mean and
# standard deviation of the transformation (unskewed variables are already
# standard normal in the population)
# Updated 20.11.2026
skew_continuous <- function(skew, data = NULL, sample_size = 1e06)
{

  # Generate data
  if(is.null(data)){
    data <- rnorm_ziggurat(sample_size)
  }

  # Check for zero skew
  if(skew == 0){
    return(data)
  }

  # Obtain sign and transformation
  skew_sign <- sign(skew)
  transformation <- skew_transformation(abs(skew))

  # Transform (reflected for negative skew)
  transformed <- sinh(
    transformation$kurtosis * (asinh(skew_sign * data) + transformation$location)
  )

  # Return standardized result (reflected back for negative skew)
  return(skew_sign * (transformed - transformation$mean) / transformation$sd)

}

#' @noRd
# Faster `ifelse` ----
# For single value replacements:
# 1.5x faster with 1 value
# 2.5x faster with 10 values
# >= 18x faster with >= 100 values
# Updated 26.02.2026
swiftelse <- function(condition, true, false)
{
  # Get condition length
  condition_length <- length(condition)

  # Check for single value
  if(condition_length == 1){

    # If TRUE
    if(condition){
      return(true)
    }else{ # Otherwise, FALSE
      return(false)
    }

  }

  # Initialize result
  result <- vector(typeof(false), condition_length)

  # Set TRUE condition
  if(length(true) == 1){
    result[condition] <- true
  }else{
    result[condition] <- true[condition]
  }

  # Set FALSE condition
  if(length(false) == 1){
    result[!condition] <- false
  }else{

    # Get opposite condition (slightly faster than repeated calls)
    opposite_condition <- !condition
    result[opposite_condition] <- false[opposite_condition]
  }


  # Return result
  return(result)

}

#' @noRd
# Determine number of digits in a number ----
# Updated 24.07.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(abs(number))) + 1)

}

#' @noRd
# Format number with certain integer ----
# Mainly for naming and printing
# Updated 14.06.2024
format_integer <- function(numbers, places)
{

  return(
    formatC(
      x = numbers,
      digits = places,
      format = "d", flag = "0"
    )
  )

}

#' @noRd
# Determine object type ----
# Updated 05.07.2023
get_object_type <- function(object)
{

  # Remove attributes and class from object
  unclassed_object <- remove_attributes(object)

  # Add object type to object's attributes
  ## In order of precedence
  if(is(object, "tbl")){ # needs precedence over 'data.frame'
    return("tibble")
  }else if(is.data.frame(object)){ # needs precedence over 'list'
    return("data.frame")
  }else if(is.list(unclassed_object)){ # needs precedence over 'vector'
    return("list")
  }else if(is.factor(object)){
    return("factor")
  }else if(is.vector(unclassed_object)){
    return("vector")
  }else if(is.matrix(unclassed_object)){
    return("matrix")
  }else if(is.array(unclassed_object)){
    return("array")
  }

}

#' @noRd
# Convert matrix to {igraph} network ----
# Updated 09.08.2023
convert2igraph <- function (A, diagonal = 0)
{

  # Convert to matrix
  A <- as.matrix(A)

  # Change diagonal (to zero)
  diag(A) <- diagonal

  # Return {igraph} network
  return(
    silent_call(
      igraph::graph_from_adjacency_matrix(
        A, weighted = TRUE, mode = "undirected",
        add.colnames = FALSE
      )
    )
  )

}

#' @noRd
# Convert {igraph} network to matrix ----
# Updated 09.08.2023
igraph2matrix <- function(I)
{
  return(as.matrix(igraph::as_adjacency_matrix(I)))
}

#' @noRd
# Faster table ----
# Skips many checks in `table`
# Updated 22.07.2023
fast_table <- function(values)
{

  # Get factor values
  factor_values <- swiftelse(
    is.factor(values), values, factor(values)
  )

  # Get factor levels
  factor_levels <- levels(factor_values)

  # Return with names
  return(
    setNames(
      tabulate(factor_values, length(factor_levels)),
      factor_levels
    )
  )

}

#' @noRd
# All-purpose symmetric checker ----
# Faster but less robust than `isSymmetric`
# Defaults to floating point tolerance
# Updated 01.07.2023
is_symmetric <- function(data, tolerance = 1e-06){

  # Get dimensions
  dimensions <- dim(data)

  # Check for whether rows equal columns
  if(dimensions[1] != dimensions[2]){ # Not a (square) matrix
    return(FALSE)
  }else{

    # Convert to matrix
    data_matrix <- as.matrix(data)

    # Check that all are equal
    return(all(abs(data_matrix - t(data_matrix)) < tolerance))

  }

}

#' @noRd
# General function to silently obtain output ----
# Updated 01.07.2023
silent_call <- function(...)
{

  # Make call
  sink <- capture.output(
    result <- suppressWarnings(
      suppressMessages(
        ...
      )
    )
  )

  # Return result
  return(result)

}

#' @noRd
# Strip classes and attributes
# Updated 09.07.2023
remove_attributes <- function(object)
{

  # Remove all except special attributes (see `?structure`)
  attributes(object) <- attributes(object)[
    names(attributes(object)) %in% c(
      "dim", "dimnames", "names", "tsp", "levels"
    )
  ]

  # Return unclassed object
  return(unclass(object))

}

#' @noRd
# Set default argument (cleaner missing function) ----
# Updated 28.07.2023
set_default <- function(argument, default, FUN, several.ok = FALSE)
{

  # Check for type error
  typeof_error(argument, c(typeof(default), "closure"), "set_default")

  # Return argument if it's a function
  if(is.function(argument)){
    return(argument)
  }

  # Return default if NULL
  if(is.null(argument)){
    return(default)
  }

  # Get argument name
  argument_name <- as.character(substitute(argument))

  # Several values are only allowed when `several.ok = TRUE`
  # (the full vector of choices is the missing-argument default)
  if(length(argument) > 1 && !several.ok){

    # Get choices
    choices <- swiftelse(
      is.function(FUN), tolower(eval(formals(FUN)[[argument_name]])), tolower(FUN)
    )

    # Return default when all choices are given (argument is missing)
    if(setequal(tolower(argument), choices)){
      return(default)
    }

    # Otherwise, error
    stop(
      paste0(
        "Only one value can be used for '", argument_name, "'. Please choose from: ",
        paste0("\"", choices, "\"", collapse = ", ")
      ), call. = FALSE
    )

  }

  # Get choices for the function
  if(!is.function(FUN)){ # choices are provided in 'FUN'
    choices <- FUN
  }else{ # get formal arguments from the function itself

    # Get formal argument choices from default call
    choices <- formals(FUN)[[argument_name]]

    # Check if choices is a call
    if(is.call(choices)){

      # Get choices
      choices <- as.character(choices)

      # Remove "c" from call
      choices <- choices[-which(choices == "c")]

    }

  }

  # Make argument lowercase (if necessary)
  argument <- tolower(argument); choices <- tolower(choices)

  # Get arguments not in choices
  if(any(!argument %in% choices)){ # Repeated calls OK since returning error

    # Get arguments not in choices
    not_in_choices <- !argument %in% choices

    # Get number not in choices
    number_not <- sum(not_in_choices)

    # Send error for bad argument
    if(number_not > 0){

      # Only one not in choices
      if(number_not == 1){
        argument_text <- paste0(" = \"", argument[not_in_choices], "\"")
      }else{ # Multiple not in choices
        argument_text <- paste0(
          " = c(",
          paste0("\"", argument[not_in_choices], "\"", collapse = ", "),
          ")"
        )
      }

      # Throw error
      stop(
        paste0(
          paste0("Invalid argument: ", argument_name, argument_text),
          "\n\nPlease choose from: ",
          paste0("\"", choices, "\"", collapse = ", ")
        ), call. = FALSE
      )

    }

  }

  # Return argument
  return(argument)

}

#' @noRd
# Error for object type ----
# Updated 04.09.2023
object_error <- function(input, expected_type, function_name){

  # Get input type
  input_type <- get_object_type(input)

  # Check for object types
  if(!input_type %in% expected_type){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is a ", paste0("'", input_type, "'", collapse = ", "), " object",
        ". Input is expected to be ", paste0("'", expected_type, "'", collapse = ", "), " object",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#object-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for `typeof` ----
# Updated 04.09.2023
typeof_error <- function(input, expected_value, function_name){

  # Switch out "closure" with "function"
  if("closure" %in% expected_value){
    expected_value[expected_value == "closure"] <- "function"
  }

  # Get type of input
  typeof_input <- typeof(input)

  # Convert "closure" to "function"
  if("closure" %in% typeof_input){
    typeof_input[typeof_input == "closure"] <- "function"
  }

  # Convert "integer" and "double" to "numeric"
  ## Input
  typeof_input <- swiftelse(
    typeof_input %in% c("integer", "double"),
    "numeric", typeof_input
  )
  ## Expected value (keeps other expected types)
  expected_value[expected_value %in% c("integer", "double")] <- "numeric"
  expected_value <- unique(expected_value)

  # Check for value
  if(!typeof_input %in% expected_value){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is ", paste("'", typeof_input, "'", sep = "", collapse = ", "), " type",
        ". Input is expected to be ",
        paste0("'", expected_value, "'", collapse = " or "), " type",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#typeof-error"
        # can use "or" because `typeof` won't ever be more than two
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for `length` ----
# Updated 04.09.2023
length_error <- function(input, expected_lengths, function_name){

  # Check for length of input in expected length
  if(!length(input) %in% expected_lengths){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Length of '", deparse(substitute(input)),
        "' (", length(input),") does not match expected length(s). Length must be: ",
        paste0("'", expected_lengths, "'", collapse = " or "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#length-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Error for `range` ----
# Updated 04.09.2023
range_error <- function(input, expected_ranges, function_name){

  # Obtain expected maximum and minimum values
  expected_maximum <- max(expected_ranges, na.rm = TRUE)
  expected_minimum <- min(expected_ranges, na.rm = TRUE)

  # Check for missing values
  if(anyNA(input)){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)), "' contains missing (NA) values"
      ),
      call = function_name
    )
  }

  # Obtain maximum and minimum values
  # (compared unrounded with a small tolerance for floating point error;
  # rounded only for the message)
  actual_maximum <- max(input)
  actual_minimum <- min(input)
  tolerance <- 1e-12

  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum + tolerance){
    actual_maximum <- signif(actual_maximum, 6)
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Maximum value of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Values must range between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#range-error"
      ),
      call = function_name
    )
  }

  # Check for minimum of input in expected range
  if(actual_minimum < expected_minimum - tolerance){
    actual_minimum <- signif(actual_minimum, 6)
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Minimum value of '", deparse(substitute(input)),
        "' (", actual_minimum,") does not match expected range(s). Values must range between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#range-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Zero-order correlations to partial correlations ----
# Updated 29.06.2023
cor2pcor <- function(correlation_matrix)
{

  # Convert to inverse correlations to partial correlations
  partial_correlations <- -cov2cor(solve(correlation_matrix))

  # Set diagonal to zero
  diag(partial_correlations) <- 0

  # Return
  return(partial_correlations)

}

#' @noRd
# Partial correlations to zero-order correlations ----
# Updated 29.06.2023
pcor2cor <- function(partial_correlations)
{

  # Set diagonal to negative 1
  diag(partial_correlations) <- -1

  # Return partial correlations as zero-order correlations
  return(cov2cor(solve(-partial_correlations)))

}

#' @noRd
# Create sparse network ----
# Updated 08.11.2023
sparse_network <- function(network)
{

  # Get number of nodes
  nodes <- dim(network)[2]

  # Get node sequence
  node_sequence <- seq_len(nodes)

  # Create data frame
  sparse_df <- data.frame(
    row = rep(node_sequence, each = nodes),
    col = rep(node_sequence, times = nodes),
    weight = as.vector(network)
  )

  # Return lower triangle
  return(sparse_df[sparse_df$row < sparse_df$col,])

}

#' @noRd
# Positive definite matrix ----
# Logical for whether a matrix is positive definite
# Updated 05.11.2024
is_positive_definite <- function(data)
{
  return(all(eigen(data, symmetric = TRUE, only.values = TRUE)$values > .Machine$double.eps))
}

#' @noRd
# Fast kappa exact ----
# Updated 11.03.2026
fast_kappa <- function(X)
{
  # Obtain eigenvalues
  E <- abs(eigen(X, symmetric = TRUE, only.values = TRUE)$values)

  # Return kappa
  return(max(E) / min(E))

}

#' @noRd
# Wrapper for `eigen` ----
# Extracts eigenvalues only
# Updated 11.03.2026
matrix_eigenvalues <- function(X)
{
  return(eigen(X, symmetric = TRUE, only.values = TRUE)$values)
}

#' @noRd
# OS and System Check ----
# Updated 28.11.2023
system.check <- function()
{

  # Get OS usage
  OS <- unname(tolower(Sys.info()["sysname"]))

  # Get RStudio usage
  RSTUDIO <- swiftelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)

  # Return list
  return(
    list(
      OS = OS,
      R = paste0(R.version$major, ".", R.version$minor),
      RSTUDIO = RSTUDIO,
      TEXT = swiftelse(!RSTUDIO & OS != "linux", FALSE, TRUE)
    )
  )

}

#' @noRd
# Style text ----
# Updated 26.07.2023
styletext <- function(
    text, defaults = c(
      "bold", "italics", "highlight",
      "underline", "strikethrough"
    )
)
{

  if(system.check()$TEXT)
  {

    if(missing(defaults)){
      number <- 0
    }else{

      number <- switch(
        defaults,
        bold = 1, italics = 3,
        underline = 4, highlight = 7,
        strikethrough = 9
      )

    }

    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))

  }else{
    return(text)
  }

}

