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
# Generate uniform data ----
# Allows adjustment of range
# Updated 26.10.2023
runif_xoshiro <- function(n, min = 0, max = 1, seed = NULL)
{

  # Get values
  values <- .Call(
    "r_xoshiro_uniform",
    as.integer(n),
    swiftelse(is.null(seed), 0, seed),
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
# Updated 30.07.2023
shuffle <- function(x, size = length(x), seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle",
      as.integer(seq_along(x)),
      swiftelse(is.null(seed), 0, seed),
      PACKAGE = "psychNetwoRk"
    )][seq_len(size)]
  )

}

#' @noRd
# Weighted Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 14.04.2026
weighted_shuffle <- function(x, size = length(x), prob, seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_weighted_shuffle",
      as.integer(seq_along(x)), prob,
      swiftelse(is.null(seed), 0, seed),
      PACKAGE = "psychNetwoRk"
    )][seq_len(size)]
  )

}

#' @noRd
# Shuffle (with replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 10.03.2026
shuffle_replace <- function(x, size = length(x), seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle_replace",
      x, size, swiftelse(is.null(seed), 0, seed),
      PACKAGE = "psychNetwoRk"
    )]
  )

}

#' @noRd
# Random normal generation with Ziggurat ----
# https://people.sc.fsu.edu/~jburkardt/cpp_src/ziggurat/ziggurat.html
# Updated 27.07.2023
rnorm_ziggurat <- function(n, seed = NULL)
{

  # Return call from C
  return(
    .Call(
      "r_ziggurat",
      as.integer(n),
      swiftelse(is.null(seed), 0, seed),
      PACKAGE = "psychNetwoRk"
    )
  )

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
# Generates skewed data for continuous data ----
# Updated 13.03.2026
skew_continuous <- function(skew, data = NULL, sample_size = 1e06, tolerance = 1e-05)
{

  # Generate data
  if(is.null(data)){
    data <- rnorm_ziggurat(sample_size)
  }

  # Check for zero skew (add minimal skew)
  if(skew == 0){
    return(data)
  }

  # Obtain sign and absolute skew
  skew_sign <- sign(skew)
  absolute <- abs(skew)

  # Compute arcsinh once
  arcsinh_data <- asinh(data)

  # Initialize kurtosis
  kurtosis <- 1

  # Optimize for result
  objective <- function(kurtosis){

    # Optimize for skew
    optimize(
      f = function(x){
        abs(absolute - psych::skew(sinh(kurtosis * (arcsinh_data + x))))
      }, interval = c(0, absolute + 2), tol = tolerance
    )

  }

  # Repeat until convergence
  repeat{

    # Get result
    result <- silent_call(objective(kurtosis))

    # Check for minimum too low
    if(is.na(result$objective)){
      return(data)
    }

    # Check for tolerance
    if(result$objective < tolerance){
      break
    }

    # Update kurtosis
    kurtosis <- kurtosis + result$objective * 0.25 # learning rate

  }

  # Return result
  return(skew_sign * scale(sinh(kurtosis * (arcsinh_data + result$minimum))))

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

  # Return default if length > 1 and several are NOT okay
  if(length(argument) > 1 & !several.ok){
    return(default)
  }

  # Get argument name
  argument_name <- as.character(substitute(argument))

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
  ## Expected value
  expected_value <- swiftelse(
    any(expected_value %in% c("integer", "double")),
    "numeric", expected_value
  )

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

  # Obtain maximum and minimum values
  actual_maximum <- round(max(input, na.rm = TRUE), 3)
  actual_minimum <- round(min(input, na.rm = TRUE), 3)

  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum){
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

  # Check for maximum of input in expected range
  if(actual_minimum < expected_minimum){
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

