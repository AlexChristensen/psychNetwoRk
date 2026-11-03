#' Simulates Small-World GGM Data
#'
#' @description
#' Simulates data from a Gaussian Graphical Model (GGM) with a small-world
#' network structure. The generative process proceeds in three stages. First,
#' a ring lattice is constructed with \code{neighbors + 1} nearest-neighbor
#' connections per node, where \code{neighbors} is derived from \code{nodes}
#' and \code{density}. Second, the lattice is randomly pruned to the target
#' \code{density}, introducing degree heterogeneity from the outset. Third,
#' edges are rewired with probability \code{rewire}, where rewired edges are
#' placed preferentially on node pairs with higher combined degree (degree-
#' weighted rewiring). This approach produces more realistic degree
#' distributions than standard Watts-Strogatz rewiring while preserving
#' the local clustering structure of the lattice. Edge weights are assigned
#' by structural priority: edges retained from the pruned lattice receive
#' larger partial correlation weights based on their original neighbor
#' distance, grounded in the empirical observation that shorter-distance
#' connections tend to carry stronger weights in psychometric networks.
#' The resulting network is used to generate multivariate normal (or skewed)
#' data. Parameters do not have default values (except
#' \code{negative_proportion}, \code{snr}, \code{target_condition},
#' \code{max_correlation}, and \code{max_iterations}) and must each be set.
#' See Details and Examples to get started.
#'
#' @param nodes Numeric (length = 1).
#' Number of nodes in the network.
#' Minimum of three nodes. The total number of nodes should be between 8
#' and 54 to remain within the range of the empirical networks used to fit
#' the Weibull parameter model; values outside this range are accepted but
#' will trigger extrapolation and may produce a warning from
#' \code{\link{weibull_parameters}}.
#'
#' @param density Numeric (length = 1).
#' Target edge density of the network after pruning the initial ring lattice.
#' Controls the number of nearest neighbors \code{k} each node is connected
#' to via \eqn{k = \mathrm{round}((\text{nodes} \times
#' \frac{\text{nodes}-1}{2} \times \text{density}) / \text{nodes})}, subject
#' to a minimum of 2 and a maximum of
#' \eqn{\lfloor (\text{nodes}-1)/2 \rfloor}. A ring lattice with
#' \code{neighbors + 1} connections per node is first generated and then
#' pruned to this density, introducing degree heterogeneity before rewiring.
#' Must be between 0 and 1. A minimum density sufficient to maintain a
#' connected graph is enforced; values below this threshold will produce an
#' informative error.
#'
#' @param rewire Numeric (length = 1).
#' Probability of rewiring each edge. Unlike the standard Watts-Strogatz
#' model, rewired edges are placed using degree-weighted random selection:
#' node pairs with higher combined degree have a greater probability of
#' receiving a rewired edge (see Details). Values near 0 preserve the
#' pruned lattice structure; values near 1 produce approximately random
#' networks. The small-world regime typically occurs at intermediate values
#' (roughly 0.01 to 0.30). Must be between 0 and 1.
#'
#' @param snr Numeric (length = 1).
#' Signal-to-noise ratio of the absolute partial correlation weights,
#' defined as \eqn{\bar{|w|} / \mathrm{SD}(|w|)}. Values less than 1
#' produce wider, more heterogeneous weight distributions; values greater
#' than 1 produce narrower, more homogeneous distributions. Empirically
#' observed SNR values ranged from 0.648 to 1.712; values outside this
#' range are accepted but will trigger a warning. Defaults to \code{1}.
#'
#' @param negative_proportion Numeric (length = 1).
#' Proportion of edges that are negative (inhibitory).
#' Must be between 0 and 0.50. The upper bound of 0.50 is a mathematical
#' constraint of the sign-flipping procedure used to assign negative edges
#' (see Details). If not provided, a value is sampled from a truncated
#' normal distribution reflecting the empirical distribution of true
#' negative partial correlations across 194 psychometric networks:
#' mean = 0.34, SD = 0.086, bounded to \eqn{[0.083, 0.50]}.
#'
#' @param sample_size Numeric (length = 1).
#' Number of observations to generate from the population multivariate
#' normal distribution. Also influences the predicted Weibull scale
#' parameter via \code{\link{weibull_parameters}}: larger samples are
#' associated with smaller, more precisely estimated edge weights.
#'
#' @param skew Numeric (length = 1 or \code{nodes}).
#' Skew applied to each variable after generation from the multivariate
#' normal. Can be a single value (applied to all variables) or one value
#' per variable. Values are rounded to the nearest increment of 0.05 in
#' the range \eqn{[-2, 2]}. Defaults to \code{0} (no skew).
#'
#' @param skew_range Numeric (length = 2).
#' If provided, a skew value is drawn uniformly from this range for each
#' variable, overriding \code{skew}. Both values must be between -2 and 2.
#'
#' @param target_condition Numeric (length = 1).
#' Target condition number (using \code{\link{kappa}} with \code{exact =
#' TRUE}) applied when ridge regularization is needed to recover a positive
#' definite precision matrix. The smallest ridge penalty \eqn{\lambda} that
#' brings the condition number to this target is found via root-finding
#' (\code{\link{uniroot}}), subject to a maximum shrinkage of approximately
#' 23\% following Peeters et al. (2020). Lower values produce
#' better-conditioned matrices. Defaults to \code{30}. Values up to
#' \code{100} are accepted but not recommended.
#'
#' @param max_correlation Numeric (length = 1).
#' Maximum allowed absolute pairwise correlation in the population
#' correlation matrix \code{R}. Any draw where
#' \code{max(abs(R[lower.tri(R)])) > max_correlation} is rejected and a
#' new attempt is made. Must be between 0 and 1. Defaults to \code{0.80}.
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of attempts to find (1) a connected network structure,
#' (2) a network satisfying the small-world screening criterion, and (3) a
#' valid set of edge weights. The error message reports a frequency table
#' of rejection reasons to assist with diagnosing convergence failures.
#' Defaults to \code{100}.
#'
#' @details
#' \strong{Lattice generation and pruning}
#'
#' The generative process begins by constructing a ring lattice with
#' \code{neighbors + 1} nearest-neighbor connections per node via
#' \code{\link[igraph]{sample_smallworld}} with \code{p = 0}. The
#' \code{+1} overshoot ensures the lattice always has more edges than the
#' target density, guaranteeing that pruning can proceed. The lattice is
#' then randomly pruned to the target \code{density} by removing edges
#' uniformly at random, subject to the constraint that the pruned graph
#' remains connected. Because removal is random, different nodes lose
#' different numbers of edges, producing a heterogeneous degree distribution
#' before any rewiring occurs. This heterogeneity provides a non-uniform
#' prior for the subsequent degree-weighted rewiring step.
#'
#' \strong{Degree-weighted rewiring}
#'
#' Each edge in the pruned lattice is independently selected for rewiring
#' with probability \code{rewire}. For each selected edge \eqn{(i, j)},
#' node \eqn{i} is kept fixed and the \eqn{j} endpoint is redirected to a
#' new target node. Valid targets are restricted to node pairs involving
#' node \eqn{i} that (1) are not currently connected, and (2) have never
#' been occupied during the current rewiring pass (i.e., were absent in
#' the original pruned lattice). Among valid targets, the new endpoint is
#' selected with probability proportional to \eqn{\sqrt{d_i + d_k}}, where
#' \eqn{d_i} and \eqn{d_k} are the current degrees of the two nodes in
#' the candidate pair. The square-root transformation moderates the
#' rich-get-richer tendency of linear preferential attachment, producing
#' degree heterogeneity consistent with the empirical range of psychometric
#' networks without generating extreme hubs. Node degrees are updated
#' incrementally after each rewire so that subsequent rewiring steps
#' reflect the current graph state.
#'
#' \strong{Edge weight assignment}
#'
#' Edge weights are assigned by ranking edges according to their distance
#' in the pruned lattice before rewiring. Edges retained from the pruned
#' lattice receive their original neighbor distance (1 = nearest neighbor,
#' 2 = second nearest, etc.); rewired edges are assigned a distance of
#' \code{neighbors + 1}, placing them at the bottom of the priority order.
#' Absolute edge weights are drawn from a Weibull distribution whose
#' parameters are predicted from the network size, sample size, and signal-
#' to-noise ratio using a Seemingly Unrelated Regression (SUR) model fitted
#' to 194 empirical psychometric networks (Huth et al., 2025). The sorted
#' (descending) Weibull weights are then mapped onto the distance ranking
#' so that the largest weights go to the shortest-distance edges. This
#' assignment is grounded in the empirical observation that shorter-distance
#' local connections tend to carry larger partial correlation weights in
#' psychometric networks.
#'
#' \strong{Negative edge assignment}
#'
#' Negative edges are introduced by flipping the signs of a subset of
#' nodes — multiplying all edges incident to those nodes by -1. The number
#' of nodes to flip is derived from \code{negative_proportion} via the
#' inversion formula \eqn{(1 - \sqrt{1 - 2p}) / 2}, where \eqn{p} is the
#' target proportion of negative edges. This formula assumes that an edge
#' is negative if and only if exactly one of its endpoints is flipped,
#' giving an expected negative proportion of \eqn{2 \cdot (k/n) \cdot
#' (1 - k/n)} for \eqn{k} flipped nodes out of \eqn{n} total. The formula
#' requires \eqn{p \leq 0.50}, which is why \code{negative_proportion} is
#' bounded at 0.50.
#'
#' @return A named list of class \code{"ggm_simulate"} with four elements
#' (see \code{\link{add_population_error}} to add population error):
#'
#' \item{data}{Numeric matrix of dimension \code{sample_size x nodes}
#' containing the simulated observations drawn from the population GGM.
#' Rows are cases; columns are variables named \code{V01}, \code{V02},
#' etc. Values are continuous (or skewed continuous when \code{skew != 0}).
#' To produce ordinal data, pass the columns through
#' \code{\link{categorize}}.}
#'
#' \item{parameters}{
#' A list of input, derived, and estimated parameters:
#' \itemize{
#'   \item \code{nodes} --- Number of nodes (as supplied)
#'   \item \code{density} --- Target edge density (as supplied)
#'   \item \code{neighbors} --- Number of nearest neighbors \code{k} derived
#'   from \code{nodes} and \code{density}; the initial ring lattice uses
#'   \code{neighbors + 1} connections before pruning
#'   \item \code{rewire} --- Edge rewiring probability (as supplied)
#'   \item \code{negative_proportion} --- Proportion of negative edges,
#'   either as supplied or as sampled from the empirical distribution
#'   \item \code{sample_size} --- Number of simulated observations
#'   \item \code{skew} --- Named numeric vector of per-variable skew values
#'   actually applied (after rounding and possible resampling)
#'   \item \code{weibull} --- Named numeric vector of length 2 giving the
#'   Weibull \code{shape} and \code{scale} parameters of the absolute edge
#'   weight distribution actually used. If ridge conditioning was applied,
#'   these are re-estimated from the conditioned network via MLE.
#'   \item \code{omega} --- Smallworldness omega statistic of the generated
#'   network (Telesford et al., 2011). Values near zero indicate small-world
#'   structure; negative values indicate lattice-like structure; positive
#'   values indicate random-like structure (see \code{\link[psychNetwoRk]{smallworldness}})
#'   \item \code{Q} --- Newman-Girvan modularity of the population network
#'   (\code{Omega}) computed via \code{igraph::modularity} on absolute edge
#'   weights, using the community partition returned by
#'   \code{igraph::cluster_leiden} with \code{objective_function =
#'   "modularity"}. Because \code{cluster_leiden} is a heuristic algorithm,
#'   \code{Q} reflects a high-quality but not necessarily optimal partition.
#'   Values near zero indicate absence of community structure, consistent
#'   with the network theory of psychopathology.
#' }
#' }
#'
#' \item{population}{
#' Population-level network parameters:
#' \itemize{
#'   \item \code{R} --- Population correlation matrix derived from the GGM
#'   via \code{pcor2cor}
#'   \item \code{Omega} --- Population partial correlation matrix (the GGM
#'   edge weight matrix), with zeros for absent edges
#' }
#' }
#'
#' \item{convergence}{
#' Iteration and conditioning diagnostics:
#' \itemize{
#'   \item \code{iterations} --- Number of sampling attempts needed to find
#'   a valid network
#'   \item \code{rejections} --- Character vector recording the rejection
#'   reason for each failed attempt. Common reasons include disconnected
#'   graph structure, condition number exceeding \code{target_condition},
#'   and maximum correlation exceeding \code{max_correlation}
#'   \item \code{lambda} --- Ridge regularization parameter \eqn{\lambda}
#'   added to the diagonal of the precision matrix to ensure positive
#'   definiteness; \code{NA} if no conditioning was required
#'   \item \code{condition} --- Condition number of the final population
#'   correlation matrix \code{R}
#' }
#' }
#'
#' @examples
#' # Basic small-world network (moderate density, moderate rewiring)
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.20,
#'   sample_size = 500
#' )
#'
#' # Lattice-like structure (low rewiring preserves local connectivity)
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.01,
#'   sample_size = 500
#' )
#'
#' # Random-like structure (high rewiring destroys lattice regularity)
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.80,
#'   sample_size = 500
#' )
#'
#' # Fix the proportion of negative edges
#' result <- simulate_smallworld(
#'   nodes = 20,
#'   density = 0.30,
#'   rewire = 0.20,
#'   sample_size = 500,
#'   negative_proportion = 0.20
#' )
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{Simulation approach for small-world networks} \cr
#' Christensen, A. P., Choi, J., Yang, H., & Xue, L. (2026).
#' Adaptive regularization via extreme value distributions for Gaussian
#' graphical models (arXiv:2609.05678). \emph{arXiv}.
#' \doi{10.48550/arXiv.2609.05678}
#'
#' \strong{Seminal introduction to the Watts-Strogatz small-world model} \cr
#' Watts, D. J., & Strogatz, S. H. (1998).
#' Collective dynamics of 'small-world' networks.
#' \emph{Nature}, \emph{393}(6684), 440--442.
#'
#' \strong{Logic for weight assignments} \cr
#' Muldoon, S. F., Bridgeford, E. W., & Bassett, D. S. (2016).
#' Small-world propensity and weighted brain networks.
#' \emph{Scientific Reports}, \emph{6}(1), 22057.
#'
#' \strong{Analytical approximation of random-graph average path length} \cr
#' Newman, M. E. J., Strogatz, S. H., & Watts, D. J. (2001).
#' Random graphs with arbitrary degree distributions and their applications.
#' \emph{Physical Review E}, \emph{64}(2), 026118.
#'
#' \strong{Empirical network data used to fit the Weibull SUR model} \cr
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., &
#' Marsman, M. (2025). Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#'
#' \strong{Maximum ridge shrinkage bound} \cr
#' Peeters, C. F., van de Wiel, M. A., & van Wieringen, W. N. (2020).
#' The spectral condition number plot for regularization parameter evaluation.
#' \emph{Computational Statistics}, \emph{35}(2), 629--646.
#'
#' @export
#'
# Simulate Small-world GGM data ----
# Updated 15.05.2026
simulate_smallworld <- function(
    nodes, density, rewire, snr = 1, negative_proportion,
    sample_size, skew = 0, skew_range = NULL,
    target_condition = 30, max_correlation = 0.80,
    max_iterations = 100
)
{

  # Check for missing negative proportion
  if(missing(negative_proportion)){

    # Compute proportion of negative edges based on empirical values
    negative_proportion <- pmin(
      pmax(
        0.34 + rnorm_ziggurat(1) * 0.086, # empirical mean +/- 1 SD
        0.083 # empirical minimum
      ), 0.50 # empirical maximum (0.54945055)
      # needs to be capped at 0.50 due to flipping variables
    )

  }

  # Check for input errors
  simulate_smallworld_errors(
    nodes, density, rewire, snr, negative_proportion, sample_size,
    skew, target_condition, max_correlation, max_iterations
  )

  # Set neighbors
  neighbors <- get_neighbors(nodes, density)

  # Create lattice network
  lattice <- lattice_generate(nodes, neighbors + 1)

  # Obtain lower triangle
  lower_triangle <- lower.tri(lattice)

  # Initialize iterations
  iter <- 0

  # Initialize condition
  search <- TRUE

  # Collect rejections
  rejections <- character(length = max_iterations + 1)

  # Generation loop
  while(search){

    # Increase iterations
    iter <- iter + 1

    # Stop on greater than max iterations
    if(iter > max_iterations){

      # Collect rejections
      rejection_table <- fast_table(rejections)

      # Remove NULL
      rejection_table <- rejection_table[names(rejection_table) != ""]

      # Return error
      stop(
        paste0(
          "Reached maximum iterations: Could not find weights that converge for the small-world structure.\n",
          "Resulting rejections were due to:\n\n",
          paste0(names(rejection_table), " = ", rejection_table, collapse = "\n")
        )
      )

    }

    # Prune lattice
    pruned_lattice <- prune2density(lattice, density, nodes, lower_triangle)

    # Small-world adjacency
    A <- smallworld_generate(pruned_lattice, rewire, lower_triangle)

    # Check adjacency
    if(!igraph::is_connected(convert2igraph(A))){

      # Add rejection reason
      rejections[iter] <- "Could not find structure where every node was connected."

      # Move to next iteration
      next

    }

    # Try to get good weights
    output <- try(
      smallworld_weights(
        A, pruned_lattice, nodes, sample_size,
        snr, negative_proportion, target_condition,
        max_correlation, lower_triangle
      ), silent = TRUE
    )

    # Check for error
    if(is(output, "try-error")){

      # Add rejection reason
      rejections[iter] <- trimws(strsplit(output[1], split = "\n")[[1]][2])

      # Move to next iteration
      next

    }

    # Update search
    search <- FALSE



  }

  # Generate data
  data_output <- simulate_data(
    n = sample_size, p = nodes, R = output$R, skew = skew, skew_range = skew_range
  )

  # Estimate modularity using optimal modularity detection
  igraph_network <- convert2igraph(abs(output$network))

  # Populate results
  results <- list(
    data = data_output$data,
    parameters = list(
      nodes = nodes,
      density = density,
      neighbors = neighbors,
      rewire = rewire,
      negative_proportion = negative_proportion,
      sample_size = sample_size,
      skew = data_output$skew,
      weibull = output$params,
      omega = smallworldness(
        network = output$network, lattice = pruned_lattice, method = "omega"
      ),
      Q = igraph::modularity(
        igraph_network, igraph::cluster_leiden(
          igraph_network, objective_function =  "modularity"
        )$membership
      )
    ),
    population = list(R = output$R, Omega = output$network),
    convergence = list(
      iterations = iter,
      rejections = rejections,
      lambda = output$lambda,
      condition = output$condition
    )
  )

  # Add class
  class(results) <- "ggm_simulate"

  # Return results
  return(results)

}

# Bug checking ----
# nodes = 20; density = 0.50; rewire = 0.10
# negative_proportion = 0.35; sample_size = 1000
# skew = 0; skew_range = NULL; snr = 1
# target_condition = 30;
# max_correlation = 0.80; max_iterations = 100

#' @noRd
# Errors ----
# Updated 16.03.2026
simulate_smallworld_errors <- function(
    nodes, density, rewire, snr, negative_proportion, sample_size,
    skew, target_condition, max_correlation, max_iterations
)
{

  # Errors for 'nodes'
  typeof_error(nodes, "numeric", "simulate_smallworld")
  length_error(nodes, 1, "simulate_smallworld")
  range_error(nodes, c(3, Inf), "simulate_smallworld")

  # Errors for 'density'
  typeof_error(density, "numeric", "simulate_smallworld")
  length_error(density, 1, "simulate_smallworld")
  range_error(density, c(0, 1), "simulate_smallworld")

  # Minimum edges required for a connected graph
  min_edges <- nodes - 1
  target_edges <- round(density * (nodes * (nodes - 1) / 2))

  # Check for whether connected graph is possible
  if(target_edges < min_edges){
    stop(
      paste0(
        "Target density (", round(density, 3), ") is too low to maintain ",
        "a connected graph for ", nodes, " nodes. ",
        "Minimum required density: ", round(min_edges / (nodes * (nodes - 1) / 2), 3)
      )
    )
  }

  # Errors for 'rewire'
  typeof_error(rewire, "numeric", "simulate_smallworld")
  length_error(rewire, 1, "simulate_smallworld")
  range_error(rewire, c(0, 1), "simulate_smallworld")

  # Error for 'snr'
  typeof_error(snr, "numeric", "simulate_smallworld")
  length_error(snr, 1, "simulate_smallworld")
  range_error(snr, c(0.01, Inf), "simulate_smallworld")

  # Send warning if beyond bounds
  if((snr < 0.64) | (snr > 1.72)){
    warning(paste0(
      "The signal-to-ratio specified (", round(snr, 2), ") is beyond ",
      "the bounds of cataloged empirical data (0.64-1.72). "
    ))
  }

  # Errors for 'negative_proportion'
  typeof_error(negative_proportion, "numeric", "simulate_smallworld")
  length_error(negative_proportion, 1, "simulate_smallworld")
  range_error(negative_proportion, c(0, 0.50), "simulate_smallworld")

  # Errors for 'sample_size'
  typeof_error(sample_size, "numeric", "simulate_smallworld")
  length_error(sample_size, 1, "simulate_smallworld")
  range_error(sample_size, c(1, Inf), "simulate_smallworld")

  # Errors for 'skew'
  typeof_error(skew, "numeric", "simulate_smallworld")
  length_error(skew, c(1, nodes), "simulate_smallworld")
  range_error(skew, c(-2, 2), "simulate_smallworld")

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric", "simulate_smallworld")
  length_error(target_condition, 1, "simulate_smallworld")
  range_error(target_condition, c(1, 100), "simulate_smallworld")

  # Errors for 'max_correlation'
  typeof_error(max_correlation, "numeric", "simulate_smallworld")
  length_error(max_correlation, 1, "simulate_smallworld")
  range_error(max_correlation, c(0, 1), "simulate_smallworld")

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric", "simulate_smallworld")
  length_error(max_iterations, 1, "simulate_smallworld")
  range_error(max_iterations, c(1, Inf), "simulate_smallworld")

}

#' @noRd
# Get number of neighbors based on density ----
# Updated 08.03.2026
get_neighbors <- function(nodes, density)
{

  # Degrees of freedom
  half_df <- (nodes - 1) / 2

  # Total edges with density
  edges <- nodes * half_df * density

  # Return with minimum neighbor check
  return(
    max(
      min(
        max(round(edges / nodes), 1), # maximum neighbors
        floor(half_df) # minimum degrees of freedom
      ), 2 # at least two neighbors
    )
  )

}

#' @noRd
# Generate lattice ----
# Updated 18.03.2026
lattice_generate <- function(nodes, neighbors)
{

  # Set node sequence
  node_sequence <- seq_len(nodes)

  # Set up distance matrix
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))
  distance_matrix <- pmin(distance_matrix, nodes - distance_matrix)

  # Convert distance matrix to zeros and ones
  distance_matrix[] <- as.numeric(distance_matrix <= neighbors)

  # Ensure zero diagonal
  diag(distance_matrix) <- 0

  # Return lattice
  return(distance_matrix)

}

#' @noRd
# Prune neighbors to density ----
# Updated 13.04.2026
prune2density <- function(lattice, density, nodes, lower_triangle)
{

  # Store original lattice
  original_lattice <- lattice

  # Initialize lower triangle
  lattice_copy <- original_lattice[lower_triangle]

  # Compute target edges
  target_edges <- round(
    density * (nodes * (nodes - 1) / 2)
  )

  # Initilaize connectedness
  connected <- FALSE

  # Connecteness loop
  while(!connected){

    # Initialize lower lattice
    lattice_lower <- lattice_copy

    # Obtain current edges
    nonzero <- which(lattice_lower != 0)
    current_edges <- length(nonzero)

    # Obtain number of nodes to remove
    remove_edges <- current_edges - target_edges

    # Check for whether edges is zero or negative
    if(remove_edges < 1){
      return(lattice)
    }

    # Select edges and remove them
    lattice_lower[shuffle(nonzero, remove_edges)] <- 0

    # Reconstruct lattice
    lattice[] <- 0
    lattice[lower_triangle] <- lattice_lower
    lattice <- lattice + t(lattice)

    # Check connectedness
    connected <- igraph::is_connected(convert2igraph(lattice))

  }

  # Return lattice
  return(lattice)

}

#' @noRd
# Generate smallworld ----
# Updated 16.03.2026
smallworld_generate <- function(lattice, rewire, lower_triangle)
{

  # Obtain degree
  degree <- colSums(lattice, na.rm = TRUE)

  # Create sparse network data frame
  sparse_lattice <- sparse_network(lattice)

  # Add degrees
  sparse_lattice$degree <- outer(degree, degree, "+")[lower_triangle]

  # Determine rewired edges
  current_edges <- sparse_lattice$weight != 0
  edge_indices <- which(current_edges)
  valid_targets <- which(!current_edges)
  rewire_edges <- edge_indices[runif_xoshiro(length(edge_indices)) < rewire]
  n_rewire <- length(rewire_edges)

  # Check for no rewiring
  if(n_rewire > 0){

    # Loop over edges to rewire
    for(edge in rewire_edges){

      # Check for no more updates
      if(length(valid_targets) == 0){
        break
      }

      # Get edge indices
      node_i <- sparse_lattice$row[edge]
      node_j <- sparse_lattice$col[edge]

      # Get candidate edges
      candidate_edges <- intersect(
        which((sparse_lattice$row == node_i) | (sparse_lattice$col == node_i)),
        valid_targets
      )

      # Check for candidate edges
      if(length(candidate_edges) == 0){
        next
      }

      # Get weights based on degree
      zero_degree <- sqrt(sparse_lattice$degree[candidate_edges])

      # Get rewire location
      location <- weighted_shuffle(x = candidate_edges, size = 1, prob = zero_degree / sum(zero_degree))

      # Update edges
      sparse_lattice$weight[edge] <- 0
      sparse_lattice$weight[location] <- 1

      # Get rewired index
      rewired_j <- swiftelse(
        sparse_lattice$row[location] == node_i,
        sparse_lattice$col[location],
        sparse_lattice$row[location]
      )

      # Obtain indices
      node_j_index <- (sparse_lattice$row == node_j) | (sparse_lattice$col == node_j)
      rewired_j_index <- (sparse_lattice$row == rewired_j) | (sparse_lattice$col == rewired_j)

      # Update degrees
      sparse_lattice$degree[node_j_index] <- sparse_lattice$degree[node_j_index] - 1
      sparse_lattice$degree[rewired_j_index] <- sparse_lattice$degree[rewired_j_index] + 1

      # Update valid_targets
      valid_targets <- setdiff(valid_targets, location)

    }

  }

  # Obtain nonzeros
  nonzero_network <- sparse_lattice[sparse_lattice$weight == 1,]

  # Rebuild network
  lattice[] <- 0
  lattice[cbind(nonzero_network$row, nonzero_network$col)] <-
  lattice[cbind(nonzero_network$col, nonzero_network$row)] <- 1


  # Return network
  return(lattice)

}

#' @noRd
# Smallworld weight generation ----
# Updated 13.04.2026
smallworld_weights <- function(
    A, lattice, nodes, sample_size,
    snr, negative_proportion, target_condition,
    max_correlation, lower_triangle
)
{

  # Set node sequence
  node_sequence <- seq_len(nodes)

  # Set up distance matrix
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))
  distance_matrix <- pmin(distance_matrix, nodes - distance_matrix)

  # Set edges
  lattice_edges <- lattice == 1
  A_edges <- A == 1

  # Set up smallworld distances
  smallworld <- matrix(0, nrow = nodes, ncol = nodes)

  # Get sorted rewired distances
  rewired_distances <- sort(distance_matrix[lower_triangle][(lattice_edges & !A_edges)[lower_triangle]])

  # Set index for rewired indices
  rewired_index <- !lattice_edges & A_edges

  # Set up lower triangle with rewired distances
  smallworld[lower_triangle][rewired_index[lower_triangle]] <- rewired_distances[
    rank(distance_matrix[lower_triangle][rewired_index[lower_triangle]], ties.method = "random")
  ]

  # Make symmetric
  smallworld <- smallworld + t(smallworld)

  # Multiply retained by distance matrix
  smallworld <- smallworld + ((lattice_edges & A_edges) * distance_matrix)

  # Obtain nonzero edges
  nonzero <- A[lower_triangle] != 0

  # Initialize network
  network <- matrix(0, nrow = nodes, ncol = nodes)

  # Generate edges (returned edges are sorted)
  edges <- generate_edges(nonzero = sum(nonzero), n = sample_size, p = nodes, snr = snr)

  # Set weights order
  # Follows: Muldoon, Bridgeford, & Bassett's (2016) implementation
  weight_order <- rank(smallworld[lower_triangle][nonzero], ties.method = "random")
  network[lower_triangle][nonzero] <- edges[weight_order]
  network <- network + t(network) # make symmetric

  # Get signs
  network <- determine_signs(network, nodes, negative_proportion)

  # Get population values
  R <- silent_call(pcor2cor(network))

  # Set lambda
  lambda <- NA

  # Check for positive definite
  if(anyNA(R) || !is_positive_definite(R)){

    # Try to condition network
    output <- try(condition_network(network, target_condition), silent = TRUE)

    # Check for whether correlation matrix is conditioned
    if(is(output, "try-error")){
      stop("Ill-conditioned")
    }

    # Store outputs
    lambda <- output$lambda

    # Update network and correlations
    P <- cor2pcor(output$R)
    P[network == 0] <- 0
    network <- P
    R <- pcor2cor(network)

    # Check for positive definite again
    if(anyNA(R) || !is_positive_definite(R)){
      stop("Not positive definite")
    }

  }

  # Check for condition
  condition <- fast_kappa(R)
  if(condition > (target_condition + 1)){
    stop("Condition number exceeds target")
  }

  # Check for maximum correlation
  if(max(abs(R[lower_triangle])) > max_correlation){
    stop("Max correlation exceeds target")
  }

  # Return results
  return(
    list(
      R = R,
      network = network,
      params = attr(edges, "params"),
      lambda = lambda,
      condition = condition
    )
  )

}

#' @noRd
# Determine balance of signs ----
# Updated 11.03.2026
determine_signs <- function(network, nodes, negative_proportion)
{

  # Check for zero negative weights
  if(negative_proportion == 0){
    return(network)
  }

  # Convert negative_proportion (edge-level) to node flip proportion
  flip_proportion <- (1 - sqrt(1 - 2 * negative_proportion)) / 2
  flip_number <- round(nodes * flip_proportion)

  # Flip at least one
  flip_number <- max(1, min(flip_number, floor(nodes / 2)))

  # Determine node indices
  index <- shuffle(seq_len(nodes), flip_number)

  # Flip nodes
  network[index,] <- -network[index,]
  network[,index] <- -network[,index]

  # Return the network
  return(network)

}