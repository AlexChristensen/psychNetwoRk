#' @title Construct a Degree-Preserving Ring Lattice via Proximity-Swap Construction
#'
#' @description Converts a network matrix into a connected ring lattice whose
#' degree sequence exactly matches the original, while maximizing the average
#' local clustering coefficient. The resulting lattice is intended as a null
#' model for small-world network analyses. An adjacency matrix is derived
#' automatically from \code{network} (non-zero entries become edges), so the
#' function accepts any weighted or binary network directly.
#'
#' Each pass randomly permutes the degree sequence onto ring positions, then
#' greedily assigns edges in strictly increasing ring-distance order (proximity
#' construction). Any residual degree deficit is resolved by a swap-repair
#' phase. Passes that yield a disconnected graph or an unsatisfied degree
#' sequence are discarded; among valid passes the one with the highest average
#' clustering coefficient is returned. \code{shuffles} independent passes are
#' attempted in total.
#'
#' @param network Matrix.
#' A square, symmetric numeric matrix representing a network (e.g., partial
#' correlations). Non-zero off-diagonal entries are treated as edges; the
#' binary adjacency is derived internally as \code{network != 0}.
#' Isolated nodes (degree zero) are supported.
#'
#' @param weighted Logical (length = 1).
#' Whether to return a weighted ring lattice. When \code{TRUE}, the edge
#' weights from \code{network} are reassigned to the lattice edges according
#' to ring distance following the implementation of Muldoon, Bridgeford, &
#' Bassett (2016): shorter-distance lattice edges receive larger weights,
#' preserving the overall weight distribution while concentrating stronger
#' connections locally. When \code{FALSE} (default), a binary adjacency
#' matrix is returned.
#'
#' @param shuffles Numeric (length = 1).
#' Number of independent random permutation passes to attempt. Each pass
#' assigns the observed degree sequence to ring positions in a new random
#' order, runs the proximity construction, and applies swap repair if needed.
#' Only passes producing a connected graph with zero degree error are
#' retained; the one with the highest clustering coefficient is returned.
#' Defaults to \code{100}.
#'
#' @return A square symmetric matrix of the same dimension as \code{network},
#' with row and column names preserved, representing the resulting ring
#' lattice. When \code{weighted = FALSE} (default), entries are binary integers
#' (\code{0L}/\code{1L}); when \code{weighted = TRUE}, entries contain the
#' reassigned edge weights from \code{network} with original signs preserved.
#' The average clustering coefficient of the returned graph is attached as the
#' attribute \code{"CC"} and can be retrieved with \code{attr(result, "CC")}.
#' When the empirical fallback is triggered (see Details), \code{"CC"} reflects
#' the empirical network's clustering coefficient rather than a lattice's.
#'
#' @details
#' ## Algorithm
#'
#' \strong{Pair precomputation.} A ring distance matrix is computed as
#' \eqn{d_{ij} = \min(|i - j|,\, n - |i - j|)}, giving true circular
#' distances bounded by \eqn{\lfloor n/2 \rfloor}. All unique unordered pairs
#' at each distance \eqn{r = 1, \ldots, \lfloor n/2 \rfloor} are extracted
#' once (retaining only entries where \eqn{i < j}) and cached for reuse
#' across all passes.
#'
#' \strong{Proximity construction.} The observed degree sequence is randomly
#' permuted onto ring positions. Starting from an empty graph, edges are added
#' in distance order. At each distance band the eligible pairs (both endpoints
#' have remaining degree budget and are not yet connected) are sorted by
#' descending \eqn{\min(\text{budget}_i, \text{budget}_j)} so that high-need
#' pairs receive short connections first. Pairs are then assigned sequentially
#' with per-pair budget re-checks. The pass exits early once all budgets reach
#' zero. This phase is implemented in compiled C via \code{proximity_pass_c}.
#'
#' \strong{Swap repair.} If any degree deficit remains after the proximity
#' pass, the highest-deficit node \eqn{i} is connected to its nearest
#' available ring neighbour, scanning clockwise and counter-clockwise
#' positions in interleaved distance order. When no direct partner with
#' remaining budget exists, an edge swap is performed: a nearby unconnected
#' node \eqn{j} is found, one of \eqn{j}'s existing edges to \eqn{k} is
#' removed, and a new edge \eqn{(i, j)} is added. Node \eqn{k} recovers its
#' budget for resolution in a subsequent iteration. This repeats until all
#' deficits are resolved or the iteration cap (\eqn{2n^2}) is reached. This
#' phase is implemented in compiled C via \code{swapping_pass_c}.
#'
#' \strong{Weight assignment.} When \code{weighted = TRUE}, edge weights are
#' reassigned after the binary topology is finalised. The observed weights are
#' sorted by \emph{absolute value} in descending order (largest magnitude first)
#' and mapped onto lattice edges ranked by ascending ring distance (i.e.,
#' \eqn{d_{ij} = \min(|i - j|,\, n - |i - j|)}), so that shorter (more local)
#' connections receive the largest-magnitude weights. Original signs are
#' preserved: the sorted weight vector — not its absolute values — is assigned
#' to the lattice edges. This directly implements the distance-weight principle
#' of Muldoon, Bridgeford, & Bassett (2016), using the ring's structural
#' distances rather than any network-derived proxy. Ties in ring distance are
#' broken at random.
#'
#' \strong{Pass selection.} A pass is valid only if the resulting graph is
#' connected and has zero residual degree error. Among valid passes, the one
#' with the highest average clustering coefficient is returned.
#'
#' \strong{Empirical fallback.} If no valid pass is found, or if the best
#' lattice clustering coefficient is lower than that of the original network,
#' the empirical adjacency (or weighted matrix, if \code{weighted = TRUE}) is
#' returned with a warning.
#'
#' @references
#' \strong{Logic for weight assignments} \cr
#' Muldoon, S. F., Bridgeford, E. W., & Bassett, D. S. (2016).
#' Small-world propensity and weighted brain networks.
#' \emph{Scientific Reports}, \emph{6}(1), 22057.
#'
#' @examples
#' # Simulate a small-world network (population partial correlations)
#' network <- simulate_smallworld(
#'   nodes = 20, density = 0.30, rewire = 0.20, sample_size = 500
#' )$population$Omega
#'
#' # Construct binary ring lattice
#' L <- proxswap_lattice(network)
#'
#' # Retrieve the attached clustering coefficient
#' attr(L, "CC")
#'
#' # Degree sequences should match exactly
#' cbind(target = colSums(network != 0), achieved = colSums(L))
#'
#' # Construct weighted ring lattice
#' L_weighted <- proxswap_lattice(network, weighted = TRUE)
#'
#' # Retrieve the attached clustering coefficient
#' attr(L_weighted, "CC")
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Proximity-swap lattice construction ----
# Updated 26.03.2026
proxswap_lattice <- function(network, weighted = FALSE, shuffles = 100)
{

  # Automatically construct adjacency
  A <- network != 0

  # Check for weighted
  if(!weighted){
    network <- A
  }

  # Initialize nodes, degree, and connectedness flag
  nodes  <- dim(network)[2]
  degree <- colSums(A)
  check_connectedness <- all(degree != 0)

  # Set up distance matrix
  node_sequence <- seq_len(nodes)
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))
  distance_matrix <- pmin(distance_matrix, nodes - distance_matrix)

  # Pre-compute unique node pairs at each ring distance, and flatten them
  # into the parallel vectors `proximity_pass` needs once here -- `pairs`
  # doesn't change across the `shuffles` loop below, so reflattening it on
  # every pass (as `proximity_pass` used to do internally) was pure waste
  pairs <- build_pairs(distance_matrix)
  pair_rows <- ulapply(pairs, function(p) p[, "row"])
  pair_cols <- ulapply(pairs, function(p) p[, "col"])
  pair_counts <- lengths(lapply(pairs, function(p) p[, "row"]))

  # Compute empirical clustering coefficient for fallback check
  empirical_CC <- igraph::transitivity(convert2igraph(network), type = "average")

  # Initialize best result trackers
  best_swap <- best_ring <- NULL
  best_CC   <- -Inf

  # Initialize to an empty graph
  ring <- matrix(FALSE, nrow = nodes, ncol = nodes)

  # Run each shuffle pass
  for(k in seq_len(shuffles)){

    # Store swap order (returns proper degree sequence)
    swap_order <- shuffle(node_sequence)

    # Run proximity construction with shuffled degrees
    result <- proximity_pass(
      nodes = nodes, ring = ring, budget = degree[swap_order],
      pair_rows = pair_rows, pair_cols = pair_cols, pair_counts = pair_counts
    )

    # Perform swapping to
    result <- swapping_pass(
      nodes = nodes, ring = result$ring, budget = result$budget,
      total_budget = result$total_budget
    )

    # Check for deficit before expensive connectedness check
    if(result$remaining > 0){
      next
    }

    # Convert ring to {igraph}
    iring <- convert2igraph(result$ring)

    # Determine whether connectedness (and then check it)
    if(check_connectedness && (!igraph::is_connected(iring))){
      next
    }

    # Check for weighted
    if(weighted){

      # Get weights
      result$ring <- assign_weights(network, result$ring, distance_matrix)

      # Update ring to {igraph}
      iring <- convert2igraph(result$ring)

    }

    # With success, compute clustering coefficient for this pass
    pass_CC <- igraph::transitivity(iring, type = "average")

    # Update best pass
    if(pass_CC > best_CC){
      best_CC <- pass_CC; best_ring <- result$ring; best_swap <- swap_order
    }

  }

  # Determine whether solution was reached
  if(is.infinite(best_CC)){

    # Send warning
    warning(paste0(
      "The lattice solution did not converge. The empirical graph was returned.\n\n",
      "Try increasing `shuffles` to find a suitable solution"
    ))

    # Set empirical flag
    empirical_flag <- TRUE

    # Check for empirical and weighted
    ring <- swiftelse(weighted, network, A)

  }else{

    # Check whether empirical beats the best lattice
    empirical_flag <- empirical_CC > best_CC

    # Warn on empirical fallback
    if(empirical_flag){
      warning(paste0(
        "The lattice solution did not produce a better CC (", round(best_CC, 3),
        ") than the empirical (", round(empirical_CC, 3), ").\n",
        "Falling back to empirical solution..."
      ))
    }

    # Select between lattice and empirical
    original_order <- order(best_swap)

    # Check for empirical and weighted
    ring <- swiftelse(
      empirical_flag,
      swiftelse(weighted, network, A),
      best_ring[original_order, original_order]
    )

  }

  # Ensure named matrix
  dimnames(ring) <- dimnames(network)

  # Attach clustering coefficient
  attr(ring, "CC") <- swiftelse(empirical_flag, empirical_CC, best_CC)

  # Return ring
  return(ring)

}

#' @noRd
# Pre-compute unique node pairs per ring distance ----
# Updated 26.03.2026
build_pairs <- function(distance_matrix)
{

  # Loop over distances
  return(
    lapply(seq_len(max(distance_matrix)), function(i){

      # Obtain pairs
      pairs <- which(distance_matrix == i, arr.ind = TRUE)

      # Return pairs
      return(pairs[pairs[,"row"] < pairs[,"col"],])

    }
    )
  )

}

#' @noRd
# Proximity construction ----
# Updated 20.07.2026
# `pair_rows`/`pair_cols`/`pair_counts` are the flattened form of `pairs`
# (see `build_pairs`), computed once by the caller since they don't change
# across the `shuffles` loop in `proxswap_lattice`
proximity_pass <- function(nodes, ring, budget, pair_rows, pair_cols, pair_counts)
{

  # # Remaining edges to assign per node
  # total_budget <- sum(budget)
  #
  # # Add edges in increasing ring distance
  # for(i in distance_sequence){
  #
  #   # Retrieve precomputed node-index vectors for this distance band
  #   rows <- pairs[[i]][,"row"]
  #   columns <- pairs[[i]][,"col"]
  #
  #   # Determine eligible
  #   eligible <- (budget[rows] > 0) & (budget[columns] > 0) & (!ring[cbind(rows, columns)])
  #
  #   # Skip this distance band if nothing is eligible
  #   if(!any(eligible)){
  #     next
  #   }
  #
  #   # Eligible candidate indices
  #   row <- rows[eligible]
  #   col <- columns[eligible]
  #
  #   # Get row length
  #   row_length  <- length(row)
  #
  #   # Sort: highest min-budget pair first so high-need pairs get short edges
  #   if(row_length > 1){
  #     target <- order(pmin(budget[row], budget[col]), decreasing = TRUE)
  #     row  <- row[target]
  #     col  <- col[target]
  #   }
  #
  #   # Assign edges one at a time; re-check both budgets before each assignment
  #   for(l in seq_len(row_length)){
  #
  #     # Set targets
  #     target_row <- row[l]
  #     target_col <- col[l]
  #
  #     # Check budget (can change with loop)
  #     if((budget[target_row] > 0) && (budget[target_col] > 0)){
  #
  #       # Add to ring
  #       ring[target_row, target_col] <- ring[target_col, target_row] <- TRUE
  #
  #       # Update budget
  #       budget[target_row] <- budget[target_row] - 1
  #       budget[target_col] <- budget[target_col] - 1
  #       total_budget <- total_budget - 2
  #
  #     }
  #   }
  #
  #   # Early exit once all budgets are satisfied
  #   if(total_budget < 1){
  #     break
  #   }
  #
  # }
  #
  # # Return ring and total absolute degree error
  # return(list(ring = ring, budget = budget, total_budget = total_budget))

  # Coerce logical ring matrix to integer storage before any .Call
  storage.mode(ring) <- "integer"

  # Collect result
  result <- .Call(
    "proximity_pass_c",
    as.integer(nodes),
    as.integer(ring),
    as.integer(budget),
    as.integer(pair_rows),
    as.integer(pair_cols),
    as.integer(pair_counts),
    PACKAGE = "psychNetwoRk"
  )

  # Set up ring matrix
  result$ring <- matrix(result$ring, nrow = nodes, ncol = nodes)

  # Return result
  return(result)

}

#' @noRd
# Swap to amend lattice ----
# Updated 24.03.2026
swapping_pass <- function(nodes, ring, budget, total_budget)
{

  # # Set maximum iterations
  # for(iter in seq_len(nodes * nodes * 2)){
  #
  #   # Stop when all budgets are satisfied
  #   if(total_budget == 0){
  #     break
  #   }
  #
  #   # Work on the node with the largest remaining deficit
  #   i <- which.max(budget)
  #   if(budget[i] == 0){
  #     break
  #   }
  #
  #   # Interleave clockwise and counter-clockwise ring positions by distance
  #   distance <- interleave(i, nodes, distance_sequence)
  #
  #   # Set avaiable ring (diagonal will be TRUE but will never be in 'distance')
  #   available_ring <- !ring
  #
  #   # Try direct connection to nearest node that has budget and no existing edge
  #   direct <- distance[(budget[distance] > 0) & (available_ring[i, distance])]
  #
  #   # Check connection
  #   if(length(direct) > 0){
  #
  #     # Set index
  #     target <- direct[1]
  #
  #     # Update ring
  #     ring[i, target] <- ring[target, i] <- TRUE
  #
  #     # Update budget
  #     budget[i] <- budget[i] - 1
  #     budget[target] <- budget[target] - 1
  #     total_budget <- total_budget - 2
  #
  #     # Continue
  #     next
  #
  #   }
  #
  #   # Initialize swap
  #   swapped <- FALSE
  #
  #   # Check nearby unconnected nodes
  #   candidates <- distance[available_ring[i, distance]]
  #
  #   # Loop over candidates
  #   for(j in candidates){
  #
  #     # Find a neighbor of the neighbor that is not connected to current node
  #     neighborhood <- which(ring[j,] & node_sequence != i & available_ring[i,])
  #
  #     # Check if any are available
  #     if(length(neighborhood) > 0){
  #
  #       # Select first available neighbor
  #       k <- neighborhood[1]
  #
  #       # Perform swap
  #       ring[j, k] <- ring[k, j] <- FALSE
  #       ring[i, j] <- ring[j, i] <- TRUE
  #
  #       # Update budgets
  #       budget[i] <- budget[i] - 1
  #       budget[k] <- budget[k] + 1
  #
  #       # Total budget is unchanged with swap
  #       swapped <- TRUE
  #
  #       # Break out of loop
  #       break
  #
  #     }
  #
  #   }
  #
  #   # If completely struck, just break
  #   if(!swapped){
  #     break
  #   }
  #
  # }
  #
  # # Return ring and remaining
  # return(list(ring = ring, remaining = total_budget))

  # Coerce output ring back to integer
  storage.mode(ring) <- "integer"

  # Collect result
  result <- .Call(
    "swapping_pass_c",
    as.integer(nodes),
    as.integer(ring),
    as.integer(budget),
    as.integer(total_budget),
    PACKAGE = "psychNetwoRk"
  )

  # Set up ring matrix
  result$ring <- matrix(result$ring, nrow = nodes, ncol = nodes)

  # Return result
  return(result)

}

# # Interleaving function ----
# # Updated 24.03.2026
# interleave <- function(index, nodes, distance_sequence)
# {
#
#   return(
#     unique(
#       c(
#         rbind(
#           (index - 1 + distance_sequence) %% nodes + 1, # clockwise
#           (index - 1 - distance_sequence + nodes * distance_sequence) %% nodes + 1
#           # counterclockwise
#         )
#       )
#     )
#   )
#
# }

#' @noRd
# Assign weights based on distance ----
# Updated 26.03.2026
assign_weights <- function(network, A, distance_matrix)
{

  # Obtain weights
  lower_triangle <- lower.tri(network)
  network_nonzero <- network[lower_triangle] != 0
  weights <- network[lower_triangle][network_nonzero][
    order(abs(network)[lower_triangle][network_nonzero], decreasing = TRUE)
  ]

  # Set weights
  A_nonzero <- A[lower_triangle] != 0

  # Set weights order
  # Follows: Muldoon, Bridgeford, & Bassett's (2016) implementation
  weight_order <- rank(distance_matrix[lower_triangle][A_nonzero], ties.method = "random")
  A[lower_triangle][A_nonzero] <- weights[weight_order]
  A[!lower_triangle] <- 0
  A <- A + t(A) # make symmetric

  # Return weighted lattice
  return(A)

}
