#' @title Computes Various Small-Worldness Metrics
#'
#' @description Computes the small-worldness of a network using one of five
#' methods. All simulation-based methods generate degree-preserving random
#' graphs as the null baseline. The lattice reference network (where required)
#' is constructed using \code{\link[psychNetwoRk]{proxswap_lattice}}, which produces a
#' degree-preserving ring lattice that maximizes the clustering coefficient
#'
#' @param network Matrix or data frame.
#' A square, symmetric numeric matrix representing a network (e.g., partial
#' correlations). Absolute values are taken internally, so signed weights are
#' handled automatically. Non-zero off-diagonal entries are treated as edges.
#'
#' @param lattice Matrix (optional).
#' A pre-computed lattice adjacency matrix to use as the regular-network
#' reference for the \code{"omega"} and \code{"SWP"} methods.
#' If not provided, a lattice is generated automatically via
#' \code{\link[psychNetwoRk]{proxswap_lattice}}. Ignored by \code{"analytical"}
#' and \code{"S"}
#'
#' @param method Character (length = 1).
#' The method used to compute small-worldness.
#' Defaults to \code{"SWP"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"analytical"} --- Computes the Humphries & Gurney (2008)
#' \eqn{S} metric using closed-form approximations for the Erdős–Rényi
#' random graph baseline:
#' \deqn{S = \frac{C / C_{\text{rand}}}{L / L_{\text{rand}}}}
#' where \eqn{C_{\text{rand}} \approx \langle k \rangle / n} and
#' \eqn{L_{\text{rand}} \approx \ln(n) / \ln(\langle k \rangle)}.
#' No random graphs are generated. The \code{weighted} argument is ignored.
#' The approximation requires an average degree greater than 1 (otherwise
#' \code{NA} is returned with a warning).
#' Values greater than 1 indicate small-world structure
#'
#' \item \code{"S"} --- Computes the Humphries & Gurney (2008) \eqn{S}
#' metric by simulation. Random graphs are generated under the
#' Erdős–Rényi model (\code{\link[igraph]{sample_gnm}}) with the same
#' number of nodes and edges as \code{network}. The empirical and random
#' clustering coefficients both use global transitivity (\eqn{C^\Delta}).
#' \eqn{S} is computed from the average clustering coefficient and path length
#' across the random graphs.
#' When \code{weighted = TRUE}, edge weights are reassigned to each random
#' graph topology via \code{assign_weights}. Values greater than 1
#' indicate small-world structure
#'
#' \item \code{"omega"} --- Computes the \eqn{\omega} metric of
#' Telesford et al. (2011):
#' \deqn{\omega = \frac{L_{\text{rand}}}{L} - \frac{C}{C_{\text{latt}}}}
#' Random graphs preserve the degree sequence via edge-switching
#' (\code{\link[igraph]{sample_degseq}}, \code{method = "edge.switching.simple"}).
#' The clustering coefficient uses average local transitivity (\eqn{\bar{C}}).
#' Note: Telesford et al. (2011) defined and validated \eqn{\omega} exclusively
#' on binary, unweighted networks. When \code{weighted = TRUE}, edge weights
#' are reassigned to each random graph topology via \code{assign_weights} and
#' the lattice is constructed with \code{weighted = TRUE}, following the
#' approach of Muldoon et al. (2016); this is an extension beyond the original
#' formulation. Values near zero indicate small-world structure; negative
#' values indicate lattice-like structure; positive values indicate
#' random-like structure. Bounded in \eqn{[-1, 1]}
#'
#' \item \code{"SWI"} --- Computes the Small-World Index of
#' Neal (2015):
#' \deqn{\text{SWI} = \frac{L - L_{\text{latt}}}{L_{\text{rand}} - L_{\text{latt}}} \times \frac{C - C_{\text{rand}}}{C_{\text{latt}} - C_{\text{rand}}}}
#' Each term captures where the observed network's path length and clustering
#' coefficient fall within the lattice-to-random range, so that SWI = 1
#' only when \eqn{L = L_{\text{rand}}} and \eqn{C = C_{\text{latt}}}
#' simultaneously. Because this ideal is mathematically unachievable in a
#' finite network, SWI = 1 is a conceptual upper bound rather than a
#' realizable value. Random graphs preserve the degree sequence via
#' edge-switching (\code{\link[igraph]{sample_degseq}},
#' \code{method = "edge.switching.simple"}). The clustering coefficient uses
#' average local transitivity (\eqn{\bar{C}}). When \code{weighted = TRUE},
#' edge weights are reassigned to each random graph topology and the lattice
#' is constructed with \code{weighted = TRUE}. Values are bounded to
#' \eqn{[0, 1]}. Note: SWI and SWP were
#' developed independently and concurrently; both apply the same
#' double-normalization framework but differ in how the two deviation terms
#' are combined (product vs. Euclidean distance). Values close to 1 indicate
#' strong small-world structure; values near 0 indicate lattice-like or
#' random-like structure. Bounded in \eqn{[0, 1]}
#'
#' \item \code{"SWP"} --- Computes the Small-World Propensity of
#' Muldoon et al. (2016):
#' \deqn{\phi = 1 - \sqrt{\frac{\Delta_C^2 + \Delta_L^2}{2}}}
#' where \eqn{\Delta_C = (C_{\text{latt}} - C) / (C_{\text{latt}} - C_{\text{rand}})}
#' and \eqn{\Delta_L = (L - L_{\text{rand}}) / (L_{\text{latt}} - L_{\text{rand}})}.
#' Both \eqn{\Delta_C} and \eqn{\Delta_L} are bounded to \eqn{[0, 1]},
#' guaranteeing \eqn{\phi \in [0, 1]}. Random graphs preserve the degree
#' sequence via edge-switching. The clustering coefficient uses average
#' local transitivity (\eqn{\bar{C}}). When \code{weighted = TRUE}, edge
#' weights are reassigned to each random graph topology via
#' \code{assign_weights} and the lattice is constructed with
#' \code{weighted = TRUE}. Values close to 1 indicate strong small-world
#' structure; a pragmatic threshold of \eqn{\phi_T = 0.60} has been suggested
#'
#' }
#'
#' @param weighted Logical (length = 1).
#' Whether to compute small-worldness on the weighted network. When
#' \code{TRUE}, edge weights from \code{network} are preserved for the
#' empirical graph and reassigned to random and lattice graph topologies
#' via \code{assign_weights}, following Muldoon et al. (2016): path lengths
#' use distances equal to the inverse of the (absolute) weights and clustering
#' uses the weighted clustering coefficient of Onnela et al. (2005) with
#' weights normalized by the largest weight. When \code{FALSE} (default), all
#' graphs are treated as binary. Ignored when \code{method = "analytical"}
#'
#' @param seed Numeric (length = 1).
#' Seed for reproducible results. The same \code{seed} gives the same results
#' (and does not change the random number stream of the R session).
#' Defaults to \code{NULL} (random seed)
#'
#' @param iter Numeric (length = 1).
#' Number of random graphs to generate when estimating the null baseline.
#' Defaults to \code{100}.
#' Ignored when \code{method = "analytical"}.
#' Higher values produce more stable estimates at the cost of computation time
#'
#' @return Numeric (length = 1).
#' The small-worldness value computed using the chosen \code{method}:
#'
#' \itemize{
#'
#' \item \code{"analytical"} and \code{"S"} --- \eqn{S > 1} indicates small-world structure
#'
#' \item \code{"omega"} --- Values near \eqn{|\omega| < 0.50} indicate small-world structure
#' (ranges between -1 and 1)
#'
#' \item \code{"SWI"} --- Values close to 1 indicate strong small-world structure;
#' values near 0 indicate lattice-like or random-like structure
#' (ranges between 0 and 1)
#'
#' \item \code{"SWP"} --- Values \eqn{\phi > 0.60} indicate strong small-world structure
#' (ranges between 0 and 1)
#'
#' }
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{omega} \cr
#' Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011).
#' The ubiquity of small-world networks.
#' \emph{Brain Connectivity}, \emph{1}(5), 367--375.
#'
#' \strong{S and analytical} \cr
#' Humphries, M. D., & Gurney, K. (2008).
#' Network 'small-world-ness': A quantitative method for determining canonical network equivalence.
#' \emph{PLoS ONE}, \emph{3}(4), e0002051.
#'
#' \strong{SWI} \cr
#' Neal, Z. P. (2015).
#' Making big communities small: Using network science to understand the ecological and behavioral requirements for community social capital.
#' \emph{American Journal of Community Psychology}, \emph{55}(3), 369--380.
#'
#' \strong{SWP} \cr
#' Muldoon, S. F., Bridgeford, E. W., & Bassett, D. S. (2016).
#' Small-world propensity and weighted brain networks.
#' \emph{Scientific Reports}, \emph{6}(1), 22057.
#'
#' \strong{Weighted clustering coefficient} \cr
#' Onnela, J.-P., Saramaki, J., Kertesz, J., & Kaski, K. (2005).
#' Intensity and coherence of motifs in weighted complex networks.
#' \emph{Physical Review E}, \emph{71}(6), 065103.
#'
#' @examples
#' # Simulate a small-world network (population partial correlations)
#' network <- simulate_smallworld(
#'   nodes = 20, density = 0.30, rewire = 0.20, sample_size = 500
#' )$population$Omega
#'
#' # Compute SWP (default)
#' swp <- smallworldness(network)
#'
#' # Compute omega
#' omega <- smallworldness(network, method = "omega")
#'
#' # Compute analytical S
#' S <- smallworldness(network, method = "analytical")
#'
#' # Compute simulated S
#' S_sim <- smallworldness(network, method = "S")
#'
#' # Compute SWI
#' swi <- smallworldness(network, method = "SWI")
#'
#' # Compute weighted SWP
#' swp_w <- smallworldness(network, weighted = TRUE)
#'
#' # Compute weighted omega
#' omega_w <- smallworldness(network, method = "omega", weighted = TRUE)
#'
#' @export
#'
# Compute smallworldness ----
# Updated 20.11.2026
smallworldness <- function(
    network, lattice, method = c("analytical", "omega", "S", "SWI", "SWP"),
    weighted = FALSE, iter = 100, seed = NULL
)
{

  # Set seed (reproducible results)
  seed_token <- seed_start(seed)
  on.exit(seed_end(seed_token), add = TRUE)

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  method <- set_default(method, "swp", smallworldness)

  # Check inputs
  network <- smallworldness_errors(network, weighted, iter)

  # Ensure network is in absolute values with a zero diagonal
  network <- abs(network)
  diag(network) <- 0

  # Obtain adjacency matrix
  A <- network != 0

  # Check for weights
  if(!weighted){
    network <- A + 0
  }

  # Obtain degree and nodes
  degree <- colSums(A)
  nodes  <- length(degree)

  # Set up distance matrix
  distance_matrix <- ring_distance(nodes)

  # Convert network to {igraph}
  I <- convert2igraph(network)

  # Set up function switch
  smallworld_FUN <- switch(
    method,
    "analytical" = smallworldness_analytical,
    "omega" = smallworldness_omega,
    "s" = smallworldness_S,
    "swi" = smallworldness_SWI,
    "swp" = smallworldness_SWP
  )

  # Obtain empirical values
  ASPL <- graph_path_length(I, weighted && method != "analytical")
  CC <- graph_clustering(
    network, swiftelse(method %in% c("s", "analytical"), "global", "average"),
    weighted && method != "analytical"
  )

  # Compute smallworldness
  return(
    smallworld_FUN(
      network = network, I = I, ASPL = ASPL, CC = CC, degree = degree,
      lattice = lattice, weighted = weighted, distance_matrix = distance_matrix,
      iter = iter
    )
  )

}

#' @noRd
# Errors ----
# Updated 20.11.2026
smallworldness_errors <- function(network, weighted, iter)
{

  # Check network
  object_error(network, c("matrix", "data.frame", "tibble"), "smallworldness")
  network <- as.matrix(network)
  typeof_error(network, "numeric", "smallworldness")
  if(!is_symmetric(network)){
    stop("`network` must be a square, symmetric matrix", call. = FALSE)
  }

  # Check weighted
  typeof_error(weighted, "logical", "smallworldness")
  length_error(weighted, 1, "smallworldness")

  # Check iterations
  typeof_error(iter, "numeric", "smallworldness")
  length_error(iter, 1, "smallworldness")
  range_error(iter, c(1, Inf), "smallworldness")

  # Return network
  return(network)

}

#' @noRd
# Average shortest path length ----
# Weighted path lengths use distances equal to the inverse of the weights
# (stronger connections are shorter; Muldoon et al., 2016)
# Updated 20.11.2026
graph_path_length <- function(I, weighted)
{
  if(weighted){
    return(igraph::mean_distance(I, weights = 1 / igraph::E(I)$weight))
  }else{
    return(igraph::mean_distance(I, weights = NA))
  }
}

#' @noRd
# Clustering coefficient ----
# Binary: global transitivity or average local transitivity (nodes with
# degree less than 2 are excluded, as in {igraph}). Weighted: the weighted
# clustering coefficient of Onnela et al. (2005) with weights normalized by
# the largest weight (C_i = (W^(1/3))^3_ii / (k_i (k_i - 1))), which equals
# the binary coefficient for binary networks
# Updated 20.11.2026
graph_clustering <- function(network, type = c("global", "average"), weighted)
{

  # Binary (as {igraph})
  if(!weighted){
    return(igraph::transitivity(convert2igraph(network), type = type))
  }

  # Normalized weights (cube root)
  network <- as.matrix(network)
  diag(network) <- 0
  W <- (abs(network) / max(abs(network)))^(1 / 3)

  # Degree and closed walks of length three
  degree <- colSums(network != 0)
  triangles <- diag(W %*% W %*% W)
  possible <- degree * (degree - 1)

  # Return clustering
  if(type == "global"){
    return(sum(triangles) / sum(possible))
  }else{
    return(mean(triangles[degree > 1] / possible[degree > 1]))
  }

}

#' @noRd
# Random reference graphs ----
# Degree-preserving (edge switching) or Erdos-Renyi graphs with the same
# number of edges; when weighted, the network's weights are reassigned
# to each random topology (`assign_weights`)
# Updated 20.11.2026
random_reference_graphs <- function(
    network, degree, iter, weighted, distance_matrix, type = c("degree", "gnm")
)
{

  # Generate random graphs
  graphs <- lapply(seq_len(iter), function(i){
    if(type == "degree"){
      igraph::sample_degseq(out.deg = degree, method = "edge.switching.simple")
    }else{
      igraph::sample_gnm(n = length(degree), m = sum(degree) / 2)
    }
  })

  # Return adjacency or weighted matrices
  return(
    lapply(graphs, function(x){
      A <- igraph2matrix(x)
      swiftelse(weighted, list(assign_weights(network, A, distance_matrix)), list(A))[[1]]
    })
  )

}

#' @noRd
# Path length and clustering of random reference graphs ----
# Updated 20.11.2026
random_reference_values <- function(graphs, weighted, type = "average")
{
  return(
    list(
      ASPL = mean(nvapply(graphs, function(x){
        graph_path_length(convert2igraph(x), weighted)
      }), na.rm = TRUE),
      CC = mean(nvapply(graphs, graph_clustering, type = type, weighted = weighted), na.rm = TRUE)
    )
  )
}

#' @noRd
# Lattice reference ----
# Supplied lattice (weighted with `assign_weights` when weighted) or
# a degree-preserving lattice from `proxswap_lattice`
# Updated 20.11.2026
lattice_reference <- function(network, lattice, weighted, distance_matrix)
{

  # Obtain lattice
  if(missing(lattice)){
    lattice <- proxswap_lattice(network, weighted = weighted)
  }else if(weighted){
    lattice <- assign_weights(network, lattice, distance_matrix)
  }

  # Return path length and clustering
  return(
    list(
      ASPL = graph_path_length(convert2igraph(lattice), weighted),
      CC = graph_clustering(lattice, "average", weighted)
    )
  )

}

#' @noRd
# Analytical ----
# Updated 20.11.2026
smallworldness_analytical <- function(network, ASPL, CC, degree, ...)
{

  # Obtain nodes and average degree
  nodes <- dim(network)[2]
  mean_degree <- mean(degree, na.rm = TRUE)

  # Approximation requires an average degree greater than one
  if(mean_degree <= 1){
    warning(
      "The analytical approximation requires an average degree greater than 1; `NA` was returned",
      call. = FALSE
    )
    return(NA_real_)
  }

  # Return analytical
  return(
    (CC / (mean_degree / nodes)) / (ASPL / (log(nodes) / log(mean_degree)))
  )

}

#' @noRd
# Telesford et al. (2011) ----
# Updated 20.11.2026
smallworldness_omega <- function(
    network, I, ASPL, CC, degree, lattice, weighted, distance_matrix, iter, ...
)
{

  # Random reference
  random <- random_reference_values(
    random_reference_graphs(network, degree, iter, weighted, distance_matrix, "degree"),
    weighted
  )

  # Lattice reference
  regular <- lattice_reference(network, lattice, weighted, distance_matrix)

  # Return omega
  return((random$ASPL / ASPL) - (CC / regular$CC))

}

#' @noRd
# Humphries & Gurney (2008) ----
# Updated 20.11.2026
smallworldness_S <- function(network, I, ASPL, CC, degree, weighted, distance_matrix, iter, ...)
{

  # Random reference (Erdos-Renyi with the same number of edges)
  random <- random_reference_values(
    random_reference_graphs(network, degree, iter, weighted, distance_matrix, "gnm"),
    weighted, type = "global"
  )

  # Return S
  return((CC / random$CC) / (ASPL / random$ASPL))

}

#' @noRd
# Neal (2015) ----
# Updated 20.11.2026
smallworldness_SWI <- function(
    network, I, ASPL, CC, degree, lattice, weighted, distance_matrix, iter, ...
)
{

  # Random reference
  random <- random_reference_values(
    random_reference_graphs(network, degree, iter, weighted, distance_matrix, "degree"),
    weighted
  )

  # Lattice reference
  regular <- lattice_reference(network, lattice, weighted, distance_matrix)

  # Compute SWI
  SWI <- ((ASPL - regular$ASPL) / (random$ASPL - regular$ASPL)) *
    ((CC - random$CC) / (regular$CC - random$CC))

  # Return SWI (bounded to [0, 1])
  return(swiftelse(is.finite(SWI), min(1, max(0, SWI)), NA_real_))

}

#' @noRd
# Boundaries for SWP ----
# Updated 22.03.2026
boundary_values <- function(empirical, random, regular)
{
  any(is.infinite(c(empirical, random, regular))) | any(is.na(c(empirical, random, regular)))
}

#' @noRd
# Muldoon et al. (2016) ----
# Updated 20.11.2026
smallworldness_SWP <- function(
    network, I, ASPL, CC, degree, lattice, weighted, distance_matrix, iter, ...
)
{

  # Random reference
  random <- random_reference_values(
    random_reference_graphs(network, degree, iter, weighted, distance_matrix, "degree"),
    weighted
  )

  # Lattice reference
  regular <- lattice_reference(network, lattice, weighted, distance_matrix)

  # Compute deltas (bounded to [0, 1]; degenerate references give 1)
  delta_ASPL <- (ASPL - random$ASPL) / (regular$ASPL - random$ASPL)
  delta_CC <- (regular$CC - CC) / (regular$CC - random$CC)
  delta_ASPL <- swiftelse(
    boundary_values(ASPL, random$ASPL, regular$ASPL) || !is.finite(delta_ASPL),
    1, min(1, max(0, delta_ASPL))
  )
  delta_CC <- swiftelse(
    boundary_values(CC, random$CC, regular$CC) || !is.finite(delta_CC),
    1, min(1, max(0, delta_CC))
  )

  # Return SWP
  return(1 - sqrt((delta_ASPL^2 + delta_CC^2) / 2))

}
