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
#' Values greater than 1 indicate small-world structure
#'
#' \item \code{"S"} --- Computes the Humphries & Gurney (2008) \eqn{S}
#' metric by simulation. Random graphs are generated under the
#' Erdős–Rényi model (\code{\link[igraph]{sample_gnm}}) with the same
#' number of nodes and edges as \code{network}. The empirical and random
#' clustering coefficients both use global transitivity (\eqn{C^\Delta}).
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
#' realisable value. Random graphs preserve the degree sequence via
#' edge-switching (\code{\link[igraph]{sample_degseq}},
#' \code{method = "edge.switching.simple"}). The clustering coefficient uses
#' average local transitivity (\eqn{\bar{C}}). When \code{weighted = TRUE},
#' edge weights are reassigned to each random graph topology and the lattice
#' is constructed with \code{weighted = TRUE}. Note: SWI and SWP were
#' developed independently and concurrently; both apply the same
#' double-normalisation framework but differ in how the two deviation terms
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
#' via \code{assign_weights}, following Muldoon et al. (2016). When
#' \code{FALSE} (default), all graphs are treated as binary. Ignored when
#' \code{method = "analytical"}
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
# Updated 25.03.2026
smallworldness <- function(
    network, lattice, method = c("analytical", "omega", "S", "SWI", "SWP"),
    weighted = FALSE, iter = 100
)
{

  # Ensure network is in absolute values
  network <- abs(network)

  # Check for missing arguments (argument, default, function)
  # Uses actual function they will be used in
  method <- set_default(method, "swp", smallworldness)

  # Obtain adjacency matrix
  A <- network != 0

  # Check for weights
  if(!weighted){
    network <- A
  }

  # Obtain edges
  degree <- colSums(A, na.rm = TRUE)

  # Get nodes and degree
  nodes  <- length(degree)

  # Set up distance matrix
  node_sequence <- seq_len(nodes)
  distance_matrix <- abs(outer(node_sequence, node_sequence, "-"))
  distance_matrix <- pmin(distance_matrix, nodes - distance_matrix)

  # Convert network to {igraph}
  I <- convert2igraph(abs(network))

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
  ASPL <- igraph::mean_distance(I)
  CC <- igraph::transitivity(I, type = swiftelse(method %in% c("s", "analytical"), "global", "average"))
  # Switch CC based on method

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
# Analytical ----
# Updated 22.03.2026
smallworldness_analytical <- function(network, ASPL, CC, degree, ...)
{

  # Obtain nodes and average degree
  nodes <- dim(network)[2]
  mean_degree <- mean(degree, na.rm = TRUE)

  # Return analytical
  return(
    (CC / (mean_degree / nodes)) / (ASPL / (log(nodes) / log(mean_degree)))
  )

}

#' @noRd
# Telesford et al. (2011) ----
# Updated 25.03.2026
smallworldness_omega <- function(
    network, I, ASPL, CC, degree, lattice, weighted, distance_matrix, iter, ...
)
{

  # Collect random graphs
  random_graphs <- lapply(
    seq_len(iter), function(i){
      igraph::sample_degseq(out.deg = degree, method = "edge.switching.simple")
    }
  )

  # Check for weighted
  if(weighted){

    # Obtain weighted graphs
    random_graphs <- lapply(random_graphs, function(x){

      # Get weighted graphs
      weighted_graph <- assign_weights(network, igraph2matrix(x), distance_matrix)

      # Return back as {igraph}
      return(convert2igraph(weighted_graph))

    })

  }

  # Compute ASPL
  random_ASPL <- mean(nvapply(random_graphs, igraph::mean_distance), na.rm = TRUE)

  # Check if lattice is missing
  if(missing(lattice)){
    lattice_CC <- attributes(proxswap_lattice(network, weighted = weighted))$CC
  }else{

    # Check for weighted
    if(weighted){
      lattice <- assign_weights(network, lattice, distance_matrix)
    }

    lattice_CC <- igraph::transitivity(convert2igraph(lattice), type = "average")
  }

  # Return omega
  return((random_ASPL / ASPL) - (CC / lattice_CC))

}

#' @noRd
# Humphries & Gurney (2008) ----
# Updated 25.03.2026
smallworldness_S <- function(network, I, ASPL, CC, weighted, distance_matrix, iter, ...)
{

  # Obtain nodes and edges
  nodes <- dim(network)[2]
  edges <- sum(network[lower.tri(network)] != 0)

  # Obtain random graphs
  random_graphs <- lapply(seq_len(iter), function(i){
    igraph::sample_gnm(n = nodes, m = edges)
  })

  # Check for weighted
  if(weighted){

    # Obtain weighted graphs
    random_graphs <- lapply(random_graphs, function(x){

      # Get weighted graphs
      weighted_graph <- assign_weights(network, igraph2matrix(x), distance_matrix)

      # Return back as {igraph}
      return(convert2igraph(weighted_graph))

    })

  }

  # Compute S
  S <- nvapply(
    random_graphs, function(x){
      (CC / igraph::transitivity(x, type = "global")) / (ASPL / igraph::mean_distance(x))
    }
  )

  # Return S
  return(mean(S, na.rm = TRUE))

}

#' @noRd
# Neal (2015) ----
# Updated 25.03.2026
smallworldness_SWI <- function(
    network, I, ASPL, CC, degree, lattice, weighted, distance_matrix, iter, ...
)
{

  # Collect random graphs
  random_graphs <- lapply(
    seq_len(iter), function(i){
      igraph::sample_degseq(out.deg = degree, method = "edge.switching.simple")
    }
  )

  # Check for weighted
  if(weighted){

    # Obtain weighted graphs
    random_graphs <- lapply(random_graphs, function(x){

      # Get weighted graphs
      weighted_graph <- assign_weights(network, igraph2matrix(x), distance_matrix)

      # Return back as {igraph}
      return(convert2igraph(weighted_graph))

    })

  }

  # Compute ASPL and CC
  random_ASPL <- mean(nvapply(random_graphs, igraph::mean_distance), na.rm = TRUE)
  random_CC <- mean(nvapply(random_graphs, igraph::transitivity, type = "average"), na.rm = TRUE)

  # Check if lattice is missing
  if(missing(lattice)){
    lattice <- proxswap_lattice(network, weighted = weighted)
  }else if(weighted){
    lattice <- assign_weights(network, lattice, distance_matrix)
  }

  # Ensure lattice is {igraph}
  lattice <- convert2igraph(lattice)

  # Compute ASPL and CC
  lattice_ASPL <- igraph::mean_distance(lattice)
  lattice_CC <- igraph::transitivity(lattice, type = "average")

  # Return SWI
  return(
    ((ASPL - lattice_ASPL) / (random_ASPL - lattice_ASPL)) *
    ((CC - random_CC) / (lattice_CC - random_CC))
  )

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
# Updated 25.03.2026
smallworldness_SWP <- function(
    network, I, ASPL, CC, degree, lattice, weighted, distance_matrix, iter, ...
)
{

  # Collect random graphs
  random_graphs <- lapply(
    seq_len(iter), function(i){
      igraph::sample_degseq(out.deg = degree, method = "edge.switching.simple")
    }
  )

  # Check for weighted
  if(weighted){

    # Obtain weighted graphs
    random_graphs <- lapply(random_graphs, function(x){

      # Get weighted graphs
      weighted_graph <- assign_weights(network, igraph2matrix(x), distance_matrix)

      # Return back as {igraph}
      return(convert2igraph(weighted_graph))

    })

  }

  # Compute ASPL and CC
  random_ASPL <- mean(nvapply(random_graphs, igraph::mean_distance), na.rm = TRUE)
  random_CC <- mean(nvapply(random_graphs, igraph::transitivity, type = "average"), na.rm = TRUE)

  # Check if lattice is missing
  if(missing(lattice)){
    lattice <- proxswap_lattice(network, weighted = weighted)
  }else if(weighted){
    lattice <- assign_weights(network, lattice, distance_matrix)
  }

  # Ensure lattice is {igraph}
  lattice <- convert2igraph(lattice)

  # Compute ASPL and CC
  lattice_ASPL <- igraph::mean_distance(lattice)
  lattice_CC <- igraph::transitivity(lattice, type = "average")

  # Compute deltas
  delta_ASPL <- max(0, (ASPL - random_ASPL)) / (lattice_ASPL - random_ASPL)
  delta_CC <- max(0, (lattice_CC - CC)) / (lattice_CC - random_CC)

  # Boundary check
  delta_ASPL <- min(1, swiftelse(boundary_values(ASPL, random_ASPL, lattice_ASPL), 1, delta_ASPL))
  delta_CC <- min(1, swiftelse(boundary_values(CC, random_CC, lattice_CC), 1, delta_CC))

  # Return SWP
  return(1 - sqrt((delta_ASPL^2 + delta_CC^2) / 2))

}