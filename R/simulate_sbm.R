#' Simulates Stochastic Block Model Data
#'
#' @description
#' Simulates data from a Gaussian Graphical Model (GGM) with a stochastic
#' block model (SBM) structure. Nodes are partitioned into communities, with
#' edge density controlled separately within and between communities via a
#' \code{blocks x blocks} density matrix. Absolute edge weights are drawn from
#' a Weibull distribution whose parameters are predicted from the network size,
#' sample size, and signal-to-noise ratio using a Seemingly Unrelated
#' Regression (SUR) model fitted to 194 empirical psychometric networks
#' (Huth et al., 2025). The resulting population partial correlation matrix is
#' used to generate multivariate normal (or skewed) data.
#'
#' A \code{diffusion} parameter controls the minimum proportion of the
#' strongest edges that are reserved for within-community positions. Because
#' the remaining edges are shuffled randomly, the effective within-community
#' advantage will generally exceed \code{1 - diffusion}; see Details.
#'
#' Parameters do not have default values (except \code{negative_proportion},
#' \code{diffusion}, \code{target_condition}, \code{max_correlation}, and
#' \code{max_iterations}) and must each be set. See Details and Examples to
#' get started.
#'
#' @param nodes Numeric (length = 1 or \code{blocks}).
#' Number of nodes per community block. Can be a single value applied to all
#' blocks, or a vector of length \code{blocks} specifying each block's size
#' individually. Minimum of three nodes per block. The total number of nodes
#' (\code{sum(nodes)}) should be between 8 and 54 to remain within the range
#' of the empirical networks used to fit the Weibull parameter model; values
#' outside this range are accepted but will trigger extrapolation and may
#' produce a warning from \code{\link{weibull_parameters}}.
#'
#' @param blocks Numeric (length = 1).
#' Number of community blocks.
#'
#' @param density_matrix Matrix (dim = \code{blocks x blocks}).
#' A symmetric numeric matrix specifying edge probabilities within and between
#' community blocks. Diagonal entry \eqn{[i,i]} gives the within-block edge
#' probability for community \eqn{i}; off-diagonal entry \eqn{[i,j]}
#' (\eqn{i \neq j}) gives the between-block edge probability for the pair of
#' communities \eqn{i} and \eqn{j}. All entries must be in \eqn{[0, 1]} and
#' the matrix must be symmetric (i.e., \code{density_matrix[i,j] ==
#' density_matrix[j,i]}). See Details for construction guidance.
#'
#' @param snr Numeric (length = 1).
#' Signal-to-noise ratio of the absolute partial correlation weights,
#' defined as \eqn{\bar{|w|} / \mathrm{SD}(|w|)}. This ratio governs the
#' shape and scale of the Weibull distribution used to generate edge weights:
#' values below 1 produce wider, more heterogeneous weight distributions;
#' values above 1 produce narrower, more homogeneous distributions. Note that
#' the same SNR can arise from different combinations of mean and standard
#' deviation (e.g., mean = 0.05, SD = 0.05 and mean = 0.15, SD = 0.15 both
#' give SNR = 1), so SNR alone does not determine the absolute magnitude of
#' edge weights, which is additionally governed by \code{nodes} and
#' \code{sample_size} via \code{\link{weibull_parameters}}. Empirically
#' observed SNR values ranged from 0.648 to 1.712; values outside this range
#' are accepted but will trigger a warning. Defaults to \code{1}.
#'
#' @param negative_proportion Numeric (length = 1).
#' Proportion of between-community edges assigned a negative sign (inhibitory
#' connections). Each between-community edge is independently signed negative
#' with this probability (i.e., a Bernoulli draw per edge), so the realized
#' proportion will vary around the specified value. Must be between 0 and 1.
#' Applies only to between-community edges; all within-community edges are
#' positive. If not provided, a value is sampled from a truncated normal
#' distribution reflecting the empirical distribution of true negative partial
#' correlations across 194 psychometric networks: mean = 0.34, SD = 0.086,
#' bounded to \eqn{[0.083, 0.55]}.
#'
#' @param diffusion Numeric (length = 1).
#' Controls the minimum proportion of the top-ranked edges (by absolute
#' weight) that are guaranteed to be placed within communities rather than
#' distributed freely across the network. Specifically, \code{1 - diffusion}
#' of the within-community edges are reserved for the highest-weight draws;
#' the remaining edges (both within- and between-community) are filled from
#' a randomly shuffled pool. Because shuffled edges can land within communities
#' by chance, the actual proportion of top edges that end up within communities
#' will on average exceed \code{1 - diffusion}. Lower values of
#' \code{diffusion} (e.g., \code{0.10}) therefore produce stronger community
#' contrast than higher values (e.g., \code{0.90}), but \code{diffusion}
#' should be interpreted as a floor on within-community weight concentration,
#' not an exact control. Defaults to \code{0.30}. Must be between 0 and 1.
#'
#' @param diffusion_range Numeric (length = 2).
#' If provided, overrides \code{diffusion} by drawing the diffusion proportion
#' uniformly from this interval on each call. Useful for introducing
#' replication-to-replication variability in community contrast. For example,
#' \code{diffusion_range = c(0.05, 0.20)} samples a value between 5\% and
#' 20\% on each draw. The same floor interpretation applies as for
#' \code{diffusion}. Both values must be between 0 and 1.
#'
#' @param sample_size Numeric (length = 1).
#' Number of observations to generate from the population multivariate
#' distribution. Also influences the predicted Weibull scale parameter via
#' \code{\link{weibull_parameters}}: larger samples are associated with
#' smaller, more precisely estimated edge weights.
#'
#' @param skew Numeric (length = 1 or \code{sum(nodes)}).
#' Skew applied to each variable after generation from the multivariate
#' normal. Can be a single value applied to all variables, or one value per
#' variable. Values are rounded to the nearest 0.05 increment and must be
#' in \eqn{[-2, 2]}. Defaults to \code{0} (no skew).
#'
#' @param skew_range Numeric (length = 2).
#' If provided, overrides \code{skew} by drawing a skew value independently
#' and uniformly from this interval for each variable. Both values must be
#' in \eqn{[-2, 2]}.
#'
#' @param target_condition Numeric (length = 1).
#' Target condition number (computed via \code{\link{kappa}} with
#' \code{exact = TRUE}) used when ridge regularization is required to recover
#' a positive definite precision matrix. The smallest ridge penalty
#' \eqn{\lambda} that brings the condition number to this target is found via
#' root-finding (\code{\link{uniroot}}), subject to a maximum shrinkage of
#' approximately 23\% following Peeters et al. (2020). After conditioning,
#' the Weibull bounds are re-verified on the updated edge weights; draws that
#' fall outside the empirical bounds after conditioning are rejected. Lower
#' values produce better-conditioned (more stable) matrices. Defaults to
#' \code{30}. Values up to \code{100} are accepted but not recommended.
#'
#' @param max_correlation Numeric (length = 1).
#' Maximum allowed absolute pairwise correlation in the population correlation
#' matrix \code{R}. Any draw where \code{max(abs(R[lower.tri(R)])) >
#' max_correlation} is rejected and a new attempt is made. Must be between
#' 0 and 1. Defaults to \code{0.80}.
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of attempts to find a connected network with valid edge
#' weights before stopping with an error. The error message reports a
#' frequency table of rejection reasons to assist with diagnosing convergence
#' failures. Defaults to \code{100}.
#'
#' @details
#' \strong{Constructing \code{density_matrix}}
#'
#' The \code{density_matrix} is a \code{blocks x blocks} symmetric matrix
#' where entry \eqn{[i,i]} gives the within-block edge probability for
#' community \eqn{i}, and entry \eqn{[i,j]} (for \eqn{i \neq j}) gives the
#' between-block edge probability for communities \eqn{i} and \eqn{j}. The
#' simplest construction uses a uniform off-diagonal density with
#' block-specific diagonals:
#'
#' \preformatted{
#' # Uniform within (0.90) and between (0.20) density for 3 blocks
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- 0.90
#' }
#'
#' For asymmetric community structure, each diagonal entry can differ:
#'
#' \preformatted{
#' # Varying within-block density per community
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- c(0.85, 0.90, 0.95)
#' }
#'
#' For full pairwise control over between-block densities, specify the
#' complete symmetric matrix directly:
#'
#' \preformatted{
#' dm <- matrix(c(
#'   0.90, 0.20, 0.10,
#'   0.20, 0.85, 0.30,
#'   0.10, 0.30, 0.95
#' ), nrow = 3, ncol = 3)
#' }
#'
#' \strong{Diffusion and within-community weight concentration}
#'
#' The \code{diffusion} parameter does not exactly fix the proportion of
#' top-ranked edges placed within communities. Instead, \code{1 - diffusion}
#' of the within-community edge slots are filled deterministically from the
#' highest-weight draws. The remaining edge slots (within- and
#' between-community alike) are filled from a randomly shuffled pool of
#' lower-ranked weights, meaning some additional high-weight edges will land
#' within communities by chance. The realized within-community weight
#' advantage will therefore always be at least as large as implied by
#' \code{1 - diffusion}, and typically larger. The \code{Q} field in the
#' returned \code{parameters} list (Newman-Girvan modularity) provides a
#' post-hoc summary of the actual community contrast achieved.
#'
#' @return A named list of class \code{"ggm_simulate"} with four elements
#' (see \code{\link{add_population_error}} to add population error):
#'
#' \item{data}{Numeric matrix of dimension \code{sample_size x sum(nodes)}
#' containing the simulated observations drawn from the population GGM.
#' Rows are cases; columns are variables named \code{V01}, \code{V02},
#' etc. Values are continuous (or skewed continuous when \code{skew != 0}).
#' To produce ordinal data, pass the columns through
#' \code{\link{categorize}}.}
#'
#' \item{parameters}{
#' A list of input, derived, and estimated parameters:
#' \itemize{
#'   \item \code{nodes} --- Integer vector of length \code{blocks} giving the
#'   number of nodes per block (scalar input is expanded to this length)
#'   \item \code{blocks} --- Number of community blocks
#'   \item \code{sample_size} --- Number of simulated observations
#'   \item \code{skew} --- Named numeric vector of per-variable skew values
#'   actually applied (after rounding and possible resampling)
#'   \item \code{density_matrix} --- The \code{blocks x blocks} density matrix
#'   as supplied
#'   \item \code{negative_proportion} --- The proportion of between-community
#'   edges assigned a negative sign, either as supplied or as sampled from the
#'   empirical distribution
#'   \item \code{weibull} --- Named numeric vector of length 2 giving the
#'   Weibull \code{shape} and \code{scale} parameters of the absolute edge
#'   weight distribution actually used. If ridge conditioning was applied,
#'   these are re-estimated from the conditioned network via MLE.
#'   \item \code{diffusion} --- Numeric vector of length 2 giving the
#'   within-block reservation range as a proportion of within-community edges,
#'   on the internal scale used by the sampler (i.e., \code{range(1 -
#'   diffusion)} or \code{range(1 - diffusion_range)}). Both values are equal
#'   when \code{diffusion} is scalar. Note this is the complement of the
#'   user-supplied \code{diffusion} value and represents the fraction of
#'   within-community slots filled from the top-ranked draws.
#'   \item \code{omega} --- Smallworldness omega statistic of the generated
#'   network (Telesford et al., 2011). Values near zero indicate small-world
#'   structure; negative values indicate lattice-like structure; positive
#'   values indicate random-like structure (see \code{\link[psychNetwoRk]{smallworldness}})
#'   \item \code{Q} --- Newman-Girvan modularity of the population network
#'   (\code{Omega}) with respect to the block membership, computed via
#'   \code{igraph::modularity} on absolute edge weights. Provides a summary
#'   of the community contrast actually achieved after weight assignment and
#'   any ridge conditioning.
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
#'   \item \code{membership} --- Named integer vector of length
#'   \code{sum(nodes)} giving the community block assignment (1 to
#'   \code{blocks}) for each node
#' }
#' }
#'
#' \item{convergence}{
#' Iteration and conditioning diagnostics:
#' \itemize{
#'   \item \code{iterations} --- Number of sampling attempts needed to find a
#'   valid network (including graph structure and edge weight draws)
#'   \item \code{rejections} --- Character vector of length
#'   \code{max_iterations + 1} recording the rejection reason for each failed
#'   attempt; entries for successful or unused iterations are empty strings.
#'   Common reasons include disconnected graph structure, condition number
#'   exceeding \code{target_condition}, maximum correlation exceeding
#'   \code{max_correlation}, and Weibull parameters falling outside empirical
#'   bounds after ridge conditioning.
#'   \item \code{lambda} --- Ridge regularization parameter \eqn{\lambda}
#'   added to the diagonal of the precision matrix to ensure positive
#'   definiteness; \code{NA} if no conditioning was required
#'   \item \code{condition} --- Condition number of the final population
#'   correlation matrix \code{R}, computed via \code{kappa} with
#'   \code{exact = TRUE}
#' }
#' }
#'
#' @examples
#' # Construct density matrix for 3 blocks with uniform densities
#' dm <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm) <- 0.90
#'
#' # Basic 3-block simulation with equal-sized communities
#' result <- simulate_sbm(
#'   nodes = 6, # 6 nodes per block = 18 total
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm
#' )
#' \donttest{
#' # Unequal block sizes
#' result <- simulate_sbm(
#'   nodes = c(4, 6, 8), # 18 total nodes
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm
#' )
#'
#' # Varying within-block density per community
#' dm_varying <- matrix(0.20, nrow = 3, ncol = 3)
#' diag(dm_varying) <- c(0.85, 0.90, 0.95)
#'
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm_varying
#' )
#'
#' # Full pairwise between-block density control
#' dm_pairwise <- matrix(c(
#'   0.90, 0.20, 0.10,
#'   0.20, 0.85, 0.30,
#'   0.10, 0.30, 0.95
#' ), nrow = 3, ncol = 3)
#'
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm_pairwise
#' )
#'
#' # Fix the proportion of negative between-community edges
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm,
#'   negative_proportion = 0.20
#' )
#'
#' # Introduce variability in diffusion across replications
#' result <- simulate_sbm(
#'   nodes = 6,
#'   blocks = 3,
#'   sample_size = 500,
#'   density_matrix = dm,
#'   diffusion_range = c(0.30, 0.70)
#' )
#' }
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @references
#' \strong{Simulation approach for stochastic block model networks} \cr
#' Christensen, A. P., & Choi, J. (2026).
#' Network accuracy across levels of analysis using stochastic block models.
#' \emph{Psychological Methods}.
#' \doi{10.1037/met0000853}
#'
#' \strong{Seminal introduction to Stochastic Block Models} \cr
#' Holland, P. W., Laskey, K. B., & Leinhardt, S. (1983).
#' Stochastic blockmodels: First steps.
#' \emph{Social Networks}, \emph{5}(2), 109--137.
#'
#' \strong{Empirical network data used to fit the Weibull SUR model} \cr
#' Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., &
#' Marsman, M. (2025). Statistical evidence in psychological networks.
#' \emph{Nature Human Behaviour}.
#' \doi{10.1038/s41562-025-02314-2}
#'
#' \strong{Maximum ridge shrinkage bound} \cr
#' Peeters, C. F., van de Wiel, M. A., & van Wieringen, W. N. (2020).
#' The spectral condition number plot for regularization parameter evaluation.
#' \emph{Computational Statistics}, \emph{35}(2), 629--646.
#'
#' @export
#'
# Simulate SBM GGM data ----
# Updated 15.05.2026
simulate_sbm <- function(
    nodes, blocks, density_matrix,
    snr = 1, diffusion = 0.30, diffusion_range = NULL,
    negative_proportion, sample_size,
    skew = 0, skew_range = NULL,
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
      ), 0.55 # empirical maximum
    )

  }

  # Check for input errors
  ## Returns updated nodes to ensure same length as blocks
  nodes <- simulate_sbm_errors(
    nodes, blocks, density_matrix,
    snr, negative_proportion, sample_size, skew, diffusion,
    target_condition, max_correlation, max_iterations
  )

  # Determine total number of nodes
  total_nodes <- sum(nodes)

  # Get community sequence
  community_sequence <- seq_len(blocks)

  # Set membership
  membership <- sort(rep(community_sequence, times = nodes))

  # Set up membership matrix
  membership_matrix <-  outer(membership, membership, "==")

  # Obtain lower triangle
  lower_triangle <- lower.tri(membership_matrix)

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
          "Reached maximum iterations: Could not find weights that converge for the SBM structure. ",
          "Resulting rejections were due to:\n\n",
          paste0(names(rejection_table), " = ", rejection_table, collapse = "\n")
        )
      )

    }

    # Get block matrix
    block_matrix <- generate_sbm(
      membership_matrix, blocks, membership,
      membership_matrix, density_matrix, lower_triangle
    )

    # Check block matrix
    if(!igraph::is_connected(convert2igraph(block_matrix))){

      # Add rejection reason
      rejections[iter] <- "Could not find structure where every node was connected."

      # Move to next iteration
      next

    }

    # Try to get good weights
    output <- try(
      sbm_weights(
        block_matrix, membership, membership_matrix, sample_size,
        snr, total_nodes, negative_proportion,
        diffusion, diffusion_range, target_condition,
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
    n = sample_size, p = total_nodes, R = output$R, skew = skew, skew_range = skew_range
  )

  # Populate results
  results <- list(
    data = data_output$data,
    parameters = list(
      nodes = nodes,
      blocks = blocks,
      sample_size = sample_size,
      skew = data_output$skew,
      density_matrix = density_matrix,
      negative_proportion = negative_proportion,
      weibull = output$params,
      diffusion = output$diffusion,
      omega = smallworldness(network = output$network, method = "omega"),
      Q = igraph::modularity(convert2igraph(abs(output$network)), output$membership)
    ),
    population = list(
      R = output$R,
      Omega = output$network,
      membership = output$membership
    ),
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
# nodes = 6; blocks = 3; sample_size = 1000
# density_matrix = matrix(0.20, nrow = blocks, ncol = blocks)
# diag(density_matrix) <- 0.90
# snr = 1; negative_proportion = 0.35
# skew = 0; skew_range = NULL
# diffusion = 0.05; diffusion_range = NULL
# target_condition = 10;
# max_correlation = 0.80; max_iterations = 100

#' @noRd
# Errors ----
# Updated 14.03.2026
simulate_sbm_errors <- function(
    nodes, blocks, density_matrix,
    snr, negative_proportion, sample_size, skew, diffusion,
    target_condition, max_correlation, max_iterations
)
{

  # Errors for 'blocks'
  typeof_error(blocks, "numeric", "simulate_sbm")
  length_error(blocks, 1, "simulate_sbm")
  range_error(blocks, c(1, Inf), "simulate_sbm")

  # Errors for 'nodes'
  typeof_error(nodes, "numeric", "simulate_sbm")
  length_error(nodes, c(1, blocks), "simulate_sbm")
  range_error(nodes, c(3, Inf), "simulate_sbm")

  # Ensure 'nodes' is the length of blocks
  if(length(nodes) != blocks){
    nodes <- rep(nodes, blocks)
  }

  # Object error for 'density_matrix'
  object_error(density_matrix, c("matrix", "data.frame", "tibble"), "simulate_sbm")

  # Check for matrix
  if(get_object_type(density_matrix) != "matrix"){
    density_matrix <- as.matrix(as.data.frame(density_matrix))
  }

  # Errors for 'density_matrix'
  typeof_error(density_matrix, "numeric", "simulate_sbm")
  length_error(density_matrix, blocks * blocks, "simulate_sbm")
  range_error(density_matrix, c(0, 1), "simulate_sbm")

  # Symmetric error for 'density_matrix'
  if(!is_symmetric(density_matrix)){
    stop("`simulate_sbm` can only handle symmetric 'density_matrix' input.")
  }

  # Error for 'snr'
  typeof_error(snr, "numeric", "simulate_sbm")
  length_error(snr, 1, "simulate_sbm")
  range_error(snr, c(0.01, Inf), "simulate_sbm")

  # Send warning if beyond bounds
  if((snr < 0.64) | (snr > 1.72)){
    warning(paste0(
      "The signal-to-ratio specified (", round(snr, 2), ") is beyond ",
      "the bounds of cataloged empirical data (0.64-1.72). "
    ))
  }

  # Errors for 'negative_proportion'
  typeof_error(negative_proportion, "numeric", "simulate_sbm")
  length_error(negative_proportion, 1, "simulate_sbm")
  range_error(negative_proportion, c(0, 1), "simulate_sbm")

  # Errors for 'sample_size'
  typeof_error(sample_size, "numeric", "simulate_sbm")
  length_error(sample_size, 1, "simulate_sbm")
  range_error(sample_size, c(1, Inf), "simulate_sbm")

  # Errors for 'skew'
  typeof_error(skew, "numeric", "simulate_sbm")
  length_error(skew, c(1, sum(nodes)), "simulate_sbm")
  range_error(skew, c(-2, 2), "simulate_sbm")

  # Errors for 'diffusion'
  typeof_error(diffusion, "numeric", "simulate_sbm")
  length_error(diffusion, 1, "simulate_sbm")
  range_error(diffusion, c(0, 1), "simulate_sbm")

  # Errors for 'target_condition'
  typeof_error(target_condition, "numeric", "simulate_sbm")
  length_error(target_condition, 1, "simulate_sbm")
  range_error(target_condition, c(1, 100), "simulate_sbm")

  # Errors for 'max_correlation'
  typeof_error(max_correlation, "numeric", "simulate_sbm")
  length_error(max_correlation, 1, "simulate_sbm")
  range_error(max_correlation, c(0, 1), "simulate_sbm")

  # Errors for 'max_iterations'
  typeof_error(max_iterations, "numeric", "simulate_sbm")
  length_error(max_iterations, 1, "simulate_sbm")
  range_error(max_iterations, c(1, Inf), "simulate_sbm")

  # Return updated nodes
  return(nodes)

}

#' @noRd
# Generate SBM ----
# Updated 14.03.2026
generate_sbm <- function(
    block_matrix, blocks, membership,
    membership_matrix, density_matrix, lower_triangle
)
{

  # Set upper triangles to FALSE
  block_matrix[!lower_triangle] <- FALSE

  # Set up lower block matrix
  for(i in seq_len(blocks)){

    # Get community index
    community_index <- membership == i

    # Get block edges
    community_block <- block_matrix[community_index, community_index]

    # Compute probabilities
    community_block[community_block] <- runif_xoshiro(sum(community_block)) < density_matrix[i,i]

    # Insert block back into block matrix
    block_matrix[community_index, community_index] <- community_block

    # Check for off-blocks
    if(i != blocks){

      # Set off-block
      for(j in (i + 1):blocks){

        # Get off-block index
        off_index <- membership == j

        # Get block edges
        off_block <- block_matrix[off_index, community_index]

        # Set index
        not_index <- !off_block

        # Compute probabilities
        off_block[not_index] <- runif_xoshiro(sum(not_index)) < density_matrix[i,j]

        # Insert block back into block matrix
        block_matrix[off_index, community_index] <- off_block

      }

    }

  }

  # Return block matrix
  return(block_matrix + t(block_matrix))

}

#' @noRd
# SBM weight generation ----
# Updated 14.03.2026
sbm_weights <- function(
    block_matrix, membership, membership_matrix, sample_size,
    snr, total_nodes, negative_proportion,
    diffusion, diffusion_range, target_condition,
    max_correlation, lower_triangle
)
{

  # Check for diffusion range
  if(!is.null(diffusion_range)){
    typeof_error(diffusion_range, "numeric") # object type error
    length_error(diffusion_range, 2) # object length error
    range_error(diffusion_range, c(0, 1)) # object range error
    diffusion <- diffusion_range
  }else{
    diffusion <- c(diffusion, diffusion)
  }

  # Reverse diffusion
  diffusion <- range(1 - diffusion)

  # Set lower triangles
  block_lower <- block_matrix[lower_triangle]
  membership_lower <- membership_matrix[lower_triangle]

  # Generate edges (returned edges are sorted)
  total_edges <- sum(block_lower)
  edges <- generate_edges(nonzero = total_edges, n = sample_size, p = total_nodes, snr = snr)

  # Obtain edges
  within_index <- membership_lower & block_lower
  total_within <- sum(within_index)

  # Set diffusion parameter
  within_reserved <- round(total_within * runif_xoshiro(1, min = diffusion[1], max = diffusion[2]))
  within_random <- total_within - within_reserved
  random_sequence <- seq_len(within_random)

  # Top edges guaranteed to go into communities
  reserved_index <- seq_len(within_reserved)

  # Identify between-community edges
  between_sequence <- total_edges - within_reserved
  between_index <- swiftelse(
    within_random > 0,
    seq_len(between_sequence)[-random_sequence],
    seq_len(between_sequence)
  )
  n_between <- length(between_index)

  # Remaining edge weights
  remaining_weights <- shuffle(edges[-reserved_index])

  # Check that there is between weights to assign
  if(n_between > 0){

    # Get signs
    signs <- swiftelse(runif_xoshiro(n_between) < negative_proportion, -1, 1)

    # Set between-community edges
    block_matrix[lower_triangle][!membership_lower & block_lower] <-
    shuffle(remaining_weights[between_index]) * signs

  }

  # Set within-community edges
  block_matrix[lower_triangle][within_index] <- shuffle(
    c(edges[reserved_index], remaining_weights[random_sequence])
  )

  # Set upper triangle to zero
  block_matrix[!lower_triangle] <- 0

  # Make symmetric
  network <- block_matrix + t(block_matrix)

  # Add names to network
  names(membership) <- row.names(network) <- colnames(network) <- paste0(
    "V", format_integer(seq_len(total_nodes), digits(total_nodes) - 1)
  )

  # Get population values
  R <- silent_call(pcor2cor(network))

  # Set condition and lambda
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
      membership = membership,
      params = attr(edges, "params"),
      diffusion = diffusion,
      lambda = lambda,
      condition = condition
    )
  )

}

