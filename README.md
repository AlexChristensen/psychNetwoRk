# psychNetwoRk
#### Data Simulation Based on Psychometric Networks

---

## Overview

This package simulates multivariate normal data from **Gaussian graphical models** (GGMs) with realistic psychometric network structure. Two topology generators are included: `simulate_smallworld()` produces small-world networks via degree-weighted Watts-Strogatz rewiring, and `simulate_sbm()` produces community-structured networks via a stochastic block model. Both are calibrated to edge weight distributions empirically observed across 194 psychometric networks (Huth et al., 2025) through `weibull_parameters()`, a Seemingly Unrelated Regression model that predicts Weibull shape and scale from network size, sample size, and signal-to-noise ratio.

Additional tools include `smallworldness()` (small-world network metrics), `proxswap_lattice()` (proximity-swap ring lattice construction), `categorize()` (categorizing continuous data with skew), and `add_population_error()` (adds population/model error so that the true network fits at a target CFI, TLI, RMSEA, or SRMR using one of four methods).

To estimate networks from the simulated data with $L_0$ penalties, see [{L0ggm}](https://github.com/AlexChristensen/L0ggm).

---

## Installation

```r
# Ensure {remotes} is installed
if(!"remotes" %in% row.names(installed.packages())){
  install.packages("remotes")
}

# Install from GitHub
remotes::install_github("AlexChristensen/psychNetwoRk")
```

---

## Usage

```r
# Load package
library(psychNetwoRk)

# Simulate smallworld data
basic_smallworld <- simulate_smallworld(
  nodes = 20, # 20 nodes in the network
  density = 0.30, # moderate initial lattice connectivity
  rewire = 0.20, # 20% rewiring probability
  sample_size = 1000 # number of cases = 1000
)

# Simulate community-structured data
dm <- matrix(0.20, nrow = 3, ncol = 3)
diag(dm) <- 0.90
sbm <- simulate_sbm(nodes = 6, blocks = 3, density_matrix = dm, sample_size = 500)

# Add population (model) error so the true network has CFI = 0.95
sbm_error <- add_population_error(sbm, method = "cb", metric = "CFI", fit_value = 0.95)
```

---

## Simulation

{psychNetwoRk} provides two functions for generating GGM data with realistic psychometric network structure. Both are calibrated to edge weight distributions observed across 194 empirical psychometric networks (Huth et al., 2025) via `weibull_parameters()`, with shared support for skewed data generation, positive-definiteness conditioning, and reproducible rejection diagnostics.

<details>
<summary><strong><code>simulate_smallworld()</code> — Small-world network structure</strong></summary>

Generates GGM data from a small-world topology via a three-stage process closely aligned with empirical psychometric network properties.

**Stage 1 — Ring lattice.** A ring lattice is constructed with `neighbors + 1` nearest-neighbor connections per node, where `neighbors` is derived from `nodes` and `density`. The extra neighbor ensures the lattice always has more edges than the target, guaranteeing pruning can proceed.

**Stage 2 — Pruning.** The lattice is randomly pruned to the target `density` by removing edges uniformly at random while maintaining graph connectivity. This creates degree heterogeneity before any rewiring, providing a non-uniform prior for the next stage.

**Stage 3 — Degree-weighted rewiring.** Each edge is independently rewired with probability `rewire`. Unlike standard Watts-Strogatz, the new endpoint is selected proportionally to $\sqrt{d_i + d_k}$ (combined degree of the candidate pair), moderating the extreme hubs of linear preferential attachment while still producing realistic degree heterogeneity.

**Edge weights.** Retained lattice edges receive weights by neighbor-distance priority (nearer = stronger); rewired edges receive lowest priority. Absolute weights are drawn from a Weibull distribution via `weibull_parameters()`.

| Argument | Description |
|---|---|
| `nodes` | Number of nodes (8–54 recommended) |
| `density` | Target edge density after pruning |
| `rewire` | Rewiring probability (small-world regime: ~0.01–0.30) |
| `snr` | Signal-to-noise ratio $\bar{\lvert w \rvert} / \mathrm{SD}(\lvert w \rvert)$; default `1` |
| `negative_proportion` | Fraction of edges signed negative; default drawn from empirical distribution |
| `sample_size` | Number of observations to generate |
| `skew` / `skew_range` | Per-variable skew in $[-2, 2]$; default `0` |

**Reference:** Christensen et al. (2026; [arXiv:2609.05678](https://doi.org/10.48550/arXiv.2609.05678)).

**Returns:** a named list with `data`, `parameters` (including `omega` smallworldness), `population` (`R`, `Omega`), and `convergence`.

```r
result <- simulate_smallworld(
  nodes = 20, density = 0.30, rewire = 0.20, sample_size = 500
)
```

</details>

<details>
<summary><strong><code>simulate_sbm()</code> — Stochastic block model structure</strong></summary>

Generates GGM data with explicit community (block) structure, where edge density is controlled separately within and between communities via a `blocks × blocks` density matrix.

**Structure.** Nodes are partitioned into `blocks` communities of size `nodes` (scalar or per-block vector). For each pair of communities $(i, j)$, edges are included independently with probability `density_matrix[i, j]`. The resulting graph is required to be connected; disconnected draws are rejected and resampled.

**Edge weights and diffusion.** Absolute weights are drawn from a Weibull distribution via `weibull_parameters()`. The `diffusion` parameter controls the minimum fraction of high-weight edges reserved for within-community positions: `1 - diffusion` of within-community slots are filled from the top-ranked weight draws, with the remaining slots drawn randomly. Because random draws can also land within communities by chance, the actual within-community weight advantage typically exceeds `1 - diffusion`. Newman-Girvan modularity $Q$ (returned in `parameters`) summarizes the community contrast actually achieved.

**Between-community sign.** Each between-community edge is independently signed negative with probability `negative_proportion`, modeling inhibitory connections. Within-community edges are always positive.

| Argument | Description |
|---|---|
| `nodes` | Nodes per block (scalar or length-`blocks` vector; min 3 per block) |
| `blocks` | Number of community blocks |
| `density_matrix` | `blocks × blocks` symmetric matrix of edge probabilities |
| `snr` | Signal-to-noise ratio; default `1` |
| `negative_proportion` | Fraction of between-community edges signed negative; default from empirical distribution |
| `diffusion` | Minimum proportion of top edges placed outside communities; default `0.30` |
| `diffusion_range` | If provided, `diffusion` is drawn uniformly from this interval each call |
| `sample_size` | Number of observations to generate |
| `skew` / `skew_range` | Per-variable skew in $[-2, 2]$; default `0` |

**Reference:** Christensen & Choi (2026; [doi:10.1037/met0000853](https://doi.org/10.1037/met0000853)).

**Returns:** a named list with `data`, `parameters` (including `Q` modularity and `weibull` edge weight parameters), `population` (`R`, `Omega`, `membership`), and `convergence`.

```r
dm <- matrix(0.20, nrow = 3, ncol = 3)
diag(dm) <- 0.90

result <- simulate_sbm(
  nodes = 6, blocks = 3, density_matrix = dm, sample_size = 500
)
```

</details>

<details>
<summary><strong><code>weibull_parameters()</code> — Empirically calibrated edge weight distribution</strong></summary>

Predicts the Weibull shape $k$ and scale $\lambda_W$ that characterize absolute partial correlation edge weights for a network of given size and sample size. Predictions come from a Seemingly Unrelated Regression (SUR) model fitted to 194 empirical psychometric networks (Huth et al., 2025), where absolute partial correlations were found to follow a Weibull distribution more consistently than Beta, Gamma, or log-normal alternatives.

**Predictors.** The two SUR equations have an asymmetric structure:

- **Shape** is predicted from `snr` and `rlp` $= 1/\log(p)$, capturing the diminishing marginal effect of network size on weight concentration. Shape governs the spread of the edge weight distribution and is not affected by sampling precision.
- **Scale** is predicted from `snr`, `rlp`, and `scaling` $= \sqrt{1/(n - p - 2)}$ (the standard error of partial correlations). Scale governs the typical magnitude of edge weights and is directly affected by estimation precision.

**Model fit.** Shape: $R^2 = 0.887$, RMSE $= 0.048$. Scale: $R^2 = 0.885$, RMSE $= 0.011$. Residual correlation between equations: $0.261$ (motivating SUR). Empirical ranges: shape $\in [0.72, 1.63]$ (M $= 1.07$, SD $= 0.14$); scale $\in [0.03, 0.19]$ (M $= 0.10$, SD $= 0.03$).

| Argument | Description |
|---|---|
| `nodes` | Number of nodes (8–54 to stay within empirical training range) |
| `sample_size` | Sample size of the dataset |
| `snr` | Signal-to-noise ratio $\bar{\lvert w \rvert} / \mathrm{SD}(\lvert w \rvert)$; default `1` |
| `bootstrap` | If `TRUE`, adds a sampled SUR residual — use in simulation for replication-to-replication variability |

```r
# Predicted parameters for a 15-node network, n = 300
weibull_parameters(nodes = 15, sample_size = 300)

# With bootstrapped residuals for Monte Carlo simulation
weibull_parameters(nodes = 15, sample_size = 300, bootstrap = TRUE)
```

</details>

<details>
<summary><strong><code>add_population_error()</code> — Population (model) error</strong></summary>

Perturbs the population correlation matrix of a `simulate_sbm()` or `simulate_smallworld()` object so that the true network no longer fits perfectly. The error is optimized until the true network, refit to the new matrix by constrained maximum likelihood, reaches a target value of a fit index. Data are then re-generated from the perturbed matrix (same sample size and skew).

| Method | Description |
|---|---|
| `"cb"` | Cudeck & Browne (1992): error on the cells the true network constrains to zero (orthogonal to the network); default |
| `"tkl"` | Tucker, Koopman, & Linn (1969): random minor factors with geometrically decaying loadings; the proportion of variance from minor factors is optimized |
| `"wb"` | Wu & Browne (2015): inverted-Wishart draw centered on the original matrix; the dispersion (approximately RMSEA$^2$) is optimized |
| `"hb"` | Hu & Bentler (1998): underparameterized misspecification; a random subset of the network's zero edges is non-zero in the population and their size is optimized |

| Argument | Description |
|---|---|
| `ggm_object` | Output from `simulate_sbm()` or `simulate_smallworld()` |
| `method` | Error method (see above); default `"cb"` |
| `metric` | Fit index to target: `"CFI"` (default), `"TLI"`, `"RMSEA"`, or `"SRMR"` |
| `fit_value` | Target value of `metric` in $(0, 1)$; default is the conventional cutoff (CFI/TLI `0.95`, RMSEA `0.06`, SRMR `0.08`) |
| `max_residual` | Largest absolute error added to any cell; default `0.10`. If the target cannot be reached within this limit, the closest solution is returned with a warning |
| `tolerance` | Maximum absolute difference between `fit_value` and the achieved `metric`; default `0.001` |
| `max_iterations` | Maximum optimization attempts (random starts or draws); default `10` |
| `control` | Method options: `minor_factors` (`max(50, 2 * nodes)`) and `decay` (`0.02`) for `"tkl"`; `omitted` (`0.10`, proportion of zero edges) for `"hb"` |

**References:** Cudeck & Browne (1992); Tucker et al. (1969); Wu & Browne (2015); Hu & Bentler (1998).

**Returns:** the input object with class `"ggm_pe"` appended and a `population_error` element containing `method`, `metric`, `fit_value`, `achieved`, `fit` (CFI, TLI, RMSEA, and SRMR achieved), `delta`, `error`, `max_error`, `max_residual`, `R_difference`, and `iterations`. `R_difference` is near zero only when the error is orthogonal to the true network (`"cb"`).

```r
dm <- matrix(0.20, nrow = 3, ncol = 3)
diag(dm) <- 0.90

sbm <- simulate_sbm(nodes = 6, blocks = 3, density_matrix = dm, sample_size = 1000)

# True network fits at RMSEA = 0.05 using Wu & Browne (2015)
result <- add_population_error(
  sbm, method = "wb", metric = "RMSEA", fit_value = 0.05, max_residual = 0.20
)
```

</details>

<details>
<summary><strong>Simulation internals (<code>helpers_simulation</code>)</strong></summary>

The shared simulation pipeline is implemented in `helpers_simulation.R` and handles tasks common to both `simulate_smallworld()` and `simulate_sbm()`.

**`generate_edges()`** draws absolute edge weights from a Weibull distribution parameterized by `weibull_parameters()`. Draws are repeated until no weight exceeds the empirical maximum ($\approx 0.845$), and are then sampled proportionally to $w^2$ (preferring larger weights for occupied edges). The Weibull parameters are restricted to empirical bounds (shape $\in [0.70, 1.70]$, scale $\in [0.03, 0.19]$); draws outside these bounds are discarded.

**`condition_network()`** recovers positive definiteness by adding a minimal ridge penalty $\lambda$ to the diagonal of the precision matrix via `uniroot()`, targeting a user-specified condition number. Maximum shrinkage is capped at approximately 23% following Peeters et al. (2020).

**`simulate_data()`** generates multivariate normal observations via a Cholesky decomposition of the population correlation matrix $\mathbf{R}$, then applies variable-level skew transformations. Skew values are rounded to the nearest 0.05 increment within $[-2, 2]$; if `skew_range` is provided, per-variable skew is drawn uniformly from the interval.

</details>

---

## References

Christensen, A. P., & Choi, J. (2026). Network accuracy across levels of analysis using stochastic block models. Psychological Methods. https://doi.org/10.1037/met0000853

Christensen, A. P., Choi, J., Yang, H., & Xue, L. (2026). Adaptive regularization via extreme value distributions for Gaussian graphical models (arXiv:2609.05678). arXiv. https://doi.org/10.48550/arXiv.2609.05678

Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value. *Psychometrika*, *57*(3), 357&ndash;369. https://doi.org/10.1007/BF02295424

Hu, L., & Bentler, P. M. (1998). Fit indices in covariance structure modeling: Sensitivity to underparameterized model misspecification. *Psychological Methods*, *3*(4), 424&ndash;453. https://doi.org/10.1037/1082-989X.3.4.424

Huth, K. B. S., Haslbeck, J. M. B., Keetelaar, S., Van Holst, R. J., & Marsman, M. (2025). Statistical evidence in psychological networks. *Nature Human Behaviour*. https://doi.org/10.1038/s41562-025-02314-2

Peeters, C. F., van de Wiel, M. A., & van Wieringen, W. N. (2020). The spectral condition number plot for regularization parameter evaluation. *Computational Statistics*, *35*(2), 629&ndash;646. https://doi.org/10.1007/s00180-019-00912-z

Tucker, L. R., Koopman, R. F., & Linn, R. L. (1969). Evaluation of factor analytic research procedures by means of simulated correlation matrices. *Psychometrika*, *34*(4), 421&ndash;459. https://doi.org/10.1007/BF02290601

Wu, H., & Browne, M. W. (2015). Quantifying adventitious error in a covariance structure as a random effect. *Psychometrika*, *80*(3), 571&ndash;600. https://doi.org/10.1007/s11336-015-9451-3
