/* Written on 24 March 2026 by Alexander P. Christensen (alexpaulchristensen@gmail.com)
 * with assistance from Claude Sonnet 4.6 for constructing logic implementation in C
 */


/*
 * proxswap_lattice.c
 * ------------------
 * C implementations of proximity_pass and swapping_pass for the
 * proximity-swap ring lattice construction algorithm.
 * 
 * License: CC0-1.0
 * To the extent possible under law, the author(s) have dedicated all
 * copyright and related and neighboring rights to this software to the
 * public domain worldwide.  This software is distributed without any
 * warranty.  See <https://creativecommons.org/publicdomain/zero/1.0/>.
 *
 * Matrices are stored and passed as flat integer vectors in column-major
 * order (R's native layout).  Element [i, j] (0-indexed) is accessed as
 * flat[i + j * n] throughout.
 *
 * Pair arrays
 * -----------
 * build_pairs() in R returns a list of matrices, one per ring-distance band.
 * Before calling proximity_pass_c the caller must flatten these into three
 * parallel vectors supplied as integer SEXPs:
 *
 *   pair_rows_r   – concatenated 1-indexed row indices, bands in distance order
 *   pair_cols_r   – corresponding 1-indexed col indices (same length)
 *   pair_counts_r – number of pairs in each band (length == number of bands)
 *
 * Ring matrix type
 * ----------------
 * Both functions expect the ring adjacency matrix passed as an integer vector
 * (0 = absent, 1 = present).  Convert a logical matrix in R with
 * storage.mode(ring) <- "integer" before calling.
 *
 * Tie-breaking note (proximity_pass_c)
 * -------------------------------------
 * Within each distance band, eligible pairs are sorted by descending
 * min(budget[row], budget[col]).  A stable insertion sort preserves the
 * original within-band order on equal keys, matching R's order() behavior.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* =========================================================================
 * proximity_pass_c
 *
 * Greedily assigns edges in strictly increasing ring-distance order.
 * Within each distance band the eligible pairs (both endpoints have
 * remaining degree budget and are not yet connected) are sorted by
 * descending min(budget_i, budget_j) so that high-need pairs receive short
 * connections first.  Pairs are then assigned one at a time with per-pair
 * budget re-checks.  Returns early once all budgets reach zero.
 *
 * Returns a named list: ring (updated adjacency), budget (updated per-node
 * degree deficits), total_budget (sum of remaining deficits).
 * ========================================================================= */
SEXP proximity_pass_c(
    SEXP nodes_r, SEXP ring_r, SEXP budget_r,
    SEXP pair_rows_r, SEXP pair_cols_r, SEXP pair_counts_r
){

    // Number of nodes in the ring
    int n = asInteger(nodes_r);

    // Number of distance bands (== max ring distance == floor(n/2))
    int n_dist = length(pair_counts_r);

    // Allocate output adjacency matrix buffer
    SEXP ring_out = PROTECT(allocVector(INTSXP, n * n));

    // Allocate output degree-deficit buffer
    SEXP budget_out = PROTECT(allocVector(INTSXP, n));

    // Get pointer to output adjacency matrix
    int *ring = INTEGER(ring_out);

    // Get pointer to output degree-deficit vector
    int *budget = INTEGER(budget_out);

    // Copy input adjacency matrix into output buffer
    int *ri = INTEGER(ring_r);
    for (int k = 0; k < n * n; k++) ring[k] = ri[k];

    // Copy input degree-deficit vector into output buffer
    int *bi = INTEGER(budget_r);
    for (int k = 0; k < n; k++) budget[k] = bi[k];

    // Get read-only pointers to flat pair arrays
    int *pair_rows   = INTEGER(pair_rows_r);
    int *pair_cols   = INTEGER(pair_cols_r);
    int *pair_counts = INTEGER(pair_counts_r);

    // Compute total remaining degree deficit across all nodes
    int total_budget = 0;
    for (int k = 0; k < n; k++) total_budget += budget[k];

    // Temporary storage: indices of eligible pairs within the current band
    // Maximum pairs per band is n (distance-1 band has exactly n upper-triangle pairs)
    int *elig  = (int *) R_alloc(n, sizeof(int));

    // Temporary storage: sorted permutation of eligible-pair indices
    int *order = (int *) R_alloc(n, sizeof(int));

    // Running offset into the flat pair arrays
    int offset = 0;

    // Process each distance band in increasing ring-distance order
    for (int di = 0; di < n_dist; di++) {

        // Number of candidate pairs in this band
        int count = pair_counts[di];

        // Save band start offset before advancing
        int base = offset;

        // Advance offset to next band
        offset += count;

        // Collect indices of eligible pairs:
        // both endpoints must have budget > 0 and not yet be connected
        int n_elig = 0;
        for (int p = 0; p < count; p++) {

            // Convert 1-indexed pair to 0-indexed node indices
            int r = pair_rows[base + p] - 1;
            int c = pair_cols[base + p] - 1;

            // Check eligibility
            if (budget[r] > 0 && budget[c] > 0 && !ring[r + c * n]) {
                elig[n_elig++] = p;
            }
        }

        // Skip band if no eligible pairs were found
        if (n_elig == 0) continue;

        // Initialize sort permutation to identity
        for (int e = 0; e < n_elig; e++) order[e] = e;

        // Stable insertion sort: descending min(budget[r], budget[c])
        // Stable sort preserves within-band order on ties, matching R's order()
        for (int e = 1; e < n_elig; e++) {

            // Element being inserted into the sorted prefix
            int tmp = order[e];
            int p_tmp = elig[tmp];
            int r_tmp = pair_rows[base + p_tmp] - 1;
            int c_tmp = pair_cols[base + p_tmp] - 1;

            // Sort key: smaller of the two endpoint budgets
            int key = (budget[r_tmp] < budget[c_tmp]) ? budget[r_tmp] : budget[c_tmp];

            // Scan backwards past elements with a strictly smaller key
            int pos = e - 1;
            while (pos >= 0) {
                int p_pos  = elig[order[pos]];
                int r_pos  = pair_rows[base + p_pos] - 1;
                int c_pos  = pair_cols[base + p_pos] - 1;
                int key_pos = (budget[r_pos] < budget[c_pos]) ? budget[r_pos] : budget[c_pos];

                // Stop once we reach an element with a key >= current key
                if (key_pos >= key) break;

                // Shift element one position right
                order[pos + 1] = order[pos];
                pos--;
            }

            // Insert element in its correct sorted position
            order[pos + 1] = tmp;
        }

        // Assign edges in sorted order, re-checking budgets before each assignment
        // (budgets can change as earlier pairs in the band consume degree)
        for (int e = 0; e < n_elig; e++) {

            // Look up the pair for this sorted position
            int p = elig[order[e]];
            int r = pair_rows[base + p] - 1;
            int c = pair_cols[base + p] - 1;

            // Re-check budgets: both may have been consumed by earlier pairs
            if (budget[r] > 0 && budget[c] > 0) {

                // Add edge (symmetric)
                ring[r + c * n] = 1;
                ring[c + r * n] = 1;

                // Decrement both endpoint budgets
                budget[r]--;
                budget[c]--;

                // Update global deficit count
                total_budget -= 2;
            }
        }

        // Early exit once all degree deficits are satisfied
        if (total_budget < 1) break;
    }

    // Build named return list: ring, budget, total_budget
    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SEXP names  = PROTECT(allocVector(STRSXP, 3));
    SEXP tb_out = PROTECT(ScalarInteger(total_budget));

    // Set list element names
    SET_STRING_ELT(names, 0, mkChar("ring"));
    SET_STRING_ELT(names, 1, mkChar("budget"));
    SET_STRING_ELT(names, 2, mkChar("total_budget"));

    // Populate list elements
    SET_VECTOR_ELT(result, 0, ring_out);
    SET_VECTOR_ELT(result, 1, budget_out);
    SET_VECTOR_ELT(result, 2, tb_out);

    // Attach names attribute
    setAttrib(result, R_NamesSymbol, names);

    // Unprotect all five PROTECT calls: ring_out, budget_out, result, names, tb_out
    UNPROTECT(5);

    // Return named list
    return result;
}

/* =========================================================================
 * swapping_pass_c
 *
 * Resolves any residual degree deficit remaining after proximity_pass_c.
 * At each iteration the highest-deficit node i is connected to its nearest
 * available ring neighbour, scanning interleaved clockwise and counter-
 * clockwise positions in distance order.  When no direct partner with
 * remaining budget exists, an edge swap is performed: a nearby unconnected
 * node j is found, one of j's existing edges to k is removed, and a new
 * edge (i, j) is added.  Node k recovers its budget for resolution in a
 * subsequent iteration.  Iteration stops when all deficits are resolved, no
 * swap can be found, or the cap (2 * n^2) is reached.
 *
 * Interleaving
 * ------------
 * Ring positions are visited by interleaving clockwise and counter-clockwise
 * offsets from i:  cw_1, ccw_1, cw_2, ccw_2, ...  Duplicate positions
 * (which arise at distance n/2 when n is even) are skipped via a visited
 * flag array.
 *
 * Tie-breaking note (swap selection)
 * ------------------------------------
 * The first eligible neighbor k of j in node-index order is used, matching
 * R's which()[1] behavior.
 *
 * Returns a named list: ring (updated adjacency), remaining (total
 * unresolved degree deficit, 0 on full success).
 * ========================================================================= */
SEXP swapping_pass_c(
    SEXP nodes_r, SEXP ring_r, SEXP budget_r, SEXP total_budget_r
){

    // Number of nodes in the ring
    int n = asInteger(nodes_r);

    // Maximum ring distance (floor(n/2))
    int n_dist = n / 2;

    // Total remaining degree deficit passed in from proximity_pass_c
    int total_budget = asInteger(total_budget_r);

    // Allocate output adjacency matrix buffer
    SEXP ring_out = PROTECT(allocVector(INTSXP, n * n));

    // Allocate output degree-deficit buffer
    SEXP budget_out = PROTECT(allocVector(INTSXP, n));

    // Get pointer to output adjacency matrix
    int *ring = INTEGER(ring_out);

    // Get pointer to output degree-deficit vector
    int *budget = INTEGER(budget_out);

    // Copy input adjacency matrix into output buffer
    int *ri = INTEGER(ring_r);
    for (int k = 0; k < n * n; k++) ring[k] = ri[k];

    // Copy input degree-deficit vector into output buffer
    int *bi = INTEGER(budget_r);
    for (int k = 0; k < n; k++) budget[k] = bi[k];

    // Interleaved ring-position buffer: at most 2 * n_dist entries before deduplication
    int *interleaved = (int *) R_alloc(2 * n_dist, sizeof(int));

    // Per-node visited flag used to deduplicate the interleaved sequence
    int *visited = (int *) R_alloc(n, sizeof(int));

    // Run at most 2 * n^2 repair iterations
    int max_iter = 2 * n * n;

    for (int iter = 0; iter < max_iter; iter++) {

        // Stop when all degree deficits are resolved
        if (total_budget == 0) break;

        // Find the node with the largest remaining degree deficit
        int i = 0;
        for (int k = 1; k < n; k++) {
            if (budget[k] > budget[i]) i = k;
        }

        // Stop if even the largest deficit is zero
        if (budget[i] == 0) break;

        // Generate interleaved clockwise / counter-clockwise positions from i,
        // sorted by ascending ring distance and deduplicated via visited flags
        for (int k = 0; k < n; k++) visited[k] = 0;

        // Mark i itself as visited so it is never a candidate
        visited[i] = 1;

        // Build interleaved sequence
        int n_pos = 0;
        for (int d = 1; d <= n_dist; d++) {

            // Clockwise neighbour at ring distance d
            int cw = (i + d) % n;
            if (!visited[cw]) {
                visited[cw] = 1;
                interleaved[n_pos++] = cw;
            }

            // Counter-clockwise neighbour at ring distance d
            int ccw = ((i - d) % n + n) % n;
            if (!visited[ccw]) {
                visited[ccw] = 1;
                interleaved[n_pos++] = ccw;
            }
        }

        // --- Attempt 1: direct connection ---
        // Connect i to the nearest ring position that has remaining budget
        // and is not yet an existing neighbour of i
        int direct = -1;
        for (int pi = 0; pi < n_pos; pi++) {
            int j = interleaved[pi];
            if (budget[j] > 0 && !ring[i + j * n]) {
                direct = j;
                break;
            }
        }

        // If a direct partner was found, connect and continue
        if (direct >= 0) {

            // Add edge (symmetric)
            ring[i + direct * n] = 1;
            ring[direct + i * n] = 1;

            // Decrement both endpoint budgets
            budget[i]--;
            budget[direct]--;

            // Update global deficit count
            total_budget -= 2;

            // Move to next iteration
            continue;
        }

        // --- Attempt 2: edge swap ---
        // No direct partner available; find a nearby unconnected node j,
        // steal one of j's existing edges to k, then connect i to j.
        // Node k recovers its budget for resolution in a later iteration.
        int swapped = 0;

        // Iterate over nearby unconnected ring positions in distance order
        for (int pi = 0; pi < n_pos && !swapped; pi++) {
            int j = interleaved[pi];

            // j must not already be connected to i
            if (ring[i + j * n]) continue;

            // Find a current neighbour k of j that is not i
            // and is not already connected to i
            for (int k = 0; k < n; k++) {

                // k must not be i itself
                if (k == i) continue;

                // k must currently be a neighbour of j
                if (!ring[j + k * n]) continue;

                // k must not already be connected to i
                if (ring[i + k * n]) continue;

                // Perform the swap: remove (j, k), add (i, j)
                ring[j + k * n] = 0;
                ring[k + j * n] = 0;
                ring[i + j * n] = 1;
                ring[j + i * n] = 1;

                // i gains an edge (budget decreases)
                budget[i]--;

                // k loses an edge (budget recovers)
                budget[k]++;

                // total_budget is unchanged

                // Record that a swap occurred and stop scanning
                swapped = 1;
                break;
            }
        }

        // If no swap could be performed the graph is stuck; exit early
        if (!swapped) break;
    }

    // Build named return list: ring, remaining
    SEXP result    = PROTECT(allocVector(VECSXP, 2));
    SEXP names     = PROTECT(allocVector(STRSXP, 2));
    SEXP remaining = PROTECT(ScalarInteger(total_budget));

    // Set list element names
    SET_STRING_ELT(names, 0, mkChar("ring"));
    SET_STRING_ELT(names, 1, mkChar("remaining"));

    // Populate list elements
    SET_VECTOR_ELT(result, 0, ring_out);
    SET_VECTOR_ELT(result, 1, remaining);

    // Attach names attribute
    setAttrib(result, R_NamesSymbol, names);

    // Unprotect all five PROTECT calls: ring_out, budget_out, result, names, remaining
    UNPROTECT(5);

    // Return named list
    return result;
}
