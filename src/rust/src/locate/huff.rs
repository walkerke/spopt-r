//! Huff Model for Spatial Interaction and Market Share Analysis
//!
//! Computes probability surfaces based on distance decay and attractiveness factors.
//! Used for retail site selection and market share analysis.

use extendr_api::prelude::*;

/// Compute Huff model probabilities
///
/// Formula: P_ij = (A_j^α × D_ij^β) / Σ_k(A_k^α × D_ik^β)
///
/// # Arguments
/// * `cost_matrix` - Distance/cost matrix (n_demand x n_stores)
/// * `attractiveness` - Attractiveness values for each store (already raised to exponents if multiple)
/// * `distance_exponent` - Distance decay exponent (typically negative, e.g., -1.5)
///
/// # Returns
/// List with probability matrix, market shares, and expected sales
pub fn compute_probabilities(
    cost_matrix: RMatrix<f64>,
    attractiveness: &[f64],
    distance_exponent: f64,
    sales_potential: Option<&[f64]>,
) -> List {
    let n_demand = cost_matrix.nrows();
    let n_stores = cost_matrix.ncols();

    if attractiveness.len() != n_stores {
        return list!(
            error = format!(
                "attractiveness length ({}) must equal number of stores ({})",
                attractiveness.len(), n_stores
            )
        );
    }

    // Compute utility matrix: U_ij = A_j × D_ij^β
    // Note: attractiveness is already raised to its exponent(s) by R wrapper
    let mut utility_matrix: Vec<f64> = Vec::with_capacity(n_demand * n_stores);

    for i in 0..n_demand {
        for j in 0..n_stores {
            let dist = cost_matrix[[i, j]];
            // Handle zero distance (same location)
            let dist_component = if dist.abs() < 1e-10 {
                // Very close - use small distance to avoid infinity
                (0.001_f64).powf(distance_exponent)
            } else {
                dist.powf(distance_exponent)
            };
            let utility = attractiveness[j] * dist_component;
            utility_matrix.push(utility);
        }
    }

    // Compute probabilities: P_ij = U_ij / Σ_k U_ik
    let mut prob_matrix: Vec<f64> = Vec::with_capacity(n_demand * n_stores);
    let mut row_sums: Vec<f64> = Vec::with_capacity(n_demand);

    for i in 0..n_demand {
        let row_start = i * n_stores;
        let row_sum: f64 = utility_matrix[row_start..(row_start + n_stores)].iter().sum();
        row_sums.push(row_sum);

        for j in 0..n_stores {
            let prob = if row_sum > 1e-10 {
                utility_matrix[row_start + j] / row_sum
            } else {
                // All stores equally inaccessible
                1.0 / n_stores as f64
            };
            prob_matrix.push(prob);
        }
    }

    // Compute market share for each store (average probability weighted by sales potential)
    let mut market_shares: Vec<f64> = vec![0.0; n_stores];
    let mut expected_sales: Vec<f64> = vec![0.0; n_stores];

    let total_potential: f64;

    match sales_potential {
        Some(potential) => {
            total_potential = potential.iter().sum();
            for j in 0..n_stores {
                let mut weighted_prob_sum = 0.0;
                let mut sales_sum = 0.0;
                for i in 0..n_demand {
                    let prob = prob_matrix[i * n_stores + j];
                    weighted_prob_sum += prob * potential[i];
                    sales_sum += prob * potential[i];
                }
                market_shares[j] = if total_potential > 0.0 {
                    weighted_prob_sum / total_potential
                } else {
                    0.0
                };
                expected_sales[j] = sales_sum;
            }
        }
        None => {
            // Unweighted - each demand point counts equally
            total_potential = n_demand as f64;
            for j in 0..n_stores {
                let mut prob_sum = 0.0;
                for i in 0..n_demand {
                    prob_sum += prob_matrix[i * n_stores + j];
                }
                market_shares[j] = prob_sum / n_demand as f64;
                expected_sales[j] = prob_sum; // Just count of "expected customers"
            }
        }
    }

    // Find primary store for each demand point (highest probability)
    let primary_store: Vec<i32> = (0..n_demand)
        .map(|i| {
            let row_start = i * n_stores;
            let mut best_j = 0;
            let mut best_prob = prob_matrix[row_start];
            for j in 1..n_stores {
                if prob_matrix[row_start + j] > best_prob {
                    best_prob = prob_matrix[row_start + j];
                    best_j = j;
                }
            }
            (best_j + 1) as i32 // 1-based for R
        })
        .collect();

    // Compute entropy (measure of competition/uncertainty)
    // H_i = -Σ_j P_ij × log(P_ij)
    let entropy: Vec<f64> = (0..n_demand)
        .map(|i| {
            let row_start = i * n_stores;
            let mut h = 0.0;
            for j in 0..n_stores {
                let p = prob_matrix[row_start + j];
                if p > 1e-10 {
                    h -= p * p.ln();
                }
            }
            h
        })
        .collect();

    list!(
        probabilities = prob_matrix,
        market_shares = market_shares,
        expected_sales = expected_sales,
        primary_store = primary_store,
        entropy = entropy,
        total_potential = total_potential
    )
}
