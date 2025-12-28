//! Coverage location problems: LSCP and MCLP
//!
//! LSCP: Location Set Covering Problem - minimize facilities to cover all demand
//! MCLP: Maximum Coverage Location Problem - maximize coverage with p facilities

use extendr_api::prelude::*;
use crate::locate::solver::solve_mip;

/// Solve LSCP (Location Set Covering Problem)
///
/// Minimize the number of facilities such that all demand points are covered
/// within the service radius.
pub fn solve_lscp(cost_matrix: RMatrix<f64>, service_radius: f64) -> List {
    let n_demand = cost_matrix.nrows();
    let n_facilities = cost_matrix.ncols();

    // Build coverage matrix: 1 if facility j covers demand i
    let mut coverage: Vec<Vec<f64>> = Vec::with_capacity(n_demand);
    for i in 0..n_demand {
        let row: Vec<f64> = (0..n_facilities)
            .map(|j| if cost_matrix[[i, j]] <= service_radius { 1.0 } else { 0.0 })
            .collect();
        coverage.push(row);
    }

    // Objective: minimize sum of facility selections
    let obj_coeffs: Vec<f64> = vec![1.0; n_facilities];

    // Constraints: each demand must be covered by at least one facility
    // sum_j (coverage[i][j] * y[j]) >= 1 for all i
    let constraint_rhs: Vec<f64> = vec![1.0; n_demand];
    let constraint_sense: Vec<char> = vec!['>'; n_demand];

    // All variables are binary
    let var_types: Vec<char> = vec!['B'; n_facilities];
    let lb: Vec<f64> = vec![0.0; n_facilities];
    let ub: Vec<f64> = vec![1.0; n_facilities];

    // Solve
    let result = solve_mip(
        &obj_coeffs,
        &coverage,
        &constraint_rhs,
        &constraint_sense,
        &var_types,
        &lb,
        &ub,
        false, // minimize
    );

    // Extract selected facilities (1-based for R)
    let selected: Vec<i32> = result
        .solution
        .iter()
        .enumerate()
        .filter(|(_, &v)| v > 0.5)
        .map(|(i, _)| (i + 1) as i32)
        .collect();

    // Compute coverage statistics
    let n_selected = selected.len() as i32;
    let total_demand = n_demand as i32;
    let covered_demand = compute_coverage(&coverage, &result.solution) as i32;
    let coverage_pct = (covered_demand as f64 / total_demand as f64) * 100.0;
    let status_str = format!("{:?}", result.status);

    // Check if problem is infeasible (some demand can't be covered)
    let mut uncoverable = 0i32;
    for row in &coverage {
        if row.iter().all(|&x| x < 0.5) {
            uncoverable += 1;
        }
    }

    list!(
        selected = selected,
        n_selected = n_selected,
        objective = result.objective,
        covered_demand = covered_demand,
        total_demand = total_demand,
        coverage_pct = coverage_pct,
        status = status_str,
        uncoverable_demand = uncoverable
    )
}

/// Solve MCLP (Maximum Coverage Location Problem)
///
/// Maximize total weighted demand covered with exactly p facilities.
pub fn solve_mclp(
    cost_matrix: RMatrix<f64>,
    weights: &[f64],
    service_radius: f64,
    n_facilities_to_locate: usize,
) -> List {
    let n_demand = cost_matrix.nrows();
    let n_facilities = cost_matrix.ncols();
    let p = n_facilities_to_locate;

    // Build coverage matrix
    let mut coverage: Vec<Vec<f64>> = Vec::with_capacity(n_demand);
    for i in 0..n_demand {
        let row: Vec<f64> = (0..n_facilities)
            .map(|j| if cost_matrix[[i, j]] <= service_radius { 1.0 } else { 0.0 })
            .collect();
        coverage.push(row);
    }

    // Variables:
    // y[j] = 1 if facility j is selected (j = 0..n_facilities-1)
    // z[i] = 1 if demand i is covered (i = n_facilities..n_facilities+n_demand-1)
    let n_vars = n_facilities + n_demand;

    // Objective: maximize sum of weighted demand covered
    // max sum_i (weights[i] * z[i])
    let mut obj_coeffs = vec![0.0; n_vars];
    for i in 0..n_demand {
        obj_coeffs[n_facilities + i] = weights[i];
    }

    // Constraints:
    // 1. sum_j y[j] = p (exactly p facilities)
    // 2. z[i] <= sum_j (coverage[i][j] * y[j]) for all i (can only be covered if facility exists)

    let n_constraints = 1 + n_demand;
    let mut constraint_matrix: Vec<Vec<f64>> = Vec::with_capacity(n_constraints);
    let mut constraint_rhs: Vec<f64> = Vec::with_capacity(n_constraints);
    let mut constraint_sense: Vec<char> = Vec::with_capacity(n_constraints);

    // Constraint 1: sum_j y[j] = p
    let mut row1 = vec![0.0; n_vars];
    for j in 0..n_facilities {
        row1[j] = 1.0;
    }
    constraint_matrix.push(row1);
    constraint_rhs.push(p as f64);
    constraint_sense.push('=');

    // Constraint 2: z[i] - sum_j (coverage[i][j] * y[j]) <= 0
    for i in 0..n_demand {
        let mut row = vec![0.0; n_vars];
        row[n_facilities + i] = 1.0; // z[i]
        for j in 0..n_facilities {
            row[j] = -coverage[i][j]; // -coverage[i][j] * y[j]
        }
        constraint_matrix.push(row);
        constraint_rhs.push(0.0);
        constraint_sense.push('<');
    }

    // Variable types and bounds
    let var_types = vec!['B'; n_vars]; // All binary
    let lb = vec![0.0; n_vars];
    let ub = vec![1.0; n_vars];

    // Solve
    let result = solve_mip(
        &obj_coeffs,
        &constraint_matrix,
        &constraint_rhs,
        &constraint_sense,
        &var_types,
        &lb,
        &ub,
        true, // maximize
    );

    // Extract selected facilities (1-based for R)
    let selected: Vec<i32> = result
        .solution
        .iter()
        .take(n_facilities)
        .enumerate()
        .filter(|(_, &v)| v > 0.5)
        .map(|(i, _)| (i + 1) as i32)
        .collect();

    // Coverage statistics
    let covered_weight: f64 = result
        .solution
        .iter()
        .skip(n_facilities)
        .zip(weights.iter())
        .filter(|(&z, _)| z > 0.5)
        .map(|(_, &w)| w)
        .sum();

    let total_weight: f64 = weights.iter().sum();
    let n_selected = selected.len() as i32;
    let coverage_pct = (covered_weight / total_weight) * 100.0;

    list!(
        selected = selected,
        n_selected = n_selected,
        objective = result.objective,
        covered_weight = covered_weight,
        total_weight = total_weight,
        coverage_pct = coverage_pct
    )
}

/// Count number of demand points covered by selected facilities
fn compute_coverage(coverage: &[Vec<f64>], solution: &[f64]) -> usize {
    coverage
        .iter()
        .filter(|row| {
            row.iter()
                .zip(solution.iter())
                .any(|(&cov, &sel)| cov > 0.5 && sel > 0.5)
        })
        .count()
}
