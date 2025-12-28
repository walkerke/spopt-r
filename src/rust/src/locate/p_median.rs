//! P-Median Problem
//!
//! Minimize total weighted distance by locating exactly p facilities.

use extendr_api::prelude::*;
use crate::locate::solver::solve_mip;

/// Solve P-Median facility location problem
pub fn solve(cost_matrix: RMatrix<f64>, weights: &[f64], n_facilities: usize) -> List {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();
    let p = n_facilities;

    // Variables:
    // y[j] = 1 if facility j is selected (j = 0..n_fac-1)
    // x[i][j] = 1 if demand i is served by facility j (flattened: n_fac..n_fac+n_demand*n_fac-1)
    let n_y = n_fac;
    let n_x = n_demand * n_fac;
    let n_vars = n_y + n_x;

    // Objective: minimize sum_i sum_j (weight[i] * cost[i][j] * x[i][j])
    let mut obj_coeffs = vec![0.0; n_vars];
    for i in 0..n_demand {
        for j in 0..n_fac {
            let x_idx = n_y + i * n_fac + j;
            obj_coeffs[x_idx] = weights[i] * cost_matrix[[i, j]];
        }
    }

    // Constraints:
    let mut constraint_matrix: Vec<Vec<f64>> = Vec::new();
    let mut constraint_rhs: Vec<f64> = Vec::new();
    let mut constraint_sense: Vec<char> = Vec::new();

    // 1. sum_j y[j] = p (exactly p facilities)
    let mut c1 = vec![0.0; n_vars];
    for j in 0..n_fac {
        c1[j] = 1.0;
    }
    constraint_matrix.push(c1);
    constraint_rhs.push(p as f64);
    constraint_sense.push('=');

    // 2. sum_j x[i][j] = 1 for all i (each demand assigned to exactly one facility)
    for i in 0..n_demand {
        let mut row = vec![0.0; n_vars];
        for j in 0..n_fac {
            row[n_y + i * n_fac + j] = 1.0;
        }
        constraint_matrix.push(row);
        constraint_rhs.push(1.0);
        constraint_sense.push('=');
    }

    // 3. x[i][j] <= y[j] for all i,j (can only assign to open facility)
    for i in 0..n_demand {
        for j in 0..n_fac {
            let mut row = vec![0.0; n_vars];
            row[n_y + i * n_fac + j] = 1.0;  // x[i][j]
            row[j] = -1.0;                    // -y[j]
            constraint_matrix.push(row);
            constraint_rhs.push(0.0);
            constraint_sense.push('<');
        }
    }

    // Variable types: y[j] binary, x[i][j] continuous (LP relaxation often integral)
    let mut var_types = vec!['B'; n_vars];
    for i in n_y..n_vars {
        var_types[i] = 'C'; // Continuous for speed
    }

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
        false, // minimize
    );

    // Extract selected facilities (1-based for R)
    let selected: Vec<i32> = result
        .solution
        .iter()
        .take(n_fac)
        .enumerate()
        .filter(|(_, &v)| v > 0.5)
        .map(|(i, _)| (i + 1) as i32)
        .collect();

    // Extract assignments (which facility serves each demand)
    let assignments: Vec<i32> = (0..n_demand)
        .map(|i| {
            let mut best_j = 0;
            let mut best_val = 0.0;
            for j in 0..n_fac {
                let x_idx = n_y + i * n_fac + j;
                if result.solution[x_idx] > best_val {
                    best_val = result.solution[x_idx];
                    best_j = j;
                }
            }
            (best_j + 1) as i32 // 1-based
        })
        .collect();

    // Compute mean distance
    let total_weighted_dist: f64 = (0..n_demand)
        .map(|i| {
            let assigned_fac = (assignments[i] - 1) as usize;
            weights[i] * cost_matrix[[i, assigned_fac]]
        })
        .sum();
    let total_weight: f64 = weights.iter().sum();
    let mean_dist = total_weighted_dist / total_weight;

    let n_selected = selected.len() as i32;

    list!(
        selected = selected,
        assignments = assignments,
        n_selected = n_selected,
        objective = result.objective,
        mean_distance = mean_dist
    )
}
