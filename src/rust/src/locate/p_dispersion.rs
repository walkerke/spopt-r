//! P-Dispersion Problem
//!
//! Maximize the minimum distance between any two selected facilities.

use extendr_api::prelude::*;
use crate::locate::solver::solve_mip;

/// Solve P-Dispersion problem
pub fn solve(distance_matrix: RMatrix<f64>, n_facilities: usize) -> List {
    let n = distance_matrix.nrows();
    let p = n_facilities;

    // Variables:
    // D = minimum inter-facility distance (continuous)
    // y[i] = 1 if facility i selected
    let n_vars = 1 + n;
    let d_idx = 0;
    let y_start = 1;

    // Objective: maximize D
    let mut obj_coeffs = vec![0.0; n_vars];
    obj_coeffs[d_idx] = 1.0;

    let mut constraint_matrix: Vec<Vec<f64>> = Vec::new();
    let mut constraint_rhs: Vec<f64> = Vec::new();
    let mut constraint_sense: Vec<char> = Vec::new();

    // 1. sum_i y[i] = p
    let mut c1 = vec![0.0; n_vars];
    for i in 0..n {
        c1[y_start + i] = 1.0;
    }
    constraint_matrix.push(c1);
    constraint_rhs.push(p as f64);
    constraint_sense.push('=');

    // Compute max distance for big-M
    let mut max_dist: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let d = distance_matrix[[i, j]];
            if d > max_dist {
                max_dist = d;
            }
        }
    }
    let big_m = max_dist + 1.0;

    // 2. D <= d[i][j] + M * (2 - y[i] - y[j]) for all i < j
    // Rearranged: D + M*y[i] + M*y[j] <= d[i][j] + 2M
    for i in 0..n {
        for j in (i + 1)..n {
            let mut row = vec![0.0; n_vars];
            row[d_idx] = 1.0;           // D
            row[y_start + i] = big_m;   // M * y[i]
            row[y_start + j] = big_m;   // M * y[j]
            constraint_matrix.push(row);
            let d_ij = distance_matrix[[i, j]];
            constraint_rhs.push(d_ij + 2.0 * big_m);
            constraint_sense.push('<');
        }
    }

    // Variable types
    let mut var_types = vec!['C'; n_vars];
    for i in 0..n {
        var_types[y_start + i] = 'B';
    }

    let lb = vec![0.0; n_vars];
    let mut ub = vec![max_dist; n_vars];
    for i in 0..n {
        ub[y_start + i] = 1.0;
    }

    // Solve (maximize)
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

    // Extract results
    let selected: Vec<i32> = result
        .solution
        .iter()
        .skip(y_start)
        .take(n)
        .enumerate()
        .filter(|(_, &v)| v > 0.5)
        .map(|(i, _)| (i + 1) as i32)
        .collect();

    let min_distance = result.solution[d_idx];
    let n_selected = selected.len() as i32;

    list!(
        selected = selected,
        n_selected = n_selected,
        objective = result.objective,
        min_distance = min_distance
    )
}
