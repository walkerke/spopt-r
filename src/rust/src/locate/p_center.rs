//! P-Center Problem
//!
//! Minimize the maximum distance from any demand point to its nearest facility.

use extendr_api::prelude::*;
use crate::locate::solver::solve_mip;

/// Solve P-Center facility location problem
pub fn solve(cost_matrix: RMatrix<f64>, n_facilities: usize) -> List {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();
    let p = n_facilities;

    // Variables:
    // W = maximum distance (continuous)
    // y[j] = 1 if facility j selected
    // x[i][j] = 1 if demand i assigned to facility j
    let n_vars = 1 + n_fac + n_demand * n_fac;
    let w_idx = 0;
    let y_start = 1;
    let x_start = 1 + n_fac;

    // Objective: minimize W
    let mut obj_coeffs = vec![0.0; n_vars];
    obj_coeffs[w_idx] = 1.0;

    let mut constraint_matrix: Vec<Vec<f64>> = Vec::new();
    let mut constraint_rhs: Vec<f64> = Vec::new();
    let mut constraint_sense: Vec<char> = Vec::new();

    // 1. sum_j y[j] = p
    let mut c1 = vec![0.0; n_vars];
    for j in 0..n_fac {
        c1[y_start + j] = 1.0;
    }
    constraint_matrix.push(c1);
    constraint_rhs.push(p as f64);
    constraint_sense.push('=');

    // 2. sum_j x[i][j] = 1 for all i
    for i in 0..n_demand {
        let mut row = vec![0.0; n_vars];
        for j in 0..n_fac {
            row[x_start + i * n_fac + j] = 1.0;
        }
        constraint_matrix.push(row);
        constraint_rhs.push(1.0);
        constraint_sense.push('=');
    }

    // 3. x[i][j] <= y[j]
    for i in 0..n_demand {
        for j in 0..n_fac {
            let mut row = vec![0.0; n_vars];
            row[x_start + i * n_fac + j] = 1.0;
            row[y_start + j] = -1.0;
            constraint_matrix.push(row);
            constraint_rhs.push(0.0);
            constraint_sense.push('<');
        }
    }

    // 4. W >= d[i][j] * x[i][j] for all i,j (minmax constraint)
    // sum_j (d[i][j] * x[i][j]) <= W for all i
    for i in 0..n_demand {
        let mut row = vec![0.0; n_vars];
        row[w_idx] = -1.0; // -W
        for j in 0..n_fac {
            row[x_start + i * n_fac + j] = cost_matrix[[i, j]];
        }
        constraint_matrix.push(row);
        constraint_rhs.push(0.0);
        constraint_sense.push('<');
    }

    // Compute max distance for upper bound
    let mut max_dist: f64 = 0.0;
    for i in 0..n_demand {
        for j in 0..n_fac {
            let d = cost_matrix[[i, j]];
            if d > max_dist {
                max_dist = d;
            }
        }
    }

    // Variable types
    let mut var_types = vec!['C'; n_vars];
    for j in 0..n_fac {
        var_types[y_start + j] = 'B'; // y binary
    }

    let lb = vec![0.0; n_vars];
    let mut ub = vec![1.0; n_vars];
    ub[w_idx] = max_dist * 2.0;

    // Solve
    let result = solve_mip(
        &obj_coeffs,
        &constraint_matrix,
        &constraint_rhs,
        &constraint_sense,
        &var_types,
        &lb,
        &ub,
        false,
    );

    // Extract results
    let selected: Vec<i32> = result
        .solution
        .iter()
        .skip(y_start)
        .take(n_fac)
        .enumerate()
        .filter(|(_, &v)| v > 0.5)
        .map(|(i, _)| (i + 1) as i32)
        .collect();

    let assignments: Vec<i32> = (0..n_demand)
        .map(|i| {
            let mut best_j = 0;
            let mut best_val = 0.0;
            for j in 0..n_fac {
                let x_idx = x_start + i * n_fac + j;
                if result.solution[x_idx] > best_val {
                    best_val = result.solution[x_idx];
                    best_j = j;
                }
            }
            (best_j + 1) as i32
        })
        .collect();

    let max_distance = result.solution[w_idx];
    let n_selected = selected.len() as i32;

    list!(
        selected = selected,
        assignments = assignments,
        n_selected = n_selected,
        objective = result.objective,
        max_distance = max_distance
    )
}
