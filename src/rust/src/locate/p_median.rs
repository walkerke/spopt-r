//! P-Median Problem
//!
//! Minimize total weighted distance by locating exactly p facilities.

use extendr_api::prelude::*;
use highs::{HighsModelStatus, Sense, RowProblem, Col};

/// Solve P-Median facility location problem
pub fn solve(cost_matrix: RMatrix<f64>, weights: &[f64], n_facilities: usize) -> List {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();
    let p = n_facilities;

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Variables:
    // y[j] = 1 if facility j is selected (j = 0..n_fac-1)
    // x[i][j] = 1 if demand i is served by facility j

    // Add y variables (binary facility selection)
    let y_cols: Vec<Col> = (0..n_fac)
        .map(|_| pb.add_integer_column(0.0, 0.0..=1.0))  // obj coeff = 0 for y
        .collect();

    // Add x variables (continuous assignment) with objective coefficients
    let mut x_cols: Vec<Vec<Col>> = Vec::with_capacity(n_demand);
    for i in 0..n_demand {
        let row_cols: Vec<Col> = (0..n_fac)
            .map(|j| {
                let obj_coeff = weights[i] * cost_matrix[[i, j]];
                pb.add_column(obj_coeff, 0.0..=1.0)  // continuous
            })
            .collect();
        x_cols.push(row_cols);
    }

    // Constraint 1: sum_j y[j] = p (exactly p facilities)
    {
        let terms: Vec<(Col, f64)> = y_cols.iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(p as f64..=p as f64, terms);
    }

    // Constraint 2: sum_j x[i][j] = 1 for all i (each demand assigned to exactly one facility)
    for i in 0..n_demand {
        let terms: Vec<(Col, f64)> = x_cols[i].iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(1.0..=1.0, terms);
    }

    // Constraint 3: x[i][j] <= y[j] for all i,j (can only assign to open facility)
    // Sparse: only 2 non-zeros per constraint!
    for i in 0..n_demand {
        for j in 0..n_fac {
            let terms = vec![(x_cols[i][j], 1.0), (y_cols[j], -1.0)];
            pb.add_row(..=0.0, terms);
        }
    }

    // Solve
    let solved = pb.optimise(Sense::Minimise).solve();
    let status = solved.status();

    match status {
        HighsModelStatus::Optimal | HighsModelStatus::ModelEmpty => {
            let sol = solved.get_solution();

            // Extract selected facilities (1-based for R)
            let selected: Vec<i32> = y_cols
                .iter()
                .enumerate()
                .filter(|(_, &c)| sol[c] > 0.5)
                .map(|(i, _)| (i + 1) as i32)
                .collect();

            // Extract assignments (which facility serves each demand)
            let assignments: Vec<i32> = (0..n_demand)
                .map(|i| {
                    let mut best_j = 0;
                    let mut best_val = 0.0;
                    for j in 0..n_fac {
                        let val = sol[x_cols[i][j]];
                        if val > best_val {
                            best_val = val;
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
                objective = solved.objective_value(),
                mean_distance = mean_dist
            )
        }
        _ => {
            eprintln!("P-Median solver returned non-optimal status: {:?}", status);
            list!(
                selected = Vec::<i32>::new(),
                assignments = vec![0i32; n_demand],
                n_selected = 0i32,
                objective = f64::NAN,
                mean_distance = f64::NAN
            )
        }
    }
}
