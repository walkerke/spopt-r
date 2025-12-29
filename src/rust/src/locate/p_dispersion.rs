//! P-Dispersion Problem
//!
//! Maximize the minimum distance between any two selected facilities.

use extendr_api::prelude::*;
use highs::{HighsModelStatus, Sense, RowProblem, Col};

/// Solve P-Dispersion problem
pub fn solve(distance_matrix: RMatrix<f64>, n_facilities: usize) -> List {
    let n = distance_matrix.nrows();
    let p = n_facilities;

    // Compute max distance for big-M and bounds
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

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Variables:
    // D = minimum inter-facility distance (continuous), objective = 1 (maximize)
    let d_col = pb.add_column(1.0, 0.0..=max_dist);

    // y[i] = 1 if facility i selected (binary)
    let y_cols: Vec<Col> = (0..n)
        .map(|_| pb.add_integer_column(0.0, 0.0..=1.0))
        .collect();

    // Constraint 1: sum_i y[i] = p
    {
        let terms: Vec<(Col, f64)> = y_cols.iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(p as f64..=p as f64, terms);
    }

    // Constraint 2: D <= d[i][j] + M * (2 - y[i] - y[j]) for all i < j
    // Rearranged: D + M*y[i] + M*y[j] <= d[i][j] + 2M
    // Each constraint has exactly 3 non-zeros (sparse!)
    for i in 0..n {
        for j in (i + 1)..n {
            let d_ij = distance_matrix[[i, j]];
            let rhs = d_ij + 2.0 * big_m;
            let terms = vec![
                (d_col, 1.0),
                (y_cols[i], big_m),
                (y_cols[j], big_m),
            ];
            pb.add_row(..=rhs, terms);
        }
    }

    // Solve (maximize)
    let solved = pb.optimise(Sense::Maximise).solve();
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

            let min_distance = sol[d_col];
            let n_selected = selected.len() as i32;

            list!(
                selected = selected,
                n_selected = n_selected,
                objective = solved.objective_value(),
                min_distance = min_distance
            )
        }
        _ => {
            eprintln!("P-Dispersion solver returned non-optimal status: {:?}", status);
            list!(
                selected = Vec::<i32>::new(),
                n_selected = 0i32,
                objective = f64::NAN,
                min_distance = f64::NAN
            )
        }
    }
}
