//! P-Center Problem
//!
//! Minimize the maximum distance from any demand point to its nearest facility.

use extendr_api::prelude::*;
use highs::{HighsModelStatus, Sense, RowProblem, Col};

/// Solve P-Center facility location problem
pub fn solve(cost_matrix: RMatrix<f64>, n_facilities: usize) -> List {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();
    let p = n_facilities;

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

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Variables:
    // W = maximum distance (continuous), objective coefficient = 1
    let w_col = pb.add_column(1.0, 0.0..=(max_dist * 2.0));

    // y[j] = 1 if facility j selected (binary)
    let y_cols: Vec<Col> = (0..n_fac)
        .map(|_| pb.add_integer_column(0.0, 0.0..=1.0))
        .collect();

    // x[i][j] = 1 if demand i assigned to facility j (continuous for LP relaxation)
    let mut x_cols: Vec<Vec<Col>> = Vec::with_capacity(n_demand);
    for _ in 0..n_demand {
        let row_cols: Vec<Col> = (0..n_fac)
            .map(|_| pb.add_column(0.0, 0.0..=1.0))
            .collect();
        x_cols.push(row_cols);
    }

    // Constraint 1: sum_j y[j] = p
    {
        let terms: Vec<(Col, f64)> = y_cols.iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(p as f64..=p as f64, terms);
    }

    // Constraint 2: sum_j x[i][j] = 1 for all i
    for i in 0..n_demand {
        let terms: Vec<(Col, f64)> = x_cols[i].iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(1.0..=1.0, terms);
    }

    // Constraint 3: x[i][j] <= y[j] for all i,j
    // Sparse: only 2 non-zeros per constraint
    for i in 0..n_demand {
        for j in 0..n_fac {
            let terms = vec![(x_cols[i][j], 1.0), (y_cols[j], -1.0)];
            pb.add_row(..=0.0, terms);
        }
    }

    // Constraint 4: W >= sum_j (d[i][j] * x[i][j]) for all i (minmax constraint)
    // Rewritten: sum_j (d[i][j] * x[i][j]) - W <= 0
    for i in 0..n_demand {
        let mut terms: Vec<(Col, f64)> = x_cols[i]
            .iter()
            .enumerate()
            .map(|(j, &c)| (c, cost_matrix[[i, j]]))
            .collect();
        terms.push((w_col, -1.0));
        pb.add_row(..=0.0, terms);
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

            // Extract assignments
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
                    (best_j + 1) as i32
                })
                .collect();

            let max_distance = sol[w_col];
            let n_selected = selected.len() as i32;

            list!(
                selected = selected,
                assignments = assignments,
                n_selected = n_selected,
                objective = solved.objective_value(),
                max_distance = max_distance
            )
        }
        _ => {
            eprintln!("P-Center solver returned non-optimal status: {:?}", status);
            list!(
                selected = Vec::<i32>::new(),
                assignments = vec![0i32; n_demand],
                n_selected = 0i32,
                objective = f64::NAN,
                max_distance = f64::NAN
            )
        }
    }
}
