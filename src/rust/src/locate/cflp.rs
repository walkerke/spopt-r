//! Capacitated Facility Location Problem (CFLP)
//!
//! Minimize total weighted distance subject to facility capacity constraints.
//! Supports either fixed number of facilities (p) or facility opening costs.

use extendr_api::prelude::*;
use highs::{HighsModelStatus, Sense, RowProblem, Col};

/// Solve Capacitated Facility Location Problem
///
/// Two modes:
/// 1. If n_facilities > 0: Select exactly n_facilities (capacitated p-median)
/// 2. If n_facilities = 0 and facility_costs provided: Minimize total cost including fixed costs
pub fn solve(
    cost_matrix: RMatrix<f64>,
    weights: &[f64],
    capacities: &[f64],
    n_facilities: usize,
    facility_costs: Option<&[f64]>,
) -> List {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();

    // Validate inputs
    if capacities.len() != n_fac {
        return list!(
            error = "capacities length must equal number of facilities"
        );
    }

    // Check total capacity is sufficient
    let total_demand: f64 = weights.iter().sum();
    let total_capacity: f64 = capacities.iter().sum();

    if total_capacity < total_demand {
        return list!(
            error = format!(
                "Total capacity ({:.2}) is less than total demand ({:.2}). Problem is infeasible.",
                total_capacity, total_demand
            )
        );
    }

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Variables:
    // y[j] = 1 if facility j is selected (binary)
    let y_cols: Vec<Col> = (0..n_fac)
        .map(|j| {
            let obj_coeff = facility_costs.map_or(0.0, |c| c[j]);
            pb.add_integer_column(obj_coeff, 0.0..=1.0)
        })
        .collect();

    // x[i][j] = fraction of demand i served by facility j (continuous)
    let mut x_cols: Vec<Vec<Col>> = Vec::with_capacity(n_demand);
    for i in 0..n_demand {
        let row_cols: Vec<Col> = (0..n_fac)
            .map(|j| {
                let obj_coeff = weights[i] * cost_matrix[[i, j]];
                pb.add_column(obj_coeff, 0.0..=1.0)
            })
            .collect();
        x_cols.push(row_cols);
    }

    // Constraint 1: If n_facilities specified: sum_j y[j] = p
    if n_facilities > 0 {
        let terms: Vec<(Col, f64)> = y_cols.iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(n_facilities as f64..=n_facilities as f64, terms);
    }

    // Constraint 2: sum_j x[i][j] = 1 for all i (each demand fully served)
    for i in 0..n_demand {
        let terms: Vec<(Col, f64)> = x_cols[i].iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(1.0..=1.0, terms);
    }

    // Constraint 3: x[i][j] <= y[j] for all i,j (can only assign to open facility)
    // Sparse: only 2 non-zeros per constraint
    for i in 0..n_demand {
        for j in 0..n_fac {
            let terms = vec![(x_cols[i][j], 1.0), (y_cols[j], -1.0)];
            pb.add_row(..=0.0, terms);
        }
    }

    // Constraint 4: CAPACITY: sum_i (weight[i] * x[i][j]) <= capacity[j] * y[j] for all j
    // Rewritten: sum_i (weight[i] * x[i][j]) - capacity[j] * y[j] <= 0
    for j in 0..n_fac {
        let mut terms: Vec<(Col, f64)> = (0..n_demand)
            .map(|i| (x_cols[i][j], weights[i]))
            .collect();
        terms.push((y_cols[j], -capacities[j]));
        pb.add_row(..=0.0, terms);
    }

    // Solve
    let solved = pb.optimise(Sense::Minimise).solve();
    let status = solved.status();
    let status_str = format!("{:?}", status);

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

            // Extract assignments - return primary assignment (facility serving most of demand i)
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

            // Extract allocation fractions for split demands
            let mut allocation_fractions: Vec<f64> = Vec::with_capacity(n_demand * n_fac);
            for i in 0..n_demand {
                for j in 0..n_fac {
                    allocation_fractions.push(sol[x_cols[i][j]]);
                }
            }

            // Check if any demand is split (allocation < 1 to primary)
            let n_split: i32 = (0..n_demand)
                .filter(|&i| {
                    let assigned_j = (assignments[i] - 1) as usize;
                    sol[x_cols[i][assigned_j]] < 0.999
                })
                .count() as i32;

            // Compute utilization of each facility
            let mut utilizations: Vec<f64> = vec![0.0; n_fac];
            for j in 0..n_fac {
                if sol[y_cols[j]] > 0.5 {
                    let total_assigned: f64 = (0..n_demand)
                        .map(|i| weights[i] * sol[x_cols[i][j]])
                        .sum();
                    utilizations[j] = total_assigned / capacities[j];
                }
            }

            // Compute mean distance (weighted)
            let total_weighted_dist: f64 = (0..n_demand)
                .map(|i| {
                    let dist: f64 = (0..n_fac)
                        .map(|j| sol[x_cols[i][j]] * cost_matrix[[i, j]])
                        .sum();
                    weights[i] * dist
                })
                .sum();
            let total_weight: f64 = weights.iter().sum();
            let mean_dist = total_weighted_dist / total_weight;

            let n_selected = selected.len() as i32;

            list!(
                selected = selected,
                assignments = assignments,
                allocation_matrix = allocation_fractions,
                n_selected = n_selected,
                n_split_demand = n_split,
                objective = solved.objective_value(),
                mean_distance = mean_dist,
                utilizations = utilizations,
                status = status_str
            )
        }
        _ => {
            list!(
                error = format!("Solver returned non-optimal status: {}", status_str),
                status = status_str
            )
        }
    }
}
