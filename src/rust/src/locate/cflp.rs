//! Capacitated Facility Location Problem (CFLP)
//!
//! Minimize total weighted distance subject to facility capacity constraints.
//! Supports either fixed number of facilities (p) or facility opening costs.

use extendr_api::prelude::*;
use crate::locate::solver::solve_mip;

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

    // Variables:
    // y[j] = 1 if facility j is selected (j = 0..n_fac-1)
    // x[i][j] = fraction of demand i served by facility j (flattened)
    let n_y = n_fac;
    let n_x = n_demand * n_fac;
    let n_vars = n_y + n_x;

    // Objective: minimize sum_i sum_j (weight[i] * cost[i][j] * x[i][j])
    // Plus facility costs if provided
    let mut obj_coeffs = vec![0.0; n_vars];

    // Facility opening costs (if provided)
    if let Some(costs) = facility_costs {
        for j in 0..n_fac {
            obj_coeffs[j] = costs[j];
        }
    }

    // Transportation costs
    for i in 0..n_demand {
        for j in 0..n_fac {
            let x_idx = n_y + i * n_fac + j;
            obj_coeffs[x_idx] = weights[i] * cost_matrix[[i, j]];
        }
    }

    // Constraints
    let mut constraint_matrix: Vec<Vec<f64>> = Vec::new();
    let mut constraint_rhs: Vec<f64> = Vec::new();
    let mut constraint_sense: Vec<char> = Vec::new();

    // 1. If n_facilities specified: sum_j y[j] = p
    if n_facilities > 0 {
        let mut c1 = vec![0.0; n_vars];
        for j in 0..n_fac {
            c1[j] = 1.0;
        }
        constraint_matrix.push(c1);
        constraint_rhs.push(n_facilities as f64);
        constraint_sense.push('=');
    }

    // 2. sum_j x[i][j] = 1 for all i (each demand fully served)
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

    // 4. CAPACITY: sum_i (weight[i] * x[i][j]) <= capacity[j] * y[j] for all j
    // Rewritten: sum_i (weight[i] * x[i][j]) - capacity[j] * y[j] <= 0
    for j in 0..n_fac {
        let mut row = vec![0.0; n_vars];
        row[j] = -capacities[j];  // -capacity[j] * y[j]
        for i in 0..n_demand {
            row[n_y + i * n_fac + j] = weights[i];  // weight[i] * x[i][j]
        }
        constraint_matrix.push(row);
        constraint_rhs.push(0.0);
        constraint_sense.push('<');
    }

    // Variable types: y[j] binary, x[i][j] continuous
    let mut var_types = vec!['B'; n_vars];
    for i in n_y..n_vars {
        var_types[i] = 'C';
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

    // Check if solution is valid
    let status_str = format!("{:?}", result.status);
    if !status_str.contains("Optimal") {
        return list!(
            error = format!("Solver returned non-optimal status: {}", status_str),
            status = status_str
        );
    }

    // Extract selected facilities (1-based for R)
    let selected: Vec<i32> = result
        .solution
        .iter()
        .take(n_fac)
        .enumerate()
        .filter(|(_, &v)| v > 0.5)
        .map(|(i, _)| (i + 1) as i32)
        .collect();

    // Extract assignments - for CFLP, demand can be split across facilities
    // Return primary assignment (facility serving most of demand i)
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

    // Extract allocation fractions for split demands
    let mut allocation_fractions: Vec<f64> = Vec::with_capacity(n_demand * n_fac);
    for i in 0..n_demand {
        for j in 0..n_fac {
            let x_idx = n_y + i * n_fac + j;
            allocation_fractions.push(result.solution[x_idx]);
        }
    }

    // Check if any demand is split (allocation < 1 to primary)
    let n_split: i32 = (0..n_demand)
        .filter(|&i| {
            let assigned_j = (assignments[i] - 1) as usize;
            let x_idx = n_y + i * n_fac + assigned_j;
            result.solution[x_idx] < 0.999
        })
        .count() as i32;

    // Compute utilization of each facility
    let mut utilizations: Vec<f64> = vec![0.0; n_fac];
    for j in 0..n_fac {
        if result.solution[j] > 0.5 {
            let mut total_assigned = 0.0;
            for i in 0..n_demand {
                let x_idx = n_y + i * n_fac + j;
                total_assigned += weights[i] * result.solution[x_idx];
            }
            utilizations[j] = total_assigned / capacities[j];
        }
    }

    // Compute mean distance (weighted)
    let total_weighted_dist: f64 = (0..n_demand)
        .map(|i| {
            let mut dist = 0.0;
            for j in 0..n_fac {
                let x_idx = n_y + i * n_fac + j;
                dist += result.solution[x_idx] * cost_matrix[[i, j]];
            }
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
        objective = result.objective,
        mean_distance = mean_dist,
        utilizations = utilizations,
        status = status_str
    )
}
