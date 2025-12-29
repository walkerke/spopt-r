//! Coverage location problems: LSCP and MCLP
//!
//! LSCP: Location Set Covering Problem - minimize facilities to cover all demand
//! MCLP: Maximum Coverage Location Problem - maximize coverage with p facilities

use extendr_api::prelude::*;
use highs::{HighsModelStatus, Sense, RowProblem, Col};

/// Solve LSCP (Location Set Covering Problem)
///
/// Minimize the number of facilities such that all demand points are covered
/// within the service radius.
pub fn solve_lscp(cost_matrix: RMatrix<f64>, service_radius: f64) -> List {
    let n_demand = cost_matrix.nrows();
    let n_facilities = cost_matrix.ncols();

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Variables: y[j] = 1 if facility j is selected (binary)
    // Objective: minimize sum of facilities
    let y_cols: Vec<Col> = (0..n_facilities)
        .map(|_| pb.add_integer_column(1.0, 0.0..=1.0))  // obj coeff = 1
        .collect();

    // Build coverage info and constraints
    // Each demand must be covered by at least one facility: sum_j (coverage[i][j] * y[j]) >= 1
    let mut uncoverable = 0i32;

    for i in 0..n_demand {
        // Find facilities that can cover this demand (sparse!)
        let covering_facilities: Vec<(Col, f64)> = (0..n_facilities)
            .filter(|&j| cost_matrix[[i, j]] <= service_radius)
            .map(|j| (y_cols[j], 1.0))
            .collect();

        if covering_facilities.is_empty() {
            uncoverable += 1;
        } else {
            pb.add_row(1.0.., covering_facilities);
        }
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

            // Compute coverage statistics
            let n_selected = selected.len() as i32;
            let total_demand = n_demand as i32;

            // Count covered demand
            let covered_demand: i32 = (0..n_demand)
                .filter(|&i| {
                    (0..n_facilities).any(|j| {
                        cost_matrix[[i, j]] <= service_radius && sol[y_cols[j]] > 0.5
                    })
                })
                .count() as i32;

            let coverage_pct = (covered_demand as f64 / total_demand as f64) * 100.0;

            list!(
                selected = selected,
                n_selected = n_selected,
                objective = solved.objective_value(),
                covered_demand = covered_demand,
                total_demand = total_demand,
                coverage_pct = coverage_pct,
                status = status_str,
                uncoverable_demand = uncoverable
            )
        }
        _ => {
            list!(
                selected = Vec::<i32>::new(),
                n_selected = 0i32,
                objective = f64::NAN,
                covered_demand = 0i32,
                total_demand = n_demand as i32,
                coverage_pct = 0.0,
                status = status_str,
                uncoverable_demand = uncoverable
            )
        }
    }
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

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Variables:
    // y[j] = 1 if facility j is selected (binary)
    let y_cols: Vec<Col> = (0..n_facilities)
        .map(|_| pb.add_integer_column(0.0, 0.0..=1.0))  // no obj coeff for y
        .collect();

    // z[i] = 1 if demand i is covered (binary)
    // Objective: maximize sum of weighted coverage
    let z_cols: Vec<Col> = (0..n_demand)
        .map(|i| pb.add_integer_column(weights[i], 0.0..=1.0))
        .collect();

    // Constraint 1: sum_j y[j] = p (exactly p facilities)
    {
        let terms: Vec<(Col, f64)> = y_cols.iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(p as f64..=p as f64, terms);
    }

    // Constraint 2: z[i] <= sum_j (coverage[i][j] * y[j]) for all i
    // Can only be covered if a covering facility is open
    // Rewritten: z[i] - sum_j (coverage[i][j] * y[j]) <= 0
    for i in 0..n_demand {
        let mut terms: Vec<(Col, f64)> = vec![(z_cols[i], 1.0)];

        // Add negative coefficients for covering facilities (sparse!)
        for j in 0..n_facilities {
            if cost_matrix[[i, j]] <= service_radius {
                terms.push((y_cols[j], -1.0));
            }
        }

        pb.add_row(..=0.0, terms);
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

            // Coverage statistics
            let covered_weight: f64 = z_cols
                .iter()
                .zip(weights.iter())
                .filter(|(&c, _)| sol[c] > 0.5)
                .map(|(_, &w)| w)
                .sum();

            let total_weight: f64 = weights.iter().sum();
            let n_selected = selected.len() as i32;
            let coverage_pct = (covered_weight / total_weight) * 100.0;

            list!(
                selected = selected,
                n_selected = n_selected,
                objective = solved.objective_value(),
                covered_weight = covered_weight,
                total_weight = total_weight,
                coverage_pct = coverage_pct
            )
        }
        _ => {
            eprintln!("MCLP solver returned non-optimal status: {:?}", status);
            list!(
                selected = Vec::<i32>::new(),
                n_selected = 0i32,
                objective = f64::NAN,
                covered_weight = 0.0,
                total_weight = weights.iter().sum::<f64>(),
                coverage_pct = 0.0
            )
        }
    }
}
