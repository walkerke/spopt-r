//! HiGHS MIP solver interface
//!
//! Provides a wrapper around the HiGHS solver for facility location problems.

use highs::{HighsModelStatus, Sense, RowProblem, Col};

/// Result of solving a MIP problem
#[derive(Debug)]
pub struct SolveResult {
    pub status: HighsModelStatus,
    pub objective: f64,
    pub solution: Vec<f64>,
}

/// Solve a MIP problem using HiGHS
///
/// This is a low-level interface. High-level functions like solve_p_median
/// construct the appropriate model.
pub fn solve_mip(
    obj_coeffs: &[f64],
    constraint_matrix: &[Vec<f64>],
    constraint_rhs: &[f64],
    constraint_sense: &[char], // '<', '>', '='
    var_types: &[char],        // 'B' for binary, 'C' for continuous, 'I' for integer
    lb: &[f64],
    ub: &[f64],
    maximize: bool,
) -> SolveResult {
    let n_vars = obj_coeffs.len();

    // Create row-based problem
    let mut pb = RowProblem::new();

    // Columns (variables) with their indices
    let mut cols: Vec<Col> = Vec::with_capacity(n_vars);

    // Add variables with integrality specified at creation time
    // In highs 1.12+, integrality must be set when adding the column
    for i in 0..n_vars {
        // Use original coefficients - the sense (Maximise/Minimise) handles direction
        let is_integer = var_types[i] == 'B' || var_types[i] == 'I';
        let col = pb.add_column_with_integrality(obj_coeffs[i], lb[i]..=ub[i], is_integer);
        cols.push(col);
    }

    // Add constraints
    for (row_idx, row) in constraint_matrix.iter().enumerate() {
        // Build terms for this row
        let terms: Vec<(Col, f64)> = row
            .iter()
            .enumerate()
            .filter(|(_, &coef)| coef.abs() > 1e-10)
            .map(|(col_idx, &coef)| (cols[col_idx], coef))
            .collect();

        let rhs = constraint_rhs[row_idx];
        match constraint_sense[row_idx] {
            '<' | 'L' => {
                pb.add_row(..=rhs, terms);
            }
            '>' | 'G' => {
                pb.add_row(rhs.., terms);
            }
            '=' | 'E' => {
                pb.add_row(rhs..=rhs, terms);
            }
            _ => {
                pb.add_row(..=rhs, terms);
            }
        }
    }

    // Solve - optimise() returns a Model, solve() returns a SolvedModel
    let sense = if maximize { Sense::Maximise } else { Sense::Minimise };
    let solved = pb.optimise(sense).solve();

    // Extract results
    let status = solved.status();
    let _status_str = format!("{:?}", status);

    match status {
        HighsModelStatus::Optimal | HighsModelStatus::ModelEmpty => {
            // Get solution and extract variable values
            let sol = solved.get_solution();
            let solution: Vec<f64> = cols.iter().map(|&c| sol[c]).collect();
            SolveResult {
                status,
                objective: solved.objective_value(),
                solution,
            }
        }
        _ => {
            // Log the status for debugging
            eprintln!("HiGHS solver returned non-optimal status: {:?}", status);
            SolveResult {
                status,
                objective: f64::NAN,
                solution: vec![0.0; n_vars],
            }
        }
    }
}
