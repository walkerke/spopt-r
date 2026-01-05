//! P-Center Problem
//!
//! Minimize the maximum distance from any demand point to its nearest facility.
//!
//! Two algorithms are provided:
//! - `solve_mip`: Direct MIP formulation (minimax objective)
//! - `solve_binary_search`: Binary search over distances with set covering subproblems

use extendr_api::prelude::*;
use highs::{HighsModelStatus, Sense, RowProblem, Col};
use std::collections::BTreeSet;

/// Solve P-Center using the specified method
pub fn solve(cost_matrix: RMatrix<f64>, n_facilities: usize, method: &str) -> List {
    match method {
        "binary_search" => solve_binary_search(cost_matrix, n_facilities),
        "mip" => solve_mip(cost_matrix, n_facilities),
        _ => solve_binary_search(cost_matrix, n_facilities), // default
    }
}

/// Solve P-Center using binary search over distances with set covering subproblems
///
/// This approach is typically much faster than the direct MIP formulation because
/// it converts the difficult minimax objective into simpler feasibility problems.
fn solve_binary_search(cost_matrix: RMatrix<f64>, n_facilities: usize) -> List {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();
    let p = n_facilities;

    // Collect all unique distances and sort them
    let mut unique_distances: BTreeSet<OrderedFloat> = BTreeSet::new();
    for i in 0..n_demand {
        for j in 0..n_fac {
            let d = cost_matrix[[i, j]];
            if d.is_finite() {
                unique_distances.insert(OrderedFloat(d));
            }
        }
    }

    let distances: Vec<f64> = unique_distances.into_iter().map(|of| of.0).collect();

    if distances.is_empty() {
        return list!(
            selected = Vec::<i32>::new(),
            assignments = vec![0i32; n_demand],
            n_selected = 0i32,
            objective = f64::NAN,
            max_distance = f64::NAN
        );
    }

    // Binary search for minimum feasible distance
    let mut lo = 0;
    let mut hi = distances.len() - 1;
    let mut best_result: Option<(Vec<usize>, f64)> = None;

    while lo <= hi {
        let mid = lo + (hi - lo) / 2;
        let threshold = distances[mid];

        if let Some(selected) = can_cover_with_p_facilities(&cost_matrix, threshold, p) {
            // Feasible - try smaller distance
            best_result = Some((selected, threshold));
            if mid == 0 {
                break;
            }
            hi = mid - 1;
        } else {
            // Infeasible - need larger distance
            lo = mid + 1;
        }
    }

    match best_result {
        Some((selected_indices, max_dist)) => {
            // Compute assignments based on selected facilities
            let assignments: Vec<i32> = (0..n_demand)
                .map(|i| {
                    let mut best_j = selected_indices[0];
                    let mut best_dist = cost_matrix[[i, best_j]];
                    for &j in &selected_indices[1..] {
                        let d = cost_matrix[[i, j]];
                        if d < best_dist {
                            best_dist = d;
                            best_j = j;
                        }
                    }
                    (best_j + 1) as i32 // 1-based for R
                })
                .collect();

            let selected: Vec<i32> = selected_indices.iter().map(|&j| (j + 1) as i32).collect();
            let n_selected = selected.len() as i32;

            list!(
                selected = selected,
                assignments = assignments,
                n_selected = n_selected,
                objective = max_dist,
                max_distance = max_dist
            )
        }
        None => {
            // Should not happen if there's at least one facility
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

/// Check if all demand can be covered with exactly p facilities at distance <= threshold
/// Returns Some(selected_facilities) if feasible, None otherwise
fn can_cover_with_p_facilities(
    cost_matrix: &RMatrix<f64>,
    threshold: f64,
    p: usize,
) -> Option<Vec<usize>> {
    let n_demand = cost_matrix.nrows();
    let n_fac = cost_matrix.ncols();

    // Build coverage matrix: which facilities can cover each demand point?
    // coverage[i] = list of facilities that can cover demand i
    let mut coverage: Vec<Vec<usize>> = vec![Vec::new(); n_demand];
    let mut any_can_cover: Vec<bool> = vec![false; n_fac]; // track which facilities can cover anyone

    for i in 0..n_demand {
        for j in 0..n_fac {
            if cost_matrix[[i, j]] <= threshold {
                coverage[i].push(j);
                any_can_cover[j] = true;
            }
        }
        // If any demand point can't be covered, infeasible
        if coverage[i].is_empty() {
            return None;
        }
    }

    // Solve set covering: can we cover all demand with exactly p facilities?
    // Use ILP for exact solution
    let mut pb = RowProblem::new();

    // y[j] = 1 if facility j selected (objective = 0, just feasibility)
    let y_cols: Vec<Col> = (0..n_fac)
        .map(|_| pb.add_integer_column(0.0, 0.0..=1.0))
        .collect();

    // Constraint: exactly p facilities
    {
        let terms: Vec<(Col, f64)> = y_cols.iter().map(|&c| (c, 1.0)).collect();
        pb.add_row(p as f64..=p as f64, terms);
    }

    // Constraint: each demand point must be covered by at least one selected facility
    for i in 0..n_demand {
        let terms: Vec<(Col, f64)> = coverage[i].iter().map(|&j| (y_cols[j], 1.0)).collect();
        pb.add_row(1.0.., terms);
    }

    // Solve (just checking feasibility since objective is 0)
    let solved = pb.optimise(Sense::Minimise).solve();
    let status = solved.status();

    match status {
        HighsModelStatus::Optimal | HighsModelStatus::ModelEmpty => {
            let sol = solved.get_solution();
            let selected: Vec<usize> = y_cols
                .iter()
                .enumerate()
                .filter(|(_, &c)| sol[c] > 0.5)
                .map(|(j, _)| j)
                .collect();

            if selected.len() == p {
                Some(selected)
            } else {
                None
            }
        }
        _ => None, // Infeasible
    }
}

/// Wrapper type for f64 to enable BTreeSet ordering
#[derive(Clone, Copy)]
struct OrderedFloat(f64);

impl PartialEq for OrderedFloat {
    fn eq(&self, other: &Self) -> bool {
        self.0.to_bits() == other.0.to_bits()
    }
}

impl Eq for OrderedFloat {}

impl PartialOrd for OrderedFloat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedFloat {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.partial_cmp(&other.0).unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Solve P-Center using direct MIP formulation (original method)
fn solve_mip(cost_matrix: RMatrix<f64>, n_facilities: usize) -> List {
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
