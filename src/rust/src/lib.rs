use extendr_api::prelude::*;

mod distance;
mod graph;
mod region;
mod locate;

/// Compute Euclidean distance matrix between two sets of points
///
/// @param x1 X coordinates of first set of points
/// @param y1 Y coordinates of first set of points
/// @param x2 X coordinates of second set of points
/// @param y2 Y coordinates of second set of points
/// @return Distance matrix (n1 x n2)
/// @export
#[extendr]
fn rust_distance_matrix_euclidean(
    x1: Vec<f64>,
    y1: Vec<f64>,
    x2: Vec<f64>,
    y2: Vec<f64>,
) -> RMatrix<f64> {
    distance::euclidean_matrix(&x1, &y1, &x2, &y2)
}

/// Compute Manhattan distance matrix between two sets of points
///
/// @param x1 X coordinates of first set of points
/// @param y1 Y coordinates of first set of points
/// @param x2 X coordinates of second set of points
/// @param y2 Y coordinates of second set of points
/// @return Distance matrix (n1 x n2)
/// @export
#[extendr]
fn rust_distance_matrix_manhattan(
    x1: Vec<f64>,
    y1: Vec<f64>,
    x2: Vec<f64>,
    y2: Vec<f64>,
) -> RMatrix<f64> {
    distance::manhattan_matrix(&x1, &y1, &x2, &y2)
}

/// Compute minimum spanning tree from adjacency matrix
///
/// @param i Row indices (0-based) of adjacency matrix non-zero entries
/// @param j Column indices (0-based) of adjacency matrix non-zero entries
/// @param weights Edge weights (distances/dissimilarities)
/// @param n Number of nodes
/// @return List with MST edges (from, to, weight)
/// @export
#[extendr]
fn rust_mst(i: Vec<i32>, j: Vec<i32>, weights: Vec<f64>, n: i32) -> List {
    graph::minimum_spanning_tree(&i, &j, &weights, n as usize)
}

/// Check if a graph is connected
///
/// @param i Row indices of adjacency matrix
/// @param j Column indices of adjacency matrix
/// @param n Number of nodes
/// @return TRUE if connected, FALSE otherwise
/// @export
#[extendr]
fn rust_is_connected(i: Vec<i32>, j: Vec<i32>, n: i32) -> bool {
    graph::is_connected(&i, &j, n as usize)
}

/// Find connected components
///
/// @param i Row indices of adjacency matrix
/// @param j Column indices of adjacency matrix
/// @param n Number of nodes
/// @return Vector of component labels (0-based)
/// @export
#[extendr]
fn rust_connected_components(i: Vec<i32>, j: Vec<i32>, n: i32) -> Vec<i32> {
    graph::connected_components(&i, &j, n as usize)
}

/// Solve SKATER regionalization
///
/// @param attrs Attribute matrix (n x p)
/// @param adj_i Row indices of adjacency
/// @param adj_j Column indices of adjacency
/// @param n_regions Number of regions to create
/// @param floor_var Optional floor variable values
/// @param floor_value Minimum floor value per region
/// @param seed Random seed
/// @return Vector of region labels (1-based)
/// @export
#[extendr]
fn rust_skater(
    attrs: RMatrix<f64>,
    adj_i: Vec<i32>,
    adj_j: Vec<i32>,
    n_regions: i32,
    floor_var: Nullable<Vec<f64>>,
    floor_value: f64,
    seed: Nullable<i64>,
) -> Vec<i32> {
    region::skater::solve(
        attrs,
        &adj_i,
        &adj_j,
        n_regions as usize,
        floor_var.into_option(),
        floor_value,
        seed.into_option().map(|s| s as u64),
    )
}

/// Solve P-Median facility location problem
///
/// @param cost_matrix Cost/distance matrix (demand x facilities)
/// @param weights Demand weights
/// @param n_facilities Number of facilities to locate (p)
/// @return List with selected facilities and assignments
/// @export
#[extendr]
fn rust_p_median(
    cost_matrix: RMatrix<f64>,
    weights: Vec<f64>,
    n_facilities: i32,
) -> List {
    locate::p_median::solve(cost_matrix, &weights, n_facilities as usize)
}

/// Solve LSCP (Location Set Covering Problem)
///
/// @param cost_matrix Cost/distance matrix (demand x facilities)
/// @param service_radius Maximum service distance
/// @return List with selected facilities and coverage
/// @export
#[extendr]
fn rust_lscp(cost_matrix: RMatrix<f64>, service_radius: f64) -> List {
    locate::coverage::solve_lscp(cost_matrix, service_radius)
}

/// Solve MCLP (Maximum Coverage Location Problem)
///
/// @param cost_matrix Cost/distance matrix (demand x facilities)
/// @param weights Demand weights
/// @param service_radius Maximum service distance
/// @param n_facilities Number of facilities to locate
/// @return List with selected facilities and coverage
/// @export
#[extendr]
fn rust_mclp(
    cost_matrix: RMatrix<f64>,
    weights: Vec<f64>,
    service_radius: f64,
    n_facilities: i32,
) -> List {
    locate::coverage::solve_mclp(cost_matrix, &weights, service_radius, n_facilities as usize)
}

/// Solve P-Center facility location problem
///
/// @param cost_matrix Cost/distance matrix (demand x facilities)
/// @param n_facilities Number of facilities to locate
/// @return List with selected facilities, assignments, and max distance
/// @export
#[extendr]
fn rust_p_center(
    cost_matrix: RMatrix<f64>,
    n_facilities: i32,
) -> List {
    locate::p_center::solve(cost_matrix, n_facilities as usize)
}

/// Solve P-Dispersion facility location problem
///
/// @param distance_matrix Distance matrix between facilities
/// @param n_facilities Number of facilities to select
/// @return List with selected facilities and min distance
/// @export
#[extendr]
fn rust_p_dispersion(
    distance_matrix: RMatrix<f64>,
    n_facilities: i32,
) -> List {
    locate::p_dispersion::solve(distance_matrix, n_facilities as usize)
}

/// Solve Max-P regionalization problem
///
/// Maximize the number of regions such that each region satisfies a minimum
/// threshold constraint on a spatial extensive attribute.
///
/// @param attrs Attribute matrix (n x p) for computing within-region dissimilarity
/// @param threshold_var Values of the threshold variable (e.g., population)
/// @param threshold Minimum sum required per region
/// @param adj_i Row indices of adjacency (0-based)
/// @param adj_j Column indices of adjacency (0-based)
/// @param n_iterations Number of construction phase iterations
/// @param n_sa_iterations Number of simulated annealing iterations
/// @param cooling_rate SA cooling rate (e.g., 0.99)
/// @param tabu_length Tabu list length for SA
/// @param seed Random seed
/// @return List with labels (1-based), n_regions, and objective
/// @export
#[extendr]
fn rust_max_p(
    attrs: RMatrix<f64>,
    threshold_var: Vec<f64>,
    threshold: f64,
    adj_i: Vec<i32>,
    adj_j: Vec<i32>,
    n_iterations: i32,
    n_sa_iterations: i32,
    cooling_rate: f64,
    tabu_length: i32,
    seed: Nullable<i64>,
) -> List {
    region::maxp::solve(
        attrs,
        &threshold_var,
        threshold,
        &adj_i,
        &adj_j,
        n_iterations as usize,
        n_sa_iterations as usize,
        cooling_rate,
        tabu_length as usize,
        seed.into_option().map(|s| s as u64),
    )
}

// Macro to generate exports
extendr_module! {
    mod spopt;
    fn rust_distance_matrix_euclidean;
    fn rust_distance_matrix_manhattan;
    fn rust_mst;
    fn rust_is_connected;
    fn rust_connected_components;
    fn rust_skater;
    fn rust_max_p;
    fn rust_p_median;
    fn rust_lscp;
    fn rust_mclp;
    fn rust_p_center;
    fn rust_p_dispersion;
}
