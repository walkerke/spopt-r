use extendr_api::prelude::*;

mod distance;
mod graph;
mod region;
mod locate;
mod route;

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

/// Solve spatially-constrained Ward clustering
///
/// @param attrs Attribute matrix (n x p)
/// @param n_regions Number of regions to create
/// @param adj_i Row indices of adjacency (0-based)
/// @param adj_j Column indices of adjacency (0-based)
/// @return List with labels, n_regions, objective
/// @export
#[extendr]
fn rust_ward_constrained(
    attrs: RMatrix<f64>,
    n_regions: i32,
    adj_i: Vec<i32>,
    adj_j: Vec<i32>,
) -> List {
    region::ward::solve(
        attrs,
        n_regions as usize,
        &adj_i,
        &adj_j,
    )
}

/// Solve SPENC regionalization problem
///
/// Spatially-encouraged spectral clustering.
///
/// @param attrs Attribute matrix (n x p)
/// @param n_regions Number of regions to create
/// @param adj_i Row indices of adjacency (0-based)
/// @param adj_j Column indices of adjacency (0-based)
/// @param gamma RBF kernel parameter
/// @param seed Random seed
/// @return List with labels, n_regions, objective
/// @export
#[extendr]
fn rust_spenc(
    attrs: RMatrix<f64>,
    n_regions: i32,
    adj_i: Vec<i32>,
    adj_j: Vec<i32>,
    gamma: f64,
    seed: Nullable<i64>,
) -> List {
    region::spenc::solve(
        attrs,
        n_regions as usize,
        &adj_i,
        &adj_j,
        gamma,
        seed.into_option().map(|s| s as u64),
    )
}

/// Solve AZP regionalization problem
///
/// Automatic Zoning Procedure with basic, tabu, and SA variants.
///
/// @param attrs Attribute matrix (n x p)
/// @param n_regions Number of regions to create
/// @param adj_i Row indices of adjacency (0-based)
/// @param adj_j Column indices of adjacency (0-based)
/// @param method "basic", "tabu", or "sa"
/// @param max_iterations Maximum iterations
/// @param tabu_length Tabu list length (for tabu method)
/// @param cooling_rate SA cooling rate (for sa method)
/// @param initial_temperature SA initial temperature (for sa method)
/// @param seed Random seed
/// @return List with labels, n_regions, objective
/// @export
#[extendr]
fn rust_azp(
    attrs: RMatrix<f64>,
    n_regions: i32,
    adj_i: Vec<i32>,
    adj_j: Vec<i32>,
    method: &str,
    max_iterations: i32,
    tabu_length: i32,
    cooling_rate: f64,
    initial_temperature: f64,
    initial_regions: Nullable<Vec<i32>>,
    seed: Nullable<i64>,
) -> List {
    region::azp::solve(
        attrs,
        n_regions as usize,
        &adj_i,
        &adj_j,
        method,
        max_iterations as usize,
        tabu_length as usize,
        cooling_rate,
        initial_temperature,
        initial_regions.into_option(),
        seed.into_option().map(|s| s as u64),
    )
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
/// @param compact Whether to optimize for compactness
/// @param compact_weight Weight for compactness vs dissimilarity (0-1)
/// @param compact_metric Compactness metric to use, either "centroid" or "nmi" (normalized moment of inertia)
/// @param centroids_x X coordinates of unit centroids (for compactness)
/// @param centroids_y Y coordinates of unit centroids (for compactness)
/// @param areas Areas of units (for compactness)
/// @param moments Second areal moments of units (for compactness)
/// @return List with labels (1-based), n_regions, objective, and compactness
/// @noRd
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
    homogeneous: bool,
    compact: bool,
    compact_weight: f64,
    compact_metric: &str,
    centroids_x: Nullable<Vec<f64>>,
    centroids_y: Nullable<Vec<f64>>,
    areas: Nullable<Vec<f64>>,
    moments: Nullable<Vec<f64>>
) -> List {

    // Parse compact_metric string to enum
    let metric = match compact_metric.to_lowercase().as_str() {
        "nmi" => region::maxp::CompactnessMetric::NMI,
        "centroid" | "centroid dispersion" | "centroid-dispersion" | "centroid_dispersion" => {
            region::maxp::CompactnessMetric::CentroidDispersion
        }
        _ => {
            panic!("Invalid `compact_metric`: '{}'. Must be 'nmi' or 'centroid'", compact_metric);
        }
    };

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
        centroids_x.into_option(),
        centroids_y.into_option(),
        compact,
        compact_weight,
        homogeneous,
        metric,
        areas.into_option(),
        moments.into_option()
    )
}

/// Solve FRLM using greedy heuristic
///
/// @param n_candidates Number of candidate facility locations
/// @param path_candidates Flat array of candidate indices for each path
/// @param path_offsets Start index for each path in path_candidates
/// @param path_distances Distances to each candidate along paths
/// @param flow_volumes Volume of each flow
/// @param vehicle_range Maximum vehicle range
/// @param n_facilities Number of facilities to place
/// @return List with selected facilities and coverage info
/// @export
#[extendr]
fn rust_frlm_greedy(
    n_candidates: i32,
    path_candidates: Vec<i32>,
    path_offsets: Vec<i32>,
    path_distances: Vec<f64>,
    flow_volumes: Vec<f64>,
    vehicle_range: f64,
    n_facilities: i32,
) -> List {
    locate::frlm::solve_greedy(
        n_candidates as usize,
        &path_candidates,
        &path_offsets,
        &path_distances,
        &flow_volumes,
        vehicle_range,
        n_facilities as usize,
    )
}

/// Compute Huff Model probabilities
///
/// Computes probability surface based on distance decay and attractiveness.
/// Formula: P_ij = (A_j × D_ij^β) / Σ_k(A_k × D_ik^β)
///
/// @param cost_matrix Cost/distance matrix (demand x stores)
/// @param attractiveness Attractiveness values for each store (pre-computed with exponents)
/// @param distance_exponent Distance decay exponent (typically negative, e.g., -1.5)
/// @param sales_potential Optional sales potential for each demand point
/// @return List with probabilities, market shares, expected sales
/// @export
#[extendr]
fn rust_huff(
    cost_matrix: RMatrix<f64>,
    attractiveness: Vec<f64>,
    distance_exponent: f64,
    sales_potential: Nullable<Vec<f64>>,
) -> List {
    locate::huff::compute_probabilities(
        cost_matrix,
        &attractiveness,
        distance_exponent,
        sales_potential.into_option().as_deref(),
    )
}

/// Solve Traveling Salesman Problem (TSP)
///
/// Solve a closed tour, open route, or fixed-end path over a square
/// cost/distance matrix using nearest-neighbor construction with optional
/// 2-opt and or-opt local search. Optional time windows and service times
/// can be supplied for each stop.
///
/// @param cost_matrix Square cost/distance matrix (n x n)
/// @param start Start index (0-based)
/// @param end End index (0-based), or NULL for an open route
/// @param method Algorithm: "nn" (nearest-neighbor only) or "2-opt" (with local search)
/// @param earliest Optional earliest service times
/// @param latest Optional latest service times
/// @param service_time Optional service times at each stop
/// @return List with tour (1-based), total_cost, nn_cost, improvement_pct, iterations,
///   arrival_time, and departure_time
/// @export
#[extendr]
fn rust_tsp(
    cost_matrix: RMatrix<f64>,
    start: i32,
    end: Nullable<i32>,
    method: &str,
    earliest: Nullable<Vec<f64>>,
    latest: Nullable<Vec<f64>>,
    service_time: Nullable<Vec<f64>>,
) -> List {
    let n = cost_matrix.nrows();
    if cost_matrix.ncols() != n {
        extendr_api::throw_r_error(format!(
            "cost_matrix must be square, got {} x {}",
            n,
            cost_matrix.ncols()
        ));
    }
    if start < 0 || start as usize >= n {
        extendr_api::throw_r_error(format!(
            "start index {} is out of range [0, {})",
            start, n
        ));
    }

    let end_opt = end.into_option();
    if let Some(end_idx) = end_opt {
        if end_idx < 0 || end_idx as usize >= n {
            extendr_api::throw_r_error(format!(
                "end index {} is out of range [0, {})",
                end_idx, n
            ));
        }
    }

    let earliest_opt = earliest.into_option();
    let latest_opt = latest.into_option();
    let service_opt = service_time.into_option();

    if earliest_opt.is_some() != latest_opt.is_some() {
        extendr_api::throw_r_error("Both `earliest` and `latest` must be supplied together.");
    }
    if let Some(values) = earliest_opt.as_ref() {
        if values.len() != n {
            extendr_api::throw_r_error(format!(
                "`earliest` must have length {}, got {}",
                n,
                values.len()
            ));
        }
    }
    if let Some(values) = latest_opt.as_ref() {
        if values.len() != n {
            extendr_api::throw_r_error(format!(
                "`latest` must have length {}, got {}",
                n,
                values.len()
            ));
        }
    }
    if let Some(values) = service_opt.as_ref() {
        if values.len() != n {
            extendr_api::throw_r_error(format!(
                "`service_time` must have length {}, got {}",
                n,
                values.len()
            ));
        }
        if values.iter().any(|&value| value < 0.0 || !value.is_finite()) {
            extendr_api::throw_r_error("`service_time` must contain finite non-negative values.");
        }
    }
    if let (Some(earliest_values), Some(latest_values)) = (earliest_opt.as_ref(), latest_opt.as_ref()) {
        for idx in 0..n {
            if !earliest_values[idx].is_finite() || !latest_values[idx].is_finite() {
                extendr_api::throw_r_error("`earliest` and `latest` must contain finite values.");
            }
            if earliest_values[idx] > latest_values[idx] {
                extendr_api::throw_r_error(format!(
                    "`earliest[{idx}]` is greater than `latest[{idx}]`."
                ));
            }
        }
    }

    match route::tsp::solve(
        cost_matrix,
        start as usize,
        end_opt.map(|value| value as usize),
        method,
        earliest_opt.as_deref(),
        latest_opt.as_deref(),
        service_opt.as_deref(),
    ) {
        Ok(value) => value,
        Err(err) => extendr_api::throw_r_error(err.to_string()),
    }
}

/// Solve Capacitated Vehicle Routing Problem (CVRP)
///
/// Find minimum-cost routes for multiple vehicles, each with a capacity limit,
/// to serve all customers from a depot. Uses Clarke-Wright savings heuristic
/// with 2-opt, relocate, and swap improvement.
///
/// @param cost_matrix Square cost/distance matrix (n x n)
/// @param depot Depot index (0-based)
/// @param demands Demand at each location (depot demand should be 0)
/// @param capacity Vehicle capacity
/// @param max_vehicles Maximum number of vehicles (NULL for unlimited)
/// @param method Algorithm: "savings" (construction only) or "2-opt" (with improvement)
/// @param service_times Optional service time at each stop (NULL for zero)
/// @param max_route_time Optional maximum total time per route (NULL for unlimited)
/// @param balance_time Whether to run route-time balancing phase
/// @param earliest Optional earliest arrival times at each stop
/// @param latest Optional latest arrival times at each stop
/// @return List with vehicle assignments, visit orders, costs, and route details
/// @export
#[extendr]
fn rust_vrp(
    cost_matrix: RMatrix<f64>,
    depot: i32,
    demands: Vec<f64>,
    capacity: f64,
    max_vehicles: Nullable<i32>,
    method: &str,
    service_times: Nullable<Vec<f64>>,
    max_route_time: Nullable<f64>,
    balance_time: bool,
    earliest: Nullable<Vec<f64>>,
    latest: Nullable<Vec<f64>>,
) -> List {
    let n = cost_matrix.nrows();
    if depot < 0 || depot as usize >= n {
        extendr_api::throw_r_error(format!(
            "depot index {} is out of range [0, {})",
            depot, n
        ));
    }
    if demands.len() != n {
        extendr_api::throw_r_error(format!(
            "`demands` must have length {}, got {}",
            n, demands.len()
        ));
    }
    if capacity <= 0.0 {
        extendr_api::throw_r_error(format!(
            "capacity must be positive, got {}",
            capacity
        ));
    }
    let max_v = max_vehicles.into_option().map(|v| {
        if v < 1 {
            extendr_api::throw_r_error(format!("max_vehicles must be positive, got {}", v));
        }
        v as usize
    });

    let svc = service_times.into_option().unwrap_or_else(|| vec![0.0; n]);
    if svc.len() != n {
        extendr_api::throw_r_error(format!(
            "`service_times` must have length {}, got {}",
            n, svc.len()
        ));
    }
    if svc.iter().any(|&v| v < 0.0 || !v.is_finite()) {
        extendr_api::throw_r_error("`service_times` must contain finite non-negative values.");
    }

    let max_t = max_route_time.into_option();
    if let Some(t) = max_t {
        if t <= 0.0 || !t.is_finite() {
            extendr_api::throw_r_error(format!(
                "`max_route_time` must be positive and finite, got {}",
                t
            ));
        }
    }

    let earliest_opt = earliest.into_option();
    let latest_opt = latest.into_option();

    if earliest_opt.is_some() != latest_opt.is_some() {
        extendr_api::throw_r_error("Both `earliest` and `latest` must be supplied together.");
    }
    if let Some(values) = earliest_opt.as_ref() {
        if values.len() != n {
            extendr_api::throw_r_error(format!(
                "`earliest` must have length {}, got {}", n, values.len()
            ));
        }
        if values.iter().any(|v| !v.is_finite()) {
            extendr_api::throw_r_error("`earliest` must contain finite values.");
        }
    }
    if let Some(values) = latest_opt.as_ref() {
        if values.len() != n {
            extendr_api::throw_r_error(format!(
                "`latest` must have length {}, got {}", n, values.len()
            ));
        }
        if values.iter().any(|v| !v.is_finite()) {
            extendr_api::throw_r_error("`latest` must contain finite values.");
        }
    }
    if let (Some(ew), Some(lw)) = (earliest_opt.as_ref(), latest_opt.as_ref()) {
        for idx in 0..n {
            if ew[idx] > lw[idx] + 1e-10 {
                extendr_api::throw_r_error(format!(
                    "`earliest[{}]` ({}) is greater than `latest[{}]` ({}).",
                    idx, ew[idx], idx, lw[idx]
                ));
            }
        }
    }

    route::vrp::solve(
        cost_matrix,
        depot as usize,
        &demands,
        capacity,
        max_v,
        method,
        &svc,
        max_t,
        balance_time,
        earliest_opt.as_deref(),
        latest_opt.as_deref(),
    )
}

/// Least-cost corridor routing on a raster grid
///
/// @param values Flattened cell values (row-major), NaN = impassable
/// @param n_rows Number of rows in the raster
/// @param n_cols Number of columns in the raster
/// @param cell_width Cell width in CRS units
/// @param cell_height Cell height in CRS units
/// @param origin_cell 0-based cell index of origin
/// @param dest_cell 0-based cell index of destination
/// @param neighbours Cell connectivity: 4, 8, or 16
/// @param method Routing algorithm: "dijkstra", "bidirectional", or "astar"
/// @return List with path_cells, total_cost, solve_time_ms, graph_build_time_ms, n_edges
/// @export
#[extendr]
fn rust_corridor(
    values: Vec<f64>,
    n_rows: i32,
    n_cols: i32,
    cell_width: f64,
    cell_height: f64,
    origin_cell: i32,
    dest_cell: i32,
    neighbours: i32,
    method: &str,
) -> List {
    if n_rows <= 0 || n_cols <= 0 {
        extendr_api::throw_r_error(format!(
            "n_rows and n_cols must be positive, got {} x {}",
            n_rows, n_cols
        ));
    }
    if !cell_width.is_finite() || cell_width <= 0.0 {
        extendr_api::throw_r_error(format!(
            "cell_width must be finite and positive, got {}",
            cell_width
        ));
    }
    if !cell_height.is_finite() || cell_height <= 0.0 {
        extendr_api::throw_r_error(format!(
            "cell_height must be finite and positive, got {}",
            cell_height
        ));
    }
    let nr = n_rows as usize;
    let nc = n_cols as usize;
    let n_cells = nr * nc;
    if values.len() != n_cells {
        extendr_api::throw_r_error(format!(
            "values length ({}) does not match n_rows * n_cols ({})",
            values.len(),
            n_cells
        ));
    }
    if neighbours != 4 && neighbours != 8 && neighbours != 16 {
        extendr_api::throw_r_error(format!(
            "neighbours must be 4, 8, or 16, got {}",
            neighbours
        ));
    }
    if method != "dijkstra" && method != "bidirectional" && method != "astar" {
        extendr_api::throw_r_error(format!(
            "method must be 'dijkstra', 'bidirectional', or 'astar', got '{}'",
            method
        ));
    }
    let oc = origin_cell as usize;
    let dc = dest_cell as usize;
    if origin_cell < 0 || oc >= n_cells {
        extendr_api::throw_r_error(format!(
            "origin_cell {} is out of range [0, {})",
            origin_cell, n_cells
        ));
    }
    if dest_cell < 0 || dc >= n_cells {
        extendr_api::throw_r_error(format!(
            "dest_cell {} is out of range [0, {})",
            dest_cell, n_cells
        ));
    }
    // Validate all friction values: must be NaN (NA/impassable) or finite and positive
    for (i, &v) in values.iter().enumerate() {
        if v.is_nan() {
            continue; // NA → impassable, allowed
        }
        if !v.is_finite() || v <= 0.0 {
            extendr_api::throw_r_error(format!(
                "All non-NA values must be finite and positive (> 0). Cell {} has value {}",
                i, v
            ));
        }
    }

    if values[oc].is_nan() {
        extendr_api::throw_r_error("origin cell falls on an impassable (NA) cell");
    }
    if values[dc].is_nan() {
        extendr_api::throw_r_error("destination cell falls on an impassable (NA) cell");
    }

    let result = route::corridor::solve(
        &values, nr, nc, cell_width, cell_height, oc, dc, neighbours, method,
    );

    if result.path_cells.is_empty() {
        extendr_api::throw_r_error(
            "No path exists between origin and destination. They may be separated by impassable (NA) cells."
        );
    }

    let path_cells_i32: Vec<i32> = result.path_cells.iter().map(|&c| c as i32).collect();
    list!(
        path_cells = path_cells_i32,
        total_cost = result.total_cost,
        solve_time_ms = result.solve_time_ms,
        graph_build_time_ms = result.graph_build_time_ms,
        n_edges = result.n_edges as i32
    )
}

/// Build a cached corridor graph from a raster grid
/// @noRd
#[extendr]
fn rust_corridor_build_graph(
    values: Vec<f64>,
    n_rows: i32,
    n_cols: i32,
    cell_width: f64,
    cell_height: f64,
    neighbours: i32,
) -> ExternalPtr<route::corridor::CorridorGraph> {
    if n_rows <= 0 || n_cols <= 0 {
        extendr_api::throw_r_error(format!(
            "n_rows and n_cols must be positive, got {} x {}", n_rows, n_cols
        ));
    }
    if !cell_width.is_finite() || cell_width <= 0.0 {
        extendr_api::throw_r_error(format!(
            "cell_width must be finite and positive, got {}", cell_width
        ));
    }
    if !cell_height.is_finite() || cell_height <= 0.0 {
        extendr_api::throw_r_error(format!(
            "cell_height must be finite and positive, got {}", cell_height
        ));
    }
    let nr = n_rows as usize;
    let nc = n_cols as usize;
    let n_cells = nr * nc;
    if values.len() != n_cells {
        extendr_api::throw_r_error(format!(
            "values length ({}) does not match n_rows * n_cols ({})",
            values.len(), n_cells
        ));
    }
    if neighbours != 4 && neighbours != 8 && neighbours != 16 {
        extendr_api::throw_r_error(format!(
            "neighbours must be 4, 8, or 16, got {}", neighbours
        ));
    }
    for (i, &v) in values.iter().enumerate() {
        if v.is_nan() { continue; }
        if !v.is_finite() || v <= 0.0 {
            extendr_api::throw_r_error(format!(
                "All non-NA values must be finite and positive (> 0). Cell {} has value {}", i, v
            ));
        }
    }

    let graph = route::corridor::build_corridor_graph(
        &values, nr, nc, cell_width, cell_height, neighbours,
    );
    ExternalPtr::new(graph)
}

/// Solve on a cached corridor graph
/// @noRd
#[extendr]
fn rust_corridor_solve_cached(
    graph: ExternalPtr<route::corridor::CorridorGraph>,
    origin_cell: i32,
    dest_cell: i32,
    method: &str,
) -> List {
    let n_cells = graph.n_rows * graph.n_cols;
    let oc = origin_cell as usize;
    let dc = dest_cell as usize;
    if origin_cell < 0 || oc >= n_cells {
        extendr_api::throw_r_error(format!(
            "origin_cell {} is out of range [0, {})", origin_cell, n_cells
        ));
    }
    if dest_cell < 0 || dc >= n_cells {
        extendr_api::throw_r_error(format!(
            "dest_cell {} is out of range [0, {})", dest_cell, n_cells
        ));
    }
    if method != "dijkstra" && method != "bidirectional" && method != "astar" {
        extendr_api::throw_r_error(format!(
            "method must be 'dijkstra', 'bidirectional', or 'astar', got '{}'", method
        ));
    }

    let result = route::corridor::solve_on_graph(&graph, oc, dc, method);

    if result.path_cells.is_empty() {
        extendr_api::throw_r_error(
            "No path exists between origin and destination. They may be separated by impassable (NA) cells."
        );
    }

    let path_cells_i32: Vec<i32> = result.path_cells.iter().map(|&c| c as i32).collect();
    list!(
        path_cells = path_cells_i32,
        total_cost = result.total_cost,
        solve_time_ms = result.solve_time_ms,
        graph_build_time_ms = result.graph_build_time_ms,
        n_edges = result.n_edges as i32
    )
}

/// Get metadata from a cached corridor graph
/// @noRd
#[extendr]
fn rust_corridor_graph_info(
    graph: ExternalPtr<route::corridor::CorridorGraph>,
) -> List {
    list!(
        n_rows = graph.n_rows as i32,
        n_cols = graph.n_cols as i32,
        cell_width = graph.cell_w,
        cell_height = graph.cell_h,
        neighbours = graph.neighbours,
        n_edges = graph.n_edges() as i32,
        build_time_ms = graph.build_time_ms,
        memory_bytes = graph.memory_bytes_estimate() as f64
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
    fn rust_azp;
    fn rust_spenc;
    fn rust_ward_constrained;
    fn rust_max_p;
    fn rust_frlm_greedy;
    fn rust_huff;
    fn rust_tsp;
    fn rust_vrp;
    fn rust_corridor;
    fn rust_corridor_build_graph;
    fn rust_corridor_solve_cached;
    fn rust_corridor_graph_info;
}
