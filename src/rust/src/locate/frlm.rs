//! FRLM (Flow Refueling Location Model)
//!
//! Implements a greedy heuristic for the Flow Refueling Location Model.
//! This model optimizes facility placement to maximize coverage of flows
//! (origin-destination paths) subject to vehicle range constraints.
//!
//! The FRLM models round-trip refueling with these key assumptions:
//! - Vehicles start at origin with HALF a tank (can travel R/2)
//! - At each open station, vehicles refuel to full (can travel R)
//! - A flow is "covered" if the round trip can be completed
//!
//! Based on:
//! - Kuby, M., & Lim, S. (2005). The flow-refueling location problem for
//!   alternative-fuel vehicles. Socio-Economic Planning Sciences, 39(2), 125-145.
//! - Capar, I., & Kuby, M. (2012). An efficient formulation of the flow refueling
//!   location model for alternative-fuel stations. IIE Transactions, 44(8), 622-636.

use extendr_api::prelude::*;
use std::collections::HashSet;

/// A flow (path) with its coverage requirements
#[derive(Debug)]
struct Flow {
    /// Flow ID (0-based)
    id: usize,
    /// Volume of flow
    volume: f64,
    /// Ordered list of candidate facility indices along the path
    candidates: Vec<usize>,
    /// Distances from origin to each candidate
    distances: Vec<f64>,
}

impl Flow {
    /// Check if a set of facilities can cover this flow given vehicle range
    ///
    /// Based on the Flow Refueling Location Model (FRLM) from:
    /// - Kuby & Lim (2005): Original FRLM formulation
    /// - Capar & Kuby (2012): Efficient FRLM formulation
    ///
    /// Key assumptions:
    /// 1. Vehicles start at origin with HALF a tank (can travel R/2)
    /// 2. At each open station, vehicles refuel to full (can travel R)
    /// 3. Flow is covered if the ROUND TRIP can be completed
    ///
    /// Coverage requirements:
    /// - First open station must be within R/2 from origin (half-tank start)
    /// - Each subsequent open station must be within R from previous
    /// - Last open station must be within R/2 of destination (round-trip back)
    fn is_covered_by(&self, facilities: &HashSet<usize>, vehicle_range: f64) -> bool {
        if self.candidates.is_empty() || self.distances.is_empty() {
            return false;
        }

        let half_range = vehicle_range / 2.0;

        // Path length is the distance from origin to destination
        // (destination should be the last candidate in the sorted list)
        let path_length = *self.distances.last().unwrap();

        // Special case: if round trip distance <= range, any single station can refuel
        // Start with half tank, refuel once, complete trip and return with fuel to spare
        if 2.0 * path_length <= vehicle_range {
            return self.candidates.iter().any(|c| facilities.contains(c));
        }

        // Collect distances to open stations on this path
        // (candidates are already sorted by distance from origin)
        let open_station_distances: Vec<f64> = self
            .candidates
            .iter()
            .zip(self.distances.iter())
            .filter(|(&c, _)| facilities.contains(&c))
            .map(|(_, &d)| d)
            .collect();

        if open_station_distances.is_empty() {
            return false;
        }

        // Check 1: First open station must be within R/2 from origin
        // (driver starts with half a tank per Kuby & Lim 2005)
        if open_station_distances[0] > half_range {
            return false;
        }

        // Check 2: Each subsequent open station must be within R of previous
        for i in 1..open_station_distances.len() {
            if open_station_distances[i] - open_station_distances[i - 1] > vehicle_range {
                return false;
            }
        }

        // Check 3: Last open station must be within R/2 of destination
        // This ensures the round trip can be completed:
        // last_station -> destination -> last_station requires 2 * (path_length - last_station_dist) <= R
        let last_station_dist = *open_station_distances.last().unwrap();
        if path_length - last_station_dist > half_range {
            return false;
        }

        true
    }
}

/// Solve FRLM using greedy heuristic
///
/// # Arguments
/// * `n_candidates` - Total number of candidate facility locations
/// * `path_candidates` - List of candidate indices for each path (flat, with separators)
/// * `path_offsets` - Start index in path_candidates for each path
/// * `path_distances` - Distances to each candidate (flat, parallel to path_candidates)
/// * `flow_volumes` - Volume of each flow/path
/// * `vehicle_range` - Maximum vehicle range
/// * `n_facilities` - Number of facilities to place
///
/// # Returns
/// List with selected facilities, coverage info, and objective
pub fn solve_greedy(
    n_candidates: usize,
    path_candidates: &[i32],
    path_offsets: &[i32],
    path_distances: &[f64],
    flow_volumes: &[f64],
    vehicle_range: f64,
    n_facilities: usize,
) -> List {
    let n_flows = flow_volumes.len();

    if n_flows == 0 || n_facilities == 0 || n_candidates == 0 {
        return list!(
            selected = Vec::<i32>::new(),
            n_selected = 0i32,
            covered_volume = 0.0f64,
            total_volume = 0.0f64,
            coverage_pct = 0.0f64
        );
    }

    // Parse flows from flat arrays
    let mut flows: Vec<Flow> = Vec::with_capacity(n_flows);

    for flow_id in 0..n_flows {
        let start_idx = path_offsets[flow_id] as usize;
        let end_idx = if flow_id + 1 < path_offsets.len() {
            path_offsets[flow_id + 1] as usize
        } else {
            path_candidates.len()
        };

        let candidates: Vec<usize> = path_candidates[start_idx..end_idx]
            .iter()
            .map(|&c| c as usize)
            .collect();

        let distances: Vec<f64> = path_distances[start_idx..end_idx].to_vec();

        flows.push(Flow {
            id: flow_id,
            volume: flow_volumes[flow_id],
            candidates,
            distances,
        });
    }

    let total_volume: f64 = flow_volumes.iter().sum();

    // Greedy facility selection
    let mut selected_facilities: HashSet<usize> = HashSet::new();
    let mut covered_flows: HashSet<usize> = HashSet::new();

    for _ in 0..n_facilities {
        let mut best_facility: Option<usize> = None;
        let mut best_marginal_coverage = 0.0;

        // Evaluate each candidate facility
        for candidate in 0..n_candidates {
            if selected_facilities.contains(&candidate) {
                continue;
            }

            // Temporarily add this facility
            let mut test_set = selected_facilities.clone();
            test_set.insert(candidate);

            // Calculate marginal coverage
            let mut marginal_coverage = 0.0;
            for flow in &flows {
                if covered_flows.contains(&flow.id) {
                    continue;
                }
                if flow.is_covered_by(&test_set, vehicle_range) {
                    marginal_coverage += flow.volume;
                }
            }

            if marginal_coverage > best_marginal_coverage {
                best_marginal_coverage = marginal_coverage;
                best_facility = Some(candidate);
            }
        }

        // Add best facility
        if let Some(facility) = best_facility {
            selected_facilities.insert(facility);

            // Update covered flows
            for flow in &flows {
                if !covered_flows.contains(&flow.id) {
                    if flow.is_covered_by(&selected_facilities, vehicle_range) {
                        covered_flows.insert(flow.id);
                    }
                }
            }
        } else {
            // No more improvement possible
            break;
        }
    }

    // Calculate final coverage
    let covered_volume: f64 = covered_flows
        .iter()
        .map(|&id| flows[id].volume)
        .sum();

    let coverage_pct = if total_volume > 0.0 {
        100.0 * covered_volume / total_volume
    } else {
        0.0
    };

    // Convert to 1-based indices for R
    let selected: Vec<i32> = selected_facilities
        .iter()
        .map(|&f| f as i32 + 1)
        .collect();

    list!(
        selected = selected,
        n_selected = selected_facilities.len() as i32,
        covered_volume = covered_volume,
        total_volume = total_volume,
        coverage_pct = coverage_pct
    )
}
