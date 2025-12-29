//! AZP (Automatic Zoning Procedure)
//!
//! Implements the AZP algorithm with basic, tabu, and simulated annealing variants.
//! Based on: Openshaw (1977), Openshaw & Rao (1995)
//!
//! The algorithm performs local search to minimize within-region heterogeneity
//! while maintaining spatial contiguity.

use extendr_api::prelude::*;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha8Rng;
use std::collections::{HashSet, VecDeque};

/// Adjacency list representation
struct AdjList {
    neighbors: Vec<Vec<usize>>,
}

impl AdjList {
    fn from_indices(adj_i: &[i32], adj_j: &[i32], n: usize) -> Self {
        let mut neighbors = vec![Vec::new(); n];
        for (&i, &j) in adj_i.iter().zip(adj_j.iter()) {
            let i = i as usize;
            let j = j as usize;
            if !neighbors[i].contains(&j) {
                neighbors[i].push(j);
            }
            if !neighbors[j].contains(&i) {
                neighbors[j].push(i);
            }
        }
        Self { neighbors }
    }

    #[inline]
    fn get(&self, node: usize) -> &[usize] {
        &self.neighbors[node]
    }
}

/// Region state for AZP
#[derive(Clone)]
struct RegionState {
    labels: Vec<i32>,
    region_members: Vec<HashSet<usize>>,
    n_regions: usize,
}

impl RegionState {
    fn from_initial(labels: &[i32], n_regions: usize) -> Self {
        let mut region_members = vec![HashSet::new(); n_regions];
        for (i, &label) in labels.iter().enumerate() {
            if label >= 0 && (label as usize) < n_regions {
                region_members[label as usize].insert(i);
            }
        }
        Self {
            labels: labels.to_vec(),
            region_members,
            n_regions,
        }
    }

    /// Get border areas for a region (areas with neighbors in other regions)
    #[allow(dead_code)]
    fn get_border_areas(&self, region_id: usize, adj: &AdjList) -> Vec<usize> {
        self.region_members[region_id]
            .iter()
            .filter(|&&area| {
                adj.get(area).iter().any(|&n| {
                    self.labels[n] >= 0 && self.labels[n] as usize != region_id
                })
            })
            .copied()
            .collect()
    }

    /// Check if removing an area would disconnect its region
    fn is_articulation_point(&self, area: usize, adj: &AdjList) -> bool {
        let region_id = self.labels[area] as usize;
        let members = &self.region_members[region_id];

        if members.len() <= 2 {
            return members.len() == 2;
        }

        // Find a starting point
        let start = members.iter().find(|&&m| m != area).copied();
        if start.is_none() {
            return true;
        }
        let start = start.unwrap();

        // DFS to count reachable
        let mut visited = HashSet::new();
        let mut stack = vec![start];
        visited.insert(start);

        while let Some(current) = stack.pop() {
            for &neighbor in adj.get(current) {
                if neighbor != area
                    && members.contains(&neighbor)
                    && !visited.contains(&neighbor)
                {
                    visited.insert(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        visited.len() < members.len() - 1
    }

    /// Execute a move
    fn execute_move(&mut self, area: usize, to_region: usize) {
        let from_region = self.labels[area] as usize;
        self.region_members[from_region].remove(&area);
        self.region_members[to_region].insert(area);
        self.labels[area] = to_region as i32;
    }
}

/// Solve AZP regionalization
pub fn solve(
    attrs: RMatrix<f64>,
    n_regions: usize,
    adj_i: &[i32],
    adj_j: &[i32],
    method: &str,
    max_iterations: usize,
    tabu_length: usize,
    cooling_rate: f64,
    initial_temperature: f64,
    seed: Option<u64>,
) -> List {
    let n = attrs.nrows();
    let n_attrs = attrs.ncols();

    if n == 0 || n_regions == 0 || n_regions > n {
        return list!(
            labels = Vec::<i32>::new(),
            n_regions = 0i32,
            objective = 0.0f64
        );
    }

    // Convert R matrix to Vec<Vec<f64>>
    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    // Build adjacency list
    let adj = AdjList::from_indices(adj_i, adj_j, n);

    let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(42));

    // Generate initial solution
    let initial_labels = generate_initial_solution(n, n_regions, &adj, &mut rng);
    let mut state = RegionState::from_initial(&initial_labels, n_regions);

    // Run appropriate method
    match method {
        "tabu" => {
            state = azp_tabu(&attr_vecs, state, &adj, max_iterations, tabu_length, &mut rng);
        }
        "sa" => {
            state = azp_sa(&attr_vecs, state, &adj, max_iterations, cooling_rate, initial_temperature, &mut rng);
        }
        _ => {
            // Basic AZP (greedy local search)
            state = azp_basic(&attr_vecs, state, &adj, max_iterations, &mut rng);
        }
    }

    // Compute final objective
    let objective = compute_objective(&attr_vecs, &state);

    // Convert to 1-based labels
    let labels: Vec<i32> = state.labels.iter().map(|&l| l + 1).collect();

    list!(
        labels = labels,
        n_regions = state.n_regions as i32,
        objective = objective
    )
}

/// Generate initial solution using tree-based partitioning
fn generate_initial_solution(
    n: usize,
    n_regions: usize,
    adj: &AdjList,
    rng: &mut ChaCha8Rng,
) -> Vec<i32> {
    let mut labels = vec![-1i32; n];

    // Start with random seeds for each region
    let mut available: Vec<usize> = (0..n).collect();
    available.shuffle(rng);

    let mut seeds: Vec<usize> = Vec::new();
    let mut assigned = HashSet::new();

    // Pick n_regions seeds that are spread out
    for &candidate in &available {
        if seeds.len() >= n_regions {
            break;
        }
        // Check if far enough from existing seeds
        let far_enough = seeds.iter().all(|&s| {
            !adj.get(s).contains(&candidate)
        });
        if far_enough || seeds.is_empty() {
            seeds.push(candidate);
            labels[candidate] = seeds.len() as i32 - 1;
            assigned.insert(candidate);
        }
    }

    // If we couldn't find enough spread seeds, just pick randomly
    while seeds.len() < n_regions {
        for &candidate in &available {
            if !assigned.contains(&candidate) {
                seeds.push(candidate);
                labels[candidate] = seeds.len() as i32 - 1;
                assigned.insert(candidate);
                break;
            }
        }
    }

    // Grow regions using BFS from seeds
    let mut frontiers: Vec<VecDeque<usize>> = seeds.iter()
        .map(|&s| {
            let mut q = VecDeque::new();
            for &neighbor in adj.get(s) {
                if !assigned.contains(&neighbor) {
                    q.push_back(neighbor);
                }
            }
            q
        })
        .collect();

    // Round-robin assignment
    let mut round = 0;
    while assigned.len() < n {
        let region_id = round % n_regions;
        round += 1;

        // Try to get next area for this region
        while let Some(area) = frontiers[region_id].pop_front() {
            if assigned.contains(&area) {
                continue;
            }

            labels[area] = region_id as i32;
            assigned.insert(area);

            // Add neighbors to frontier
            for &neighbor in adj.get(area) {
                if !assigned.contains(&neighbor) {
                    frontiers[region_id].push_back(neighbor);
                }
            }
            break;
        }

        // If this region's frontier is empty but others have areas, steal from neighbors
        if frontiers[region_id].is_empty() && assigned.len() < n {
            // Find unassigned areas adjacent to this region
            for &area in &available {
                if !assigned.contains(&area) {
                    let adjacent_to_region = adj.get(area).iter().any(|&neighbor| {
                        labels[neighbor] == region_id as i32
                    });
                    if adjacent_to_region {
                        frontiers[region_id].push_back(area);
                    }
                }
            }
        }
    }

    labels
}

/// Basic AZP: greedy local search
fn azp_basic(
    attrs: &[Vec<f64>],
    mut state: RegionState,
    adj: &AdjList,
    max_iterations: usize,
    rng: &mut ChaCha8Rng,
) -> RegionState {
    let mut no_improve_count = 0;
    let max_no_improve = state.labels.len().max(50);

    for _iter in 0..max_iterations {
        if no_improve_count >= max_no_improve {
            break;
        }

        let mut improved = false;

        // Iterate through regions in random order
        let mut region_order: Vec<usize> = (0..state.n_regions).collect();
        region_order.shuffle(rng);

        for &recipient_region in &region_order {
            // Get candidate areas (border areas of neighboring regions)
            let mut candidates: Vec<(usize, usize)> = Vec::new(); // (area, donor_region)

            for &area in &state.region_members[recipient_region] {
                for &neighbor in adj.get(area) {
                    let donor = state.labels[neighbor];
                    if donor >= 0 && donor as usize != recipient_region {
                        candidates.push((neighbor, donor as usize));
                    }
                }
            }

            candidates.sort_unstable();
            candidates.dedup();
            candidates.shuffle(rng);

            for (area, donor_region) in candidates {
                // Check constraints
                if state.region_members[donor_region].len() <= 1 {
                    continue;
                }
                if state.is_articulation_point(area, adj) {
                    continue;
                }

                // Compute delta
                let delta = compute_move_delta(attrs, &state, area, donor_region, recipient_region);

                if delta < -1e-10 {
                    // Accept improving move
                    state.execute_move(area, recipient_region);
                    improved = true;
                    break;
                }
            }

            if improved {
                break;
            }
        }

        if improved {
            no_improve_count = 0;
        } else {
            no_improve_count += 1;
        }
    }

    state
}

/// AZP with Tabu search
fn azp_tabu(
    attrs: &[Vec<f64>],
    mut state: RegionState,
    adj: &AdjList,
    max_iterations: usize,
    tabu_length: usize,
    _rng: &mut ChaCha8Rng,
) -> RegionState {
    let mut best_state = state.clone();
    let mut best_objective = compute_objective(attrs, &state);
    let mut current_objective = best_objective;

    // Tabu list: stores (area, from_region, to_region)
    let mut tabu_list: VecDeque<(usize, usize, usize)> = VecDeque::new();

    let mut no_improve_count = 0;
    let max_no_improve = state.labels.len().max(100);

    for _iter in 0..max_iterations {
        if no_improve_count >= max_no_improve {
            break;
        }

        // Find all valid moves
        let mut moves: Vec<(usize, usize, usize, f64)> = Vec::new(); // (area, from, to, delta)

        for recipient in 0..state.n_regions {
            for &area in &state.region_members[recipient] {
                for &neighbor in adj.get(area) {
                    let donor = state.labels[neighbor];
                    if donor >= 0 && donor as usize != recipient {
                        let donor = donor as usize;

                        if state.region_members[donor].len() <= 1 {
                            continue;
                        }
                        if state.is_articulation_point(neighbor, adj) {
                            continue;
                        }

                        let delta = compute_move_delta(attrs, &state, neighbor, donor, recipient);
                        moves.push((neighbor, donor, recipient, delta));
                    }
                }
            }
        }

        if moves.is_empty() {
            break;
        }

        // Sort by delta (best first)
        moves.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap());

        // Find best non-tabu move (or aspiration)
        let mut selected_move: Option<(usize, usize, usize, f64)> = None;

        for &(area, from, to, delta) in &moves {
            let reverse = (area, to, from);
            let is_tabu = tabu_list.iter().any(|&m| m == reverse);

            if !is_tabu {
                selected_move = Some((area, from, to, delta));
                break;
            } else if current_objective + delta < best_objective {
                // Aspiration: accept tabu move if it's best ever
                selected_move = Some((area, from, to, delta));
                break;
            }
        }

        if let Some((area, from, to, delta)) = selected_move {
            state.execute_move(area, to);
            current_objective += delta;

            // Update tabu list
            tabu_list.push_back((area, from, to));
            if tabu_list.len() > tabu_length {
                tabu_list.pop_front();
            }

            // Update best
            if current_objective < best_objective {
                best_objective = current_objective;
                best_state = state.clone();
                no_improve_count = 0;
            } else {
                no_improve_count += 1;
            }
        } else {
            no_improve_count += 1;
        }
    }

    best_state
}

/// AZP with Simulated Annealing
fn azp_sa(
    attrs: &[Vec<f64>],
    mut state: RegionState,
    adj: &AdjList,
    max_iterations: usize,
    cooling_rate: f64,
    initial_temperature: f64,
    rng: &mut ChaCha8Rng,
) -> RegionState {
    let mut best_state = state.clone();
    let mut best_objective = compute_objective(attrs, &state);
    let mut current_objective = best_objective;

    // Set initial temperature
    let mut temperature = if initial_temperature > 0.0 {
        initial_temperature
    } else {
        best_objective / state.n_regions as f64
    };

    let mut no_improve_count = 0;
    let max_no_improve = state.labels.len().max(100);

    for _iter in 0..max_iterations {
        if no_improve_count >= max_no_improve || temperature < 1e-10 {
            break;
        }

        // Collect all valid moves
        let mut moves: Vec<(usize, usize, usize)> = Vec::new();

        for recipient in 0..state.n_regions {
            for &area in &state.region_members[recipient] {
                for &neighbor in adj.get(area) {
                    let donor = state.labels[neighbor];
                    if donor >= 0 && donor as usize != recipient {
                        let donor = donor as usize;

                        if state.region_members[donor].len() <= 1 {
                            continue;
                        }
                        if state.is_articulation_point(neighbor, adj) {
                            continue;
                        }

                        moves.push((neighbor, donor, recipient));
                    }
                }
            }
        }

        if moves.is_empty() {
            break;
        }

        // Pick random move
        let &(area, from, to) = moves.choose(rng).unwrap();
        let delta = compute_move_delta(attrs, &state, area, from, to);

        // Accept/reject
        let accept = if delta < 0.0 {
            true
        } else {
            let prob = (-delta / temperature).exp();
            rng.gen::<f64>() < prob
        };

        if accept {
            state.execute_move(area, to);
            current_objective += delta;

            if current_objective < best_objective {
                best_objective = current_objective;
                best_state = state.clone();
                no_improve_count = 0;
            } else {
                no_improve_count += 1;
            }
        } else {
            no_improve_count += 1;
        }

        // Cool down
        temperature *= cooling_rate;
    }

    best_state
}

/// Compute objective (total within-region SSD)
fn compute_objective(attrs: &[Vec<f64>], state: &RegionState) -> f64 {
    state.region_members
        .iter()
        .map(|members| compute_region_ssd(attrs, members))
        .sum()
}

/// Compute SSD for a single region
fn compute_region_ssd(attrs: &[Vec<f64>], members: &HashSet<usize>) -> f64 {
    if members.is_empty() {
        return 0.0;
    }

    let members_vec: Vec<usize> = members.iter().copied().collect();
    let n_attrs = attrs[0].len();
    let n = members_vec.len() as f64;

    // Compute centroid
    let mut centroid = vec![0.0; n_attrs];
    for &m in &members_vec {
        for (k, val) in attrs[m].iter().enumerate() {
            centroid[k] += val;
        }
    }
    for c in &mut centroid {
        *c /= n;
    }

    // Compute SSD
    let mut ssd = 0.0;
    for &m in &members_vec {
        for (k, val) in attrs[m].iter().enumerate() {
            let diff = val - centroid[k];
            ssd += diff * diff;
        }
    }

    ssd
}

/// Compute change in objective for a move
fn compute_move_delta(
    attrs: &[Vec<f64>],
    state: &RegionState,
    area: usize,
    from_region: usize,
    to_region: usize,
) -> f64 {
    let old_from_ssd = compute_region_ssd(attrs, &state.region_members[from_region]);
    let old_to_ssd = compute_region_ssd(attrs, &state.region_members[to_region]);

    // Simulate move
    let mut new_from = state.region_members[from_region].clone();
    new_from.remove(&area);
    let mut new_to = state.region_members[to_region].clone();
    new_to.insert(area);

    let new_from_ssd = compute_region_ssd(attrs, &new_from);
    let new_to_ssd = compute_region_ssd(attrs, &new_to);

    (new_from_ssd + new_to_ssd) - (old_from_ssd + old_to_ssd)
}
