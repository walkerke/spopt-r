//! Max-P Regions
//!
//! Maximize the number of regions such that each region satisfies a minimum threshold constraint.
//! Based on: Duque, Anselin & Rey (2012) / Wei, Rey, and Knaap (2020)
//! Extended with compactness support based on: Feng, Rey, & Wei (2022)
//!
//! This implementation is optimized for speed:
//! - Parallel construction phase using rayon
//! - Efficient union-find for connectivity checking
//! - Articulation point detection for move eligibility
//! - Incremental threshold tracking
//! - Early termination when max_p stabilizes
//! - Optional compactness optimization using either NMI or centroid dispersion

use extendr_api::prelude::*;
use rand::prelude::*;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::local_search::{anneal, greedy, Adj, Floor, Partition};

/// Compactness metric to use for region evaluation
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CompactnessMetric {
    CentroidDispersion,
    NMI,
}

/// Centroid coordinates for compactness calculation
#[derive(Clone)]
struct Centroids {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl Centroids {
    fn new(x: Option<Vec<f64>>, y: Option<Vec<f64>>) -> Option<Self> {
        match (x, y) {
            (Some(x), Some(y)) if x.len() == y.len() => Some(Self { x, y }),
            _ => None,
        }
    }
}

/// Adjacency list representation for fast neighbor lookup
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

/// Region state for efficient tracking
#[derive(Clone)]
struct RegionState {
    /// Which region each area belongs to (-1 = unassigned/enclave)
    labels: Vec<i32>,
    /// Areas in each region
    region_members: Vec<Vec<usize>>,
    /// Sum of threshold variable for each region
    region_thresholds: Vec<f64>,
    /// Number of valid regions
    n_regions: usize,
}

impl RegionState {
    fn new(n: usize) -> Self {
        Self {
            labels: vec![-1; n],
            region_members: Vec::new(),
            region_thresholds: Vec::new(),
            n_regions: 0,
        }
    }

    /// Add a new region with the given seed area
    fn add_region(&mut self, seed: usize, threshold_var: &[f64]) -> usize {
        let region_id = self.n_regions;
        self.labels[seed] = region_id as i32;
        self.region_members.push(vec![seed]);
        self.region_thresholds.push(threshold_var[seed]);
        self.n_regions += 1;
        region_id
    }

    /// Add area to existing region
    #[inline]
    fn add_to_region(&mut self, area: usize, region_id: usize, threshold_var: &[f64]) {
        self.labels[area] = region_id as i32;
        self.region_members[region_id].push(area);
        self.region_thresholds[region_id] += threshold_var[area];
    }

    /// Check if region meets threshold
    #[inline]
    fn region_meets_threshold(&self, region_id: usize, threshold: f64) -> bool {
        self.region_thresholds[region_id] >= threshold
    }

    /// Get all enclaves (unassigned areas)
    fn get_enclaves(&self) -> Vec<usize> {
        self.labels
            .iter()
            .enumerate()
            .filter(|(_, &l)| l < 0)
            .map(|(i, _)| i)
            .collect()
    }
}

/// Solve Max-P regionalization
///
/// Returns labels (1-based for R), n_regions, objective value, and compactness metrics
pub fn solve(
    attrs: RMatrix<f64>,
    threshold_var: &[f64],
    threshold: f64,
    adj_i: &[i32],
    adj_j: &[i32],
    n_iterations: usize,
    n_sa_iterations: usize,
    cooling_rate: f64,
    tabu_length: usize,
    seed: Option<u64>,
    centroids_x: Option<Vec<f64>>,
    centroids_y: Option<Vec<f64>>,
    compact: bool,
    compact_weight: f64,
    homogeneous: bool,
    metric: CompactnessMetric,
    areas: Option<Vec<f64>>,
    moments: Option<Vec<f64>>
) -> List {
    let n = attrs.nrows();
    let n_attrs = attrs.ncols();

    if n == 0 {
        return list!(
            labels = Vec::<i32>::new(),
            n_regions = 0i32,
            objective = 0.0f64
        );
    }

    // Check if threshold is achievable
    let total_threshold: f64 = threshold_var.iter().sum();
    if total_threshold < threshold {
        // Impossible to create even one region
        return list!(
            labels = vec![1i32; n],
            n_regions = 1i32,
            objective = f64::MAX
        );
    }

    // Convert R matrix to Vec<Vec<f64>> (row-major)
    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    // Build adjacency list
    let adj = AdjList::from_indices(adj_i, adj_j, n);

    // Parse centroids for compactness (if provided)
    let centroids = Centroids::new(centroids_x, centroids_y);
    let use_compactness = compact && centroids.is_some();
    let centroids_ref = centroids.as_ref();
    let areas_ref = areas.as_deref();
    let moments_ref = moments.as_deref();

    // Determine number of parallel workers
    let n_workers = rayon::current_num_threads().min(n_iterations).max(1);
    let iterations_per_worker = (n_iterations + n_workers - 1) / n_workers;

    // Track best result across all workers
    let best_n_regions = AtomicUsize::new(0);

    // Clone centroids for parallel access
    let centroids_clone = centroids.clone();
    let attrs_ref = &attr_vecs;

    // Run parallel construction phase
    let base_seed = seed.unwrap_or(42);
    let results: Vec<_> = (0..n_workers)
        .into_par_iter()
        .map(|worker_id| {
            let mut rng = ChaCha8Rng::seed_from_u64(base_seed.wrapping_add(worker_id as u64));
            let mut worker_best: Option<RegionState> = None;
            let mut worker_best_p = 0usize;
            let mut worker_best_obj = f64::INFINITY;

            for _iter in 0..iterations_per_worker {
                // Construction is cheap (one BFS per region), so every requested
                // iteration is run: the maximum p is often found late, and an
                // early stop after a fixed number of non-improving tries was
                // leaving regions on the table.

                // Check global best - if another worker found much better, we can stop too
                let global_best = best_n_regions.load(Ordering::Relaxed);
                if worker_best_p > 0 && global_best > worker_best_p + 2 {
                    break;
                }

                // Construction phase
                let state = construct_solution(
                    n,
                    attrs_ref,
                    threshold_var,
                    threshold,
                    &adj,
                    &mut rng,
                    use_compactness,
                    centroids_clone.as_ref(),
                    metric,
                    areas_ref,
                    moments_ref,
                );

                if state.n_regions > worker_best_p {
                    worker_best_p = state.n_regions;
                    worker_best_obj = compute_objective(attrs_ref, &state);
                    worker_best = Some(state);

                    // Update global best
                    let _ = best_n_regions.fetch_max(worker_best_p, Ordering::Relaxed);
                } else if state.n_regions == worker_best_p {
                    // Same p: keep the more homogeneous construction
                    let obj = compute_objective(attrs_ref, &state);
                    if obj < worker_best_obj {
                        worker_best_obj = obj;
                        worker_best = Some(state);
                    }
                }
            }

            worker_best
        })
        .collect();

    // Find best construction result
    let mut best_state = results
        .into_iter()
        .flatten()
        .map(|s| (s.n_regions, compute_objective(&attr_vecs, &s), s))
        .fold(None::<(usize, f64, RegionState)>, |acc, cur| match acc {
            None => Some(cur),
            Some(b) => {
                if cur.0 > b.0 || (cur.0 == b.0 && cur.1 < b.1) {
                    Some(cur)
                } else {
                    Some(b)
                }
            }
        })
        .map(|(_, _, s)| s);

    if best_state.is_none() {
        // Fallback: single region containing everything
        let mut state = RegionState::new(n);
        state.add_region(0, threshold_var);
        for i in 1..n {
            state.add_to_region(i, 0, threshold_var);
        }
        best_state = Some(state);
    }

    let mut state = best_state.unwrap();

    // Improvement phase. With the default objective (homogeneous regions, no
    // compactness term) this is a floor-constrained AZP: a greedy pass to a
    // local optimum, simulated annealing over the border, and a final greedy
    // polish. The compactness and heterogeneous objectives keep the original
    // annealing routine.
    if n_sa_iterations > 0 && state.n_regions > 1 && homogeneous && !use_compactness {
        let floor = Floor { var: threshold_var, min: threshold };
        let ls_adj = Adj::from_indices(adj_i, adj_j, n);
        let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(42).wrapping_add(999));
        let mut part = Partition::new(&state.labels, state.n_regions, &attr_vecs, Some(floor));
        greedy(&mut part, &attr_vecs, &ls_adj, Some(floor), 1000, &mut rng);
        part = anneal(
            part,
            &attr_vecs,
            &ls_adj,
            Some(floor),
            n_sa_iterations,
            cooling_rate,
            0.0,
            tabu_length,
            &mut rng,
        );
        greedy(&mut part, &attr_vecs, &ls_adj, Some(floor), 1000, &mut rng);
        state = RegionState {
            labels: part.labels.clone(),
            region_members: part.members.clone(),
            region_thresholds: part.floor_sums.clone(),
            n_regions: part.k(),
        };
    } else if n_sa_iterations > 0 && state.n_regions > 1 {
        state = simulated_annealing(
            state,
            &attr_vecs,
            threshold_var,
            threshold,
            &adj,
            n_sa_iterations,
            cooling_rate,
            tabu_length,
            seed.map(|s| s.wrapping_add(999)),
            use_compactness,
            centroids_ref,
            compact_weight,
            homogeneous,
            metric,
            areas_ref,
            moments_ref,
        );
    }

    // Compute final objective (total within-region dissimilarity)
    let objective = compute_objective(&attr_vecs, &state);

    // Compute compactness metrics
    let region_compactness: Vec<f64> = if let Some(ref cents) = centroids {
        compute_all_region_compactness(
            &state, 
            cents,
            metric,
            areas_ref,
            moments_ref
        )
    } else {
        vec![0.0; state.n_regions]
    };
    let mean_compactness: f64 = if !region_compactness.is_empty() {
        region_compactness.iter().sum::<f64>() / region_compactness.len() as f64
    } else {
        0.0
    };

    // Convert to 1-based labels for R
    let labels: Vec<i32> = state.labels.iter().map(|&l| l + 1).collect();

    list!(
        labels = labels,
        n_regions = state.n_regions as i32,
        objective = objective,
        mean_compactness = mean_compactness,
        region_compactness = region_compactness
    )
}

/// Construct a single solution using greedy region growing
fn construct_solution(
    n: usize,
    attrs: &[Vec<f64>],
    threshold_var: &[f64],
    threshold: f64,
    adj: &AdjList,
    rng: &mut ChaCha8Rng,
    use_compactness: bool,
    centroids: Option<&Centroids>,
    metric: CompactnessMetric,
    areas: Option<&[f64]>,
    moments: Option<&[f64]>,
) -> RegionState {
    let mut state = RegionState::new(n);

    // Random order for seed selection
    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(rng);

    // Phase 1: Create regions that meet threshold
    for &seed in &order {
        if state.labels[seed] >= 0 {
            continue; // Already assigned
        }

        // Start a new region from this seed
        let region_id = state.add_region(seed, threshold_var);

        // Grow region until threshold is met
        grow_region(&mut state, region_id, threshold_var, threshold, adj, rng,
                    use_compactness, centroids, metric, areas, moments);

        // If region doesn't meet threshold, dissolve it back to enclaves
        if !state.region_meets_threshold(region_id, threshold) {
            // Mark all members of this failed region as enclaves
            for &member in &state.region_members[region_id].clone() {
                state.labels[member] = -1;
            }
            state.region_members[region_id].clear();
            state.region_thresholds[region_id] = 0.0;
        }
    }

    // Phase 2: Assign all enclaves to valid regions
    assign_enclaves_to_valid_regions(&mut state, attrs, threshold_var, threshold, adj, rng);

    // Compact region IDs (remove empty regions)
    compact_regions(&mut state);

    state
}

/// Grow a region by greedily adding neighbors until threshold is met
fn grow_region(
    state: &mut RegionState,
    region_id: usize,
    threshold_var: &[f64],
    threshold: f64,
    adj: &AdjList,
    rng: &mut ChaCha8Rng,
    use_compactness: bool,
    centroids: Option<&Centroids>,
    metric: CompactnessMetric,
    areas: Option<&[f64]>,
    moments: Option<&[f64]>,
) {
    // Use BFS-like expansion with randomization
    let mut frontier: Vec<usize> = Vec::new();

    // Initialize frontier with neighbors of seed
    let seed = state.region_members[region_id][0];
    for &neighbor in adj.get(seed) {
        if state.labels[neighbor] < 0 {
            frontier.push(neighbor);
        }
    }

    while !state.region_meets_threshold(region_id, threshold) && !frontier.is_empty() {
        let candidate = if use_compactness && centroids.is_some() {
            // Select from top N candidates that optimize compactness
            // (per Feng et al. 2022: pick from top 3 most compact choices)
            let cents = centroids.unwrap();
            let top_n = 3.min(frontier.len());

            // Score each frontier candidate by how compact the region would be after adding it
            let mut scored: Vec<(usize, f64)> = frontier.iter()
                .filter(|&&c| state.labels[c] < 0)
                .map(|&c| {
                    let score = compute_region_compactness_if_added(
                        state, 
                        region_id, 
                        c, 
                        cents,
                        metric,
                        areas,
                        moments
                    );
                    (c, score)
                })
                .collect();

            // Sort by compactness (higher = better)
            scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

            if scored.is_empty() {
                break;
            }

            // Pick randomly from top N
            let pick_idx = rng.gen_range(0..top_n.min(scored.len()));
            let chosen = scored[pick_idx].0;

            // Remove from frontier
            if let Some(pos) = frontier.iter().position(|&x| x == chosen) {
                frontier.swap_remove(pos);
            }
            chosen
        } else {
            // If one frontier area is enough to reach the threshold, take the
            // one that overshoots it least, so the remaining areas can form
            // more regions (this is what max-p maximizes). Otherwise pick a
            // random frontier area, as in Duque, Anselin & Rey (2012).
            let need = threshold - state.region_thresholds[region_id];
            let mut finisher: Option<(usize, f64)> = None;
            for (pos, &c) in frontier.iter().enumerate() {
                if state.labels[c] < 0 && threshold_var[c] >= need {
                    let over = threshold_var[c] - need;
                    if finisher.map_or(true, |(_, best)| over < best) {
                        finisher = Some((pos, over));
                    }
                }
            }
            let idx = match finisher {
                Some((pos, _)) => pos,
                None => rng.gen_range(0..frontier.len()),
            };
            frontier.swap_remove(idx)
        };

        if state.labels[candidate] >= 0 {
            continue; // Already assigned
        }

        // Add to region
        state.add_to_region(candidate, region_id, threshold_var);

        // Expand frontier
        for &neighbor in adj.get(candidate) {
            if state.labels[neighbor] < 0 && !frontier.contains(&neighbor) {
                frontier.push(neighbor);
            }
        }
    }
}

/// Increase in a region's SSD if `area` joins it: n/(n+1) * |x - mean|^2.
fn ssd_increase(attrs: &[Vec<f64>], members: &[usize], area: usize) -> f64 {
    if members.is_empty() {
        return 0.0;
    }
    let n = members.len() as f64;
    let p = attrs[area].len();
    let mut d = 0.0;
    for c in 0..p {
        let mean = members.iter().map(|&m| attrs[m][c]).sum::<f64>() / n;
        let diff = attrs[area][c] - mean;
        d += diff * diff;
    }
    n / (n + 1.0) * d
}

/// Assign enclaves to valid regions (regions that meet threshold).
///
/// An enclave touching several valid regions joins the one whose
/// heterogeneity it raises least (Duque, Anselin & Rey 2012); the order in
/// which enclaves are processed stays random.
fn assign_enclaves_to_valid_regions(
    state: &mut RegionState,
    attrs: &[Vec<f64>],
    threshold_var: &[f64],
    threshold: f64,
    adj: &AdjList,
    rng: &mut ChaCha8Rng,
) {
    // Get valid regions (those meeting threshold)
    let valid_regions: Vec<usize> = (0..state.n_regions)
        .filter(|&r| state.region_meets_threshold(r, threshold))
        .collect();

    if valid_regions.is_empty() {
        // No valid regions - put everything in region 0
        // This is a fallback for impossible cases
        return;
    }

    let mut enclaves = state.get_enclaves();
    let mut max_attempts = enclaves.len() * 3 + 100;

    while !enclaves.is_empty() && max_attempts > 0 {
        max_attempts -= 1;
        enclaves.shuffle(rng);

        let mut assigned_any = false;
        let mut new_enclaves = Vec::new();

        for enclave in enclaves {
            // Find neighboring VALID regions
            let mut neighbor_regions: Vec<usize> = adj.get(enclave)
                .iter()
                .filter_map(|&n| {
                    let label = state.labels[n];
                    if label >= 0 {
                        let r = label as usize;
                        if valid_regions.contains(&r) {
                            return Some(r);
                        }
                    }
                    None
                })
                .collect();

            neighbor_regions.sort_unstable();
            neighbor_regions.dedup();

            if neighbor_regions.is_empty() {
                // No valid neighbor - check if we have any assigned neighbor
                let any_neighbor: Vec<usize> = adj.get(enclave)
                    .iter()
                    .filter_map(|&n| {
                        let label = state.labels[n];
                        if label >= 0 {
                            Some(label as usize)
                        } else {
                            None
                        }
                    })
                    .collect();

                if any_neighbor.is_empty() {
                    new_enclaves.push(enclave);
                } else {
                    // Assign to any neighboring region (it will be in a valid chain eventually)
                    let region_id = *any_neighbor.choose(rng).unwrap();
                    state.add_to_region(enclave, region_id, threshold_var);
                    assigned_any = true;
                }
                continue;
            }

            // Join the most similar valid neighbouring region
            let region_id = neighbor_regions
                .iter()
                .copied()
                .min_by(|&a, &b| {
                    ssd_increase(attrs, &state.region_members[a], enclave)
                        .partial_cmp(&ssd_increase(attrs, &state.region_members[b], enclave))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap();
            state.add_to_region(enclave, region_id, threshold_var);
            assigned_any = true;
        }

        enclaves = new_enclaves;

        if !assigned_any && !enclaves.is_empty() {
            // Stuck - force assignment to first valid region through path
            break;
        }
    }

    // Handle any remaining truly disconnected enclaves
    // Assign them to the closest valid region via BFS
    for enclave in enclaves {
        if state.labels[enclave] < 0 {
            if let Some(target_region) = find_closest_valid_region(state, enclave, &valid_regions, adj) {
                state.add_to_region(enclave, target_region, threshold_var);
            } else if !valid_regions.is_empty() {
                // Last resort: assign to first valid region
                state.add_to_region(enclave, valid_regions[0], threshold_var);
            }
        }
    }
}

/// Find the closest valid region via BFS through the graph
fn find_closest_valid_region(
    state: &RegionState,
    start: usize,
    valid_regions: &[usize],
    adj: &AdjList,
) -> Option<usize> {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();

    visited.insert(start);
    queue.push_back(start);

    while let Some(current) = queue.pop_front() {
        for &neighbor in adj.get(current) {
            if !visited.contains(&neighbor) {
                visited.insert(neighbor);

                let label = state.labels[neighbor];
                if label >= 0 && valid_regions.contains(&(label as usize)) {
                    return Some(label as usize);
                }

                queue.push_back(neighbor);
            }
        }
    }

    None
}

/// Compact region IDs by removing empty regions
fn compact_regions(state: &mut RegionState) {
    // Find non-empty regions
    let non_empty: Vec<usize> = (0..state.n_regions)
        .filter(|&r| !state.region_members[r].is_empty())
        .collect();

    if non_empty.len() == state.n_regions {
        return; // No compaction needed
    }

    // Create mapping from old to new IDs
    let mut old_to_new: HashMap<usize, usize> = HashMap::new();
    for (new_id, &old_id) in non_empty.iter().enumerate() {
        old_to_new.insert(old_id, new_id);
    }

    // Update labels
    for label in &mut state.labels {
        if *label >= 0 {
            *label = old_to_new[&(*label as usize)] as i32;
        }
    }

    // Rebuild region_members and region_thresholds
    let mut new_members = Vec::new();
    let mut new_thresholds = Vec::new();

    for &old_id in &non_empty {
        new_members.push(state.region_members[old_id].clone());
        new_thresholds.push(state.region_thresholds[old_id]);
    }

    state.region_members = new_members;
    state.region_thresholds = new_thresholds;
    state.n_regions = non_empty.len();
}

/// Simulated annealing to minimize within-region dissimilarity (and maximize compactness)
fn simulated_annealing(
    mut state: RegionState,
    attrs: &[Vec<f64>],
    threshold_var: &[f64],
    threshold: f64,
    adj: &AdjList,
    n_iterations: usize,
    cooling_rate: f64,
    tabu_length: usize,
    seed: Option<u64>,
    use_compactness: bool,
    centroids: Option<&Centroids>,
    compact_weight: f64,
    homogeneous: bool,
    metric: CompactnessMetric,
    areas: Option<&[f64]>,
    moments: Option<&[f64]>,
) -> RegionState {
    let n = state.labels.len();
    let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(12345));

    // Initial temperature based on typical objective values
    let mut temperature = compute_objective(attrs, &state) / (state.n_regions as f64).max(1.0);
    if temperature < 1e-10 {
        temperature = 1.0;
    }

    // Tabu list: (area, from_region, to_region)
    let mut tabu_list: VecDeque<(usize, usize, usize)> = VecDeque::new();

    // Compute initial combined objective
    let dissimilarity_obj = compute_objective(attrs, &state);
    let compactness_obj = if use_compactness && centroids.is_some() {
        let comps = compute_all_region_compactness(
            &state, 
            centroids.unwrap(),
            metric,
            areas,
            moments
        );
        comps.iter().sum::<f64>() / comps.len().max(1) as f64
    } else {
        0.0
    };

    // Combined objective: minimize dissimilarity (or maximize if homogeneous = false), maximize compactness
    // We negate compactness since we're minimizing
    // dissim_sign flips the objective when we want heterogeneous regions
    let dissim_sign = if homogeneous { 1.0 } else { -1.0 };
    let mut current_objective = (1.0 - compact_weight) * dissim_sign * dissimilarity_obj - compact_weight * compactness_obj * 100.0;
    let mut best_state = state.clone();
    let mut best_objective = current_objective;

    let mut no_improve_count = 0;
    let max_no_improve = n.max(100);

    for _iter in 0..n_iterations {
        if no_improve_count > max_no_improve {
            break;
        }

        // Find moveable areas (border areas that can move without breaking connectivity or threshold)
        let moveable = find_moveable_areas(&state, threshold_var, threshold, adj);

        if moveable.is_empty() {
            no_improve_count += 1;
            continue;
        }

        // Pick random moveable area
        let &(area, from_region, to_region) = moveable.choose(&mut rng).unwrap();

        // Check tabu
        let move_tuple = (area, from_region, to_region);
        let reverse_tuple = (area, to_region, from_region);

        // Compute objective change (dissimilarity component)
        let dissim_delta = compute_move_delta(attrs, &state, area, from_region, to_region);

        // Compute compactness change if enabled
        let compact_delta = if use_compactness && centroids.is_some() {
            let cents = centroids.unwrap();
            let old_from = compute_region_compactness(
                &state.region_members[from_region], 
                cents,
                metric,
                areas,
                moments
            );
            let old_to = compute_region_compactness(
                &state.region_members[to_region], 
                cents,
                metric,
                areas,
                moments
            );

            // Simulate move
            let new_from: Vec<usize> = state.region_members[from_region].iter()
                .copied().filter(|&m| m != area).collect();
            let mut new_to = state.region_members[to_region].clone();
            new_to.push(area);

            let new_from_comp = compute_region_compactness(
                &new_from, 
                cents,
                metric,
                areas,
                moments
            );
            let new_to_comp = compute_region_compactness(
                &new_to, 
                cents,
                metric,
                areas,
                moments
            );

            // Return change (positive = improved compactness)
            (new_from_comp + new_to_comp) - (old_from + old_to)
        } else {
            0.0
        };

        // Combined delta: minimize dissimilarity (or maximize if homogeneous = false),
        //                 maximize compactness
        let delta = (1.0 - compact_weight) * dissim_sign * dissim_delta - compact_weight * compact_delta * 100.0;

        if tabu_list.contains(&reverse_tuple) {
            // Skip tabu move unless it's aspiration (would be best ever)
            if current_objective + delta >= best_objective {
                no_improve_count += 1;
                continue;
            }
        }

        // Accept/reject
        let accept = if delta < 0.0 {
            true
        } else {
            let prob = (-delta / temperature).exp();
            rng.gen::<f64>() < prob
        };

        if accept {
            // Execute move
            execute_move(&mut state, area, from_region, to_region, threshold_var);
            current_objective += delta;

            // Update tabu list
            tabu_list.push_back(move_tuple);
            if tabu_list.len() > tabu_length {
                tabu_list.pop_front();
            }

            // Track best
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

/// Find areas that can move between regions without violating constraints
fn find_moveable_areas(
    state: &RegionState,
    threshold_var: &[f64],
    threshold: f64,
    adj: &AdjList,
) -> Vec<(usize, usize, usize)> {
    let mut moveable = Vec::new();
    let n = state.labels.len();

    for area in 0..n {
        let from_region = state.labels[area];
        if from_region < 0 {
            continue;
        }
        let from_region = from_region as usize;

        // Check if removal would violate threshold
        let new_threshold = state.region_thresholds[from_region] - threshold_var[area];
        if new_threshold < threshold {
            continue;
        }

        // Check if area is an articulation point (would disconnect region)
        if is_articulation_point(state, area, from_region, adj) {
            continue;
        }

        // Find neighboring regions this area could move to
        for &neighbor in adj.get(area) {
            let to_region = state.labels[neighbor];
            if to_region >= 0 && to_region as usize != from_region {
                moveable.push((area, from_region, to_region as usize));
            }
        }
    }

    moveable
}

/// Check if removing an area would disconnect its region
/// Uses DFS from a remaining member to check if all others are reachable
fn is_articulation_point(
    state: &RegionState,
    area: usize,
    region_id: usize,
    adj: &AdjList,
) -> bool {
    let members = &state.region_members[region_id];
    if members.len() <= 2 {
        return members.len() == 2; // Can't remove from size-2 region
    }

    // Find a starting point that's not the area being removed
    let start = members.iter().find(|&&m| m != area).copied();
    if start.is_none() {
        return true;
    }
    let start = start.unwrap();

    // DFS to count reachable members
    let mut visited = HashSet::new();
    let mut stack = vec![start];
    visited.insert(start);

    while let Some(current) = stack.pop() {
        for &neighbor in adj.get(current) {
            if neighbor != area
                && state.labels[neighbor] == region_id as i32
                && !visited.contains(&neighbor)
            {
                visited.insert(neighbor);
                stack.push(neighbor);
            }
        }
    }

    // If we can't reach all other members, it's an articulation point
    visited.len() < members.len() - 1
}

/// Compute change in objective if area moves from one region to another
fn compute_move_delta(
    attrs: &[Vec<f64>],
    state: &RegionState,
    area: usize,
    from_region: usize,
    to_region: usize,
) -> f64 {
    // Compute SSD contribution before and after
    let old_from_ssd = compute_region_ssd(attrs, &state.region_members[from_region]);
    let old_to_ssd = compute_region_ssd(attrs, &state.region_members[to_region]);

    // Create temporary modified member lists
    let new_from: Vec<usize> = state.region_members[from_region]
        .iter()
        .copied()
        .filter(|&m| m != area)
        .collect();
    let mut new_to = state.region_members[to_region].clone();
    new_to.push(area);

    let new_from_ssd = compute_region_ssd(attrs, &new_from);
    let new_to_ssd = compute_region_ssd(attrs, &new_to);

    (new_from_ssd + new_to_ssd) - (old_from_ssd + old_to_ssd)
}

/// Execute a move
fn execute_move(
    state: &mut RegionState,
    area: usize,
    from_region: usize,
    to_region: usize,
    threshold_var: &[f64],
) {
    // Remove from old region
    state.region_members[from_region].retain(|&m| m != area);
    state.region_thresholds[from_region] -= threshold_var[area];

    // Add to new region
    state.labels[area] = to_region as i32;
    state.region_members[to_region].push(area);
    state.region_thresholds[to_region] += threshold_var[area];
}

/// Compute total objective (sum of within-region SSD)
fn compute_objective(attrs: &[Vec<f64>], state: &RegionState) -> f64 {
    state.region_members
        .iter()
        .map(|members| compute_region_ssd(attrs, members))
        .sum()
}

/// Compute SSD for a single region
fn compute_region_ssd(attrs: &[Vec<f64>], members: &[usize]) -> f64 {
    if members.is_empty() {
        return 0.0;
    }

    let n_attrs = attrs[0].len();
    let n = members.len() as f64;

    // Compute centroid
    let mut centroid = vec![0.0; n_attrs];
    for &m in members {
        for (k, val) in attrs[m].iter().enumerate() {
            centroid[k] += val;
        }
    }
    for c in &mut centroid {
        *c /= n;
    }

    // Compute SSD
    let mut ssd = 0.0;
    for &m in members {
        for (k, val) in attrs[m].iter().enumerate() {
            let diff = val - centroid[k];
            ssd += diff * diff;
        }
    }

    ssd
}

// ============================================================================
// Compactness computation functions
// ============================================================================

/// Compute compactness for a region using centroid dispersion
///
/// Returns a value between 0 and 1, where 1 is most compact (all units at same location)
/// Uses the formula: C = 1 / (1 + sqrt(mean_squared_distance_to_centroid))
fn compute_region_centroid_dispersion(members: &[usize], centroids: &Centroids) -> f64 {
    if members.len() <= 1 {
        return 1.0; // Single unit is maximally compact
    }

    // Compute geographic centroid of the region
    let n = members.len() as f64;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for &m in members {
        cx += centroids.x[m];
        cy += centroids.y[m];
    }
    cx /= n;
    cy /= n;

    // Compute mean squared distance to centroid
    let mut sum_sq_dist = 0.0;
    for &m in members {
        let dx = centroids.x[m] - cx;
        let dy = centroids.y[m] - cy;
        sum_sq_dist += dx * dx + dy * dy;
    }
    let mean_sq_dist = sum_sq_dist / n;

    // Convert to compactness (0-1 scale)
    // Higher mean_sq_dist = less compact = lower score
    1.0 / (1.0 + mean_sq_dist.sqrt())
}

/// Compute compactness for a region using normalized moment of inertia (NMI)
/// 
/// Returns a value between 0 and 1, where 1 is most compact
/// Uses the formula C = (area^2) / (2 * pi * I) where I is the second areal moment
fn compute_region_nmi(
    members: &[usize],
    areas: &[f64],
    moments: &[f64],
    centroids: &Centroids,
) -> f64 {
    let n = members.len();
    
    // Validate inputs
    if n == 0 {
        return 0.0;
    }
    
    // Calculate the combined centroid as area-weighted average
    let mut total_area = 0.0;
    let mut combined_x = 0.0;
    let mut combined_y = 0.0;
    
    for &m in members {
        total_area += areas[m];
        combined_x += areas[m] * centroids.x[m];
        combined_y += areas[m] * centroids.y[m];
    }
    
    if total_area <= 0.0 {
        return 0.0;
    }
    
    combined_x /= total_area;
    combined_y /= total_area;
    
    // Calculate distance from each centroid to combined centroid
    // and apply parallel axis theorem: I_combined = sum(I_i + A_i * d_i^2)
    let mut combined_moment = 0.0;
    
    for &m in members {
        let dx = centroids.x[m] - combined_x;
        let dy = centroids.y[m] - combined_y;
        let distance_squared = dx * dx + dy * dy;
        
        combined_moment += moments[m] + areas[m] * distance_squared;
    }
    
    // Calculate total area for members
    let total_area: f64 = members.iter().map(|&m| areas[m]).sum();
    
    if combined_moment <= 0.0 {
        return 0.0;
    }
    
    let pi = std::f64::consts::PI;
    (total_area * total_area) / (2.0 * pi * combined_moment)
}

/// Compute compactness for a region using a specified metric
///
/// Returns a value between 0 and 1, where 1 is most compact (all units at same location)
/// Uses the formula: C = 1 / (1 + sqrt(mean_squared_distance_to_centroid))
fn compute_region_compactness(
    members: &[usize],
    centroids: &Centroids,
    metric: CompactnessMetric,
    areas: Option<&[f64]>,
    moments: Option<&[f64]>,
) -> f64 {
    match metric {
        CompactnessMetric::CentroidDispersion => {
            compute_region_centroid_dispersion(members, centroids)
        }
        CompactnessMetric::NMI => {
            match (areas, moments) {
                (Some(a), Some(m)) => compute_region_nmi(members, a, m, centroids),
                _ => {
                    panic!("NMI compactness metric requires both areas and moments to be provided");
                }
            }
        }
    }
}

/// Compute region centroid dispersion if a new member were added to a region
fn compute_region_centroid_dispersion_if_added(
    state: &RegionState,
    region_id: usize,
    new_member: usize,
    centroids: &Centroids,
) -> f64 {
    let members = &state.region_members[region_id];

    // Create temporary member list with new member
    let n = (members.len() + 1) as f64;

    // Compute new centroid
    let mut cx = centroids.x[new_member];
    let mut cy = centroids.y[new_member];
    for &m in members {
        cx += centroids.x[m];
        cy += centroids.y[m];
    }
    cx /= n;
    cy /= n;

    // Compute mean squared distance to new centroid
    let mut sum_sq_dist = 0.0;
    let dx = centroids.x[new_member] - cx;
    let dy = centroids.y[new_member] - cy;
    sum_sq_dist += dx * dx + dy * dy;

    for &m in members {
        let dx = centroids.x[m] - cx;
        let dy = centroids.y[m] - cy;
        sum_sq_dist += dx * dx + dy * dy;
    }
    let mean_sq_dist = sum_sq_dist / n;

    1.0 / (1.0 + mean_sq_dist.sqrt())
}

/// Compute NMI if a new member were added to a region
fn compute_region_nmi_if_added(
    state: &RegionState,
    region_id: usize,
    new_member: usize,
    areas: &[f64],
    moments: &[f64],
    centroids: &Centroids,
) -> f64 {
    let members = &state.region_members[region_id];

    // Create temporary member list with new member
    let mut temp_members = members.clone();
    temp_members.push(new_member);

    // Compute NMI for the temporary region
    compute_region_nmi(&temp_members, areas, moments, centroids)
}

/// Compute compactness if a new member were added to a region
fn compute_region_compactness_if_added(
    state: &RegionState,
    region_id: usize,
    new_member: usize,
    centroids: &Centroids,
    metric: CompactnessMetric,
    areas: Option<&[f64]>,
    moments: Option<&[f64]>,
) -> f64 {
    match metric {
        CompactnessMetric::CentroidDispersion => {
            compute_region_centroid_dispersion_if_added(state, region_id, new_member, centroids)
        }
        CompactnessMetric::NMI => {
            match (areas, moments) {
                (Some(a), Some(m)) => {
                    compute_region_nmi_if_added(state, region_id, new_member, a, m, centroids)
                }
                _ => {
                    panic!("NMI compactness metric requires both areas and moments to be provided");
                }
            }
        }
    }
}

/// Compute compactness for all regions
fn compute_all_region_compactness(
    state: &RegionState,
    centroids: &Centroids,
    metric: CompactnessMetric,
    areas: Option<&[f64]>,
    moments: Option<&[f64]>
) -> Vec<f64> {
    state.region_members
        .iter()
        .map(|members| compute_region_compactness(members, centroids, metric, areas, moments))
        .collect()
}
