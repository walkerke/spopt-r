//! Spatially-Constrained Ward Hierarchical Clustering
//!
//! Implements Ward's minimum variance clustering with spatial contiguity constraint.
//! Only clusters that are spatially adjacent can be merged.

use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};
use std::cmp::Ordering;

/// A merge candidate with Ward distance (reserved for future priority queue optimization)
#[derive(Clone)]
#[allow(dead_code)]
struct MergeCandidate {
    cluster1: usize,
    cluster2: usize,
    ward_distance: f64,
}

impl PartialEq for MergeCandidate {
    fn eq(&self, other: &Self) -> bool {
        self.ward_distance == other.ward_distance
    }
}

impl Eq for MergeCandidate {}

impl PartialOrd for MergeCandidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Reverse ordering for min-heap behavior
        other.ward_distance.partial_cmp(&self.ward_distance)
    }
}

impl Ord for MergeCandidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

/// Cluster state
struct ClusterState {
    /// Members of each cluster
    members: HashMap<usize, Vec<usize>>,
    /// Centroid of each cluster (sum of attributes)
    centroids: HashMap<usize, Vec<f64>>,
    /// Which cluster each observation belongs to
    labels: Vec<usize>,
    /// Adjacent clusters for each cluster
    adjacency: HashMap<usize, HashSet<usize>>,
    /// Number of active clusters
    n_clusters: usize,
}

impl ClusterState {
    fn new(n: usize, attrs: &[Vec<f64>], adj_i: &[i32], adj_j: &[i32]) -> Self {
        // Initialize each observation as its own cluster
        let mut members = HashMap::new();
        let mut centroids = HashMap::new();
        let labels: Vec<usize> = (0..n).collect();

        for i in 0..n {
            members.insert(i, vec![i]);
            centroids.insert(i, attrs[i].clone());
        }

        // Build initial adjacency
        let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
        for i in 0..n {
            adjacency.insert(i, HashSet::new());
        }

        for (&i, &j) in adj_i.iter().zip(adj_j.iter()) {
            let i = i as usize;
            let j = j as usize;
            adjacency.get_mut(&i).unwrap().insert(j);
            adjacency.get_mut(&j).unwrap().insert(i);
        }

        Self {
            members,
            centroids,
            labels,
            adjacency,
            n_clusters: n,
        }
    }

    /// Get cluster size
    fn cluster_size(&self, cluster_id: usize) -> usize {
        self.members.get(&cluster_id).map(|m| m.len()).unwrap_or(0)
    }

    /// Compute Ward distance between two clusters
    fn ward_distance(&self, c1: usize, c2: usize, n_attrs: usize) -> f64 {
        let n1 = self.cluster_size(c1) as f64;
        let n2 = self.cluster_size(c2) as f64;

        if n1 == 0.0 || n2 == 0.0 {
            return f64::MAX;
        }

        let cent1 = self.centroids.get(&c1).unwrap();
        let cent2 = self.centroids.get(&c2).unwrap();

        // Squared distance between centroids
        let sq_dist: f64 = (0..n_attrs)
            .map(|k| {
                let mean1 = cent1[k] / n1;
                let mean2 = cent2[k] / n2;
                (mean1 - mean2).powi(2)
            })
            .sum();

        // Ward's criterion: (n1 * n2 / (n1 + n2)) * ||c1 - c2||^2
        (n1 * n2 / (n1 + n2)) * sq_dist
    }

    /// Merge two clusters
    fn merge(&mut self, c1: usize, c2: usize) {
        // c1 will absorb c2
        let members2 = self.members.remove(&c2).unwrap_or_default();
        let centroid2 = self.centroids.remove(&c2).unwrap_or_default();
        let adj2 = self.adjacency.remove(&c2).unwrap_or_default();

        // Update members
        if let Some(members1) = self.members.get_mut(&c1) {
            for m in &members2 {
                members1.push(*m);
                self.labels[*m] = c1;
            }
        }

        // Update centroid (sum)
        if let Some(centroid1) = self.centroids.get_mut(&c1) {
            for (k, val) in centroid2.iter().enumerate() {
                centroid1[k] += val;
            }
        }

        // Update adjacency
        // c1 gets all neighbors of c2 (except c1 itself)
        if let Some(adj1) = self.adjacency.get_mut(&c1) {
            for &neighbor in &adj2 {
                if neighbor != c1 {
                    adj1.insert(neighbor);
                }
            }
            adj1.remove(&c2);
        }

        // Update all neighbors of c2 to point to c1 instead
        for &neighbor in &adj2 {
            if neighbor != c1 {
                if let Some(neighbor_adj) = self.adjacency.get_mut(&neighbor) {
                    neighbor_adj.remove(&c2);
                    neighbor_adj.insert(c1);
                }
            }
        }

        self.n_clusters -= 1;
    }

    /// Get all adjacent cluster pairs
    fn get_adjacent_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for (&c1, neighbors) in &self.adjacency {
            for &c2 in neighbors {
                if c1 < c2 {
                    pairs.push((c1, c2));
                }
            }
        }
        pairs
    }
}

/// Solve spatially-constrained Ward clustering
pub fn solve(
    attrs: RMatrix<f64>,
    n_regions: usize,
    adj_i: &[i32],
    adj_j: &[i32],
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

    // Initialize cluster state
    let mut state = ClusterState::new(n, &attr_vecs, adj_i, adj_j);

    // Agglomerative clustering until we reach n_regions
    while state.n_clusters > n_regions {
        // Find best merge (minimum Ward distance among adjacent pairs)
        let pairs = state.get_adjacent_pairs();

        if pairs.is_empty() {
            // No more adjacent pairs - graph is disconnected
            break;
        }

        let mut best_merge: Option<(usize, usize, f64)> = None;

        for (c1, c2) in pairs {
            let dist = state.ward_distance(c1, c2, n_attrs);
            if best_merge.is_none() || dist < best_merge.as_ref().unwrap().2 {
                best_merge = Some((c1, c2, dist));
            }
        }

        if let Some((c1, c2, _)) = best_merge {
            state.merge(c1, c2);
        } else {
            break;
        }
    }

    // Relabel to contiguous 1-based labels
    let unique_clusters: Vec<usize> = state.labels.iter().copied().collect::<HashSet<_>>().into_iter().collect();
    let mut label_map: HashMap<usize, i32> = HashMap::new();
    for (new_label, &old_label) in unique_clusters.iter().enumerate() {
        label_map.insert(old_label, new_label as i32 + 1);
    }

    let labels: Vec<i32> = state.labels.iter().map(|&l| label_map[&l]).collect();
    let actual_n_regions = unique_clusters.len();

    // Compute objective (total within-cluster SSD)
    let objective = compute_objective(&attr_vecs, &labels, actual_n_regions);

    list!(
        labels = labels,
        n_regions = actual_n_regions as i32,
        objective = objective
    )
}

/// Compute total within-cluster sum of squared deviations
fn compute_objective(attrs: &[Vec<f64>], labels: &[i32], n_regions: usize) -> f64 {
    let n_attrs = attrs[0].len();
    let mut total_ssd = 0.0;

    for region in 1..=n_regions {
        let members: Vec<usize> = labels
            .iter()
            .enumerate()
            .filter(|(_, &l)| l == region as i32)
            .map(|(i, _)| i)
            .collect();

        if members.is_empty() {
            continue;
        }

        let n = members.len() as f64;

        // Compute centroid
        let mut centroid = vec![0.0; n_attrs];
        for &m in &members {
            for (k, val) in attrs[m].iter().enumerate() {
                centroid[k] += val;
            }
        }
        for c in &mut centroid {
            *c /= n;
        }

        // Compute SSD
        for &m in &members {
            for (k, val) in attrs[m].iter().enumerate() {
                let diff = val - centroid[k];
                total_ssd += diff * diff;
            }
        }
    }

    total_ssd
}
