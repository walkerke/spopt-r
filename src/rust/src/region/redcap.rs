//! REDCAP (REgionalization with Dynamically Constrained Agglomerative
//! clustering and Partitioning)
//!
//! Guo, D. (2008). Regionalization with dynamically constrained agglomerative
//! clustering and partitioning. International Journal of Geographical
//! Information Science, 22(7), 801-823.
//!
//! REDCAP runs a contiguity-constrained hierarchical clustering, records for
//! each merge the shortest contiguity edge between the two clusters joined,
//! and treats those n-1 edges as a spanning tree. The tree is then partitioned
//! exactly as in SKATER (`skater::cut_tree`).
//!
//! The five variants exposed by GeoDa, rgeoda and pygeoda are implemented
//! with the same semantics as libgeoda:
//!
//! * `firstorder-singlelinkage`: only contiguity edges are considered and
//!   clusters are joined by the shortest one; the tree is the minimum
//!   spanning tree, i.e. SKATER's tree.
//! * `fullorder-singlelinkage`, `fullorder-averagelinkage`,
//!   `fullorder-completelinkage`: the distance between two clusters is the
//!   minimum, mean or maximum attribute distance over *all* pairs of areas
//!   between them (not only adjacent pairs), but only contiguous clusters may
//!   merge. Cluster distances live in an n x n matrix updated with the
//!   Lance-Williams rules, so these need 8 n^2 bytes.
//! * `fullorder-wardlinkage`: Ward's minimum-variance criterion between
//!   clusters, computed from running sums (`ward::merge_sequence`).
//!
//! libgeoda uses squared Euclidean distances for single and complete linkage
//! and Euclidean distances for average linkage. Minimum and maximum are
//! invariant to the square root, so Euclidean distances are used throughout.

use extendr_api::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

use super::skater::{compute_edge_dissimilarity, cut_tree, mst_edges};
use super::ward::{compute_objective, merge_sequence};

#[derive(Clone, Copy, PartialEq)]
enum Linkage {
    Single,
    Average,
    Complete,
}

#[derive(Clone, Copy, PartialEq)]
struct Candidate {
    distance: f64,
    c1: usize,
    c2: usize,
    v1: u32,
    v2: u32,
}

impl Eq for Candidate {}

impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .partial_cmp(&self.distance)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.c1.cmp(&self.c1))
            .then_with(|| other.c2.cmp(&self.c2))
    }
}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[inline]
fn key(a: usize, b: usize) -> (usize, usize) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Full-order agglomeration with single, average or complete linkage.
/// Returns the recorded tree edges in merge order.
fn fullorder_tree(
    attr_vecs: &[Vec<f64>],
    adj_i: &[i32],
    adj_j: &[i32],
    edge_costs: &[f64],
    linkage: Linkage,
) -> Vec<(usize, usize)> {
    let n = attr_vecs.len();
    // Cluster distance matrix, initialised with pairwise Euclidean distances.
    let mut dist: Vec<f64> = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d: f64 = attr_vecs[i]
                .iter()
                .zip(attr_vecs[j].iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            dist[i * n + j] = d;
            dist[j * n + i] = d;
        }
    }
    let mut size: Vec<f64> = vec![1.0; n];
    let mut alive: Vec<bool> = vec![true; n];
    let mut version: Vec<u32> = vec![0; n];
    let mut adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    let mut best_edge: HashMap<(usize, usize), (f64, usize, usize)> = HashMap::new();
    for (e, (&i, &j)) in adj_i.iter().zip(adj_j.iter()).enumerate() {
        let (i, j) = (i as usize, j as usize);
        if i == j || i >= n || j >= n {
            continue;
        }
        adjacency[i].insert(j);
        adjacency[j].insert(i);
        let k = key(i, j);
        let cost = edge_costs.get(e).copied().unwrap_or(0.0);
        let entry = best_edge.entry(k).or_insert((f64::INFINITY, k.0, k.1));
        if cost < entry.0 {
            *entry = (cost, k.0, k.1);
        }
    }

    let mut heap: BinaryHeap<Candidate> = BinaryHeap::new();
    for i in 0..n {
        for &j in &adjacency[i] {
            if i < j {
                heap.push(Candidate { distance: dist[i * n + j], c1: i, c2: j, v1: 0, v2: 0 });
            }
        }
    }

    let mut tree: Vec<(usize, usize)> = Vec::with_capacity(n.saturating_sub(1));
    while let Some(cand) = heap.pop() {
        let (c1, c2) = (cand.c1, cand.c2);
        if !alive[c1] || !alive[c2] || version[c1] != cand.v1 || version[c2] != cand.v2 {
            continue;
        }
        let (_, ea, eb) = best_edge.remove(&key(c1, c2)).unwrap_or((f64::INFINITY, c1, c2));
        tree.push((ea, eb));

        // Lance-Williams update of the merged cluster's distance to every
        // other live cluster (full order: adjacency does not matter here).
        let (n1, n2) = (size[c1], size[c2]);
        for c in 0..n {
            if !alive[c] || c == c1 || c == c2 {
                continue;
            }
            let d1 = dist[c1 * n + c];
            let d2 = dist[c2 * n + c];
            let d = match linkage {
                Linkage::Single => d1.min(d2),
                Linkage::Complete => d1.max(d2),
                Linkage::Average => (n1 * d1 + n2 * d2) / (n1 + n2),
            };
            dist[c1 * n + c] = d;
            dist[c * n + c1] = d;
        }
        alive[c2] = false;
        size[c1] += size[c2];
        let adj2 = std::mem::take(&mut adjacency[c2]);
        for &nb in &adj2 {
            if nb != c1 {
                adjacency[nb].remove(&c2);
                adjacency[nb].insert(c1);
                adjacency[c1].insert(nb);
                if let Some(from_c2) = best_edge.remove(&key(c2, nb)) {
                    let slot = best_edge
                        .entry(key(c1, nb))
                        .or_insert((f64::INFINITY, from_c2.1, from_c2.2));
                    if from_c2.0 < slot.0 {
                        *slot = from_c2;
                    }
                }
            }
        }
        adjacency[c1].remove(&c2);
        version[c1] += 1;
        for &nb in &adjacency[c1] {
            let (a, b) = if c1 < nb { (c1, nb) } else { (nb, c1) };
            heap.push(Candidate {
                distance: dist[c1 * n + nb],
                c1: a,
                c2: b,
                v1: version[a],
                v2: version[b],
            });
        }
    }
    tree
}

/// Solve REDCAP regionalization
pub fn solve(
    attrs: RMatrix<f64>,
    n_regions: usize,
    adj_i: &[i32],
    adj_j: &[i32],
    method: &str,
    floor_var: Option<Vec<f64>>,
    floor_value: f64,
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

    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    let tree: Vec<(usize, usize)> = match method {
        "firstorder-singlelinkage" => mst_edges(&attr_vecs, adj_i, adj_j),
        "fullorder-wardlinkage" => {
            let costs = compute_edge_dissimilarity(&attr_vecs, adj_i, adj_j);
            merge_sequence(&attr_vecs, adj_i, adj_j, &costs)
                .into_iter()
                .map(|m| m.edge)
                .collect()
        }
        "fullorder-singlelinkage" | "fullorder-averagelinkage" | "fullorder-completelinkage" => {
            let linkage = match method {
                "fullorder-singlelinkage" => Linkage::Single,
                "fullorder-averagelinkage" => Linkage::Average,
                _ => Linkage::Complete,
            };
            let costs = compute_edge_dissimilarity(&attr_vecs, adj_i, adj_j);
            fullorder_tree(&attr_vecs, adj_i, adj_j, &costs, linkage)
        }
        other => panic!(
            "Unknown REDCAP method '{}'; expected one of firstorder-singlelinkage, \
             fullorder-singlelinkage, fullorder-averagelinkage, fullorder-completelinkage, \
             fullorder-wardlinkage",
            other
        ),
    };

    let labels0 = cut_tree(&attr_vecs, &tree, n_regions, floor_var.as_deref(), floor_value);
    let labels: Vec<i32> = labels0.iter().map(|&l| l + 1).collect();
    let actual_n_regions = *labels.iter().max().unwrap_or(&0) as usize;
    let objective = compute_objective(&attr_vecs, &labels, actual_n_regions);

    list!(
        labels = labels,
        n_regions = actual_n_regions as i32,
        objective = objective
    )
}
