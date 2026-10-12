//! Spatially-Constrained Ward Hierarchical Clustering
//!
//! Implements Ward's minimum variance clustering with spatial contiguity constraint.
//! Only clusters that are spatially adjacent can be merged.
//!
//! Candidate merges live in a min-heap keyed by Ward distance. Each cluster
//! carries a version number that is bumped when it changes; a heap entry is
//! valid only if both of its clusters still have the version it was pushed
//! with, so stale entries are discarded lazily when popped. After a merge only
//! the merged cluster's neighbours get fresh entries. Ties on distance go to
//! the smaller cluster id pair, which makes the result deterministic (the
//! previous full rescan iterated hash maps, whose order varies between runs).
//!
//! The full merge sequence is also the basis of REDCAP (full-order Ward
//! linkage): every merge records the cheapest contiguity edge between the two
//! clusters it joins, and those n-1 edges form the spanning tree that REDCAP
//! then cuts.

use extendr_api::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

/// Heap entry: a candidate merge of two clusters at the versions recorded.
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
    // Reverse order so BinaryHeap behaves as a min-heap on (distance, c1, c2).
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

/// One agglomeration step: cluster `c2` was merged into `c1`, and `edge` is
/// the cheapest contiguity edge (by attribute distance) between them.
#[derive(Clone, Copy)]
pub struct Merge {
    pub c1: usize,
    pub c2: usize,
    pub edge: (usize, usize),
}

/// Ward's criterion: (n1 * n2 / (n1 + n2)) * ||mean1 - mean2||^2
#[inline]
fn ward_distance(n1: f64, sum1: &[f64], n2: f64, sum2: &[f64]) -> f64 {
    let sq: f64 = sum1
        .iter()
        .zip(sum2.iter())
        .map(|(a, b)| {
            let d = a / n1 - b / n2;
            d * d
        })
        .sum();
    (n1 * n2 / (n1 + n2)) * sq
}

#[inline]
fn key(a: usize, b: usize) -> (usize, usize) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Run contiguity-constrained Ward agglomeration to completion and return the
/// merge sequence. `edge_costs[e]` is the attribute distance of contiguity
/// edge `e` (parallel to `adj_i`/`adj_j`); it decides which edge between two
/// clusters is recorded for each merge. A disconnected graph yields fewer than
/// n-1 merges.
pub fn merge_sequence(
    attr_vecs: &[Vec<f64>],
    adj_i: &[i32],
    adj_j: &[i32],
    edge_costs: &[f64],
) -> Vec<Merge> {
    let n = attr_vecs.len();
    let mut size: Vec<f64> = vec![1.0; n];
    let mut sum: Vec<Vec<f64>> = attr_vecs.to_vec();
    let mut alive: Vec<bool> = vec![true; n];
    let mut version: Vec<u32> = vec![0; n];
    let mut adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    // Cheapest contiguity edge between each pair of adjacent clusters.
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
                heap.push(Candidate {
                    distance: ward_distance(size[i], &sum[i], size[j], &sum[j]),
                    c1: i,
                    c2: j,
                    v1: 0,
                    v2: 0,
                });
            }
        }
    }

    let mut merges: Vec<Merge> = Vec::with_capacity(n.saturating_sub(1));
    while let Some(cand) = heap.pop() {
        let (c1, c2) = (cand.c1, cand.c2);
        if !alive[c1] || !alive[c2] || version[c1] != cand.v1 || version[c2] != cand.v2 {
            continue; // stale
        }
        let (_, ea, eb) = best_edge
            .remove(&key(c1, c2))
            .unwrap_or((f64::INFINITY, c1, c2));

        // Merge c2 into c1.
        alive[c2] = false;
        size[c1] += size[c2];
        let sum2 = std::mem::take(&mut sum[c2]);
        for (a, b) in sum[c1].iter_mut().zip(sum2.iter()) {
            *a += b;
        }
        let adj2 = std::mem::take(&mut adjacency[c2]);
        for &nb in &adj2 {
            if nb != c1 {
                adjacency[nb].remove(&c2);
                adjacency[nb].insert(c1);
                adjacency[c1].insert(nb);
                // The cheapest edge between the merged cluster and nb is the
                // cheaper of the two it inherits.
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
        merges.push(Merge { c1, c2, edge: (ea, eb) });

        for &nb in &adjacency[c1] {
            let (a, b) = if c1 < nb { (c1, nb) } else { (nb, c1) };
            heap.push(Candidate {
                distance: ward_distance(size[c1], &sum[c1], size[nb], &sum[nb]),
                c1: a,
                c2: b,
                v1: version[a],
                v2: version[b],
            });
        }
    }
    merges
}

/// Apply the first `n - n_regions` merges and return 1-based labels numbered
/// in order of first appearance.
pub fn labels_from_merges(n: usize, merges: &[Merge], n_regions: usize) -> Vec<i32> {
    let mut owner: Vec<usize> = (0..n).collect();
    fn find(owner: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while owner[r] != r {
            r = owner[r];
        }
        let mut cur = x;
        while owner[cur] != r {
            let next = owner[cur];
            owner[cur] = r;
            cur = next;
        }
        r
    }
    let to_apply = n.saturating_sub(n_regions).min(merges.len());
    for m in &merges[..to_apply] {
        let (r1, r2) = (find(&mut owner, m.c1), find(&mut owner, m.c2));
        if r1 != r2 {
            owner[r2] = r1;
        }
    }
    let mut new_id: Vec<i32> = vec![0; n];
    let mut next = 0i32;
    let mut labels: Vec<i32> = Vec::with_capacity(n);
    for i in 0..n {
        let root = find(&mut owner, i);
        if new_id[root] == 0 {
            next += 1;
            new_id[root] = next;
        }
        labels.push(new_id[root]);
    }
    labels
}

/// Solve spatially-constrained Ward clustering
pub fn solve(attrs: RMatrix<f64>, n_regions: usize, adj_i: &[i32], adj_j: &[i32]) -> List {
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

    // Edge costs only decide which contiguity edge each merge records; the
    // Ward partition itself does not depend on them.
    let edge_costs = vec![0.0; adj_i.len()];
    let merges = merge_sequence(&attr_vecs, adj_i, adj_j, &edge_costs);
    let labels = labels_from_merges(n, &merges, n_regions);
    let actual_n_regions = *labels.iter().max().unwrap_or(&0) as usize;
    let objective = compute_objective(&attr_vecs, &labels, actual_n_regions);

    list!(
        labels = labels,
        n_regions = actual_n_regions as i32,
        objective = objective
    )
}

/// Compute total within-cluster sum of squared deviations (1-based labels)
pub fn compute_objective(attrs: &[Vec<f64>], labels: &[i32], n_regions: usize) -> f64 {
    if attrs.is_empty() || n_regions == 0 {
        return 0.0;
    }
    let n_attrs = attrs[0].len();
    let mut count = vec![0.0f64; n_regions];
    let mut sums = vec![vec![0.0f64; n_attrs]; n_regions];
    for (i, &l) in labels.iter().enumerate() {
        let r = (l - 1) as usize;
        count[r] += 1.0;
        for (k, v) in attrs[i].iter().enumerate() {
            sums[r][k] += v;
        }
    }
    let mut total = 0.0;
    for (i, &l) in labels.iter().enumerate() {
        let r = (l - 1) as usize;
        for (k, v) in attrs[i].iter().enumerate() {
            let d = v - sums[r][k] / count[r];
            total += d * d;
        }
    }
    total
}
