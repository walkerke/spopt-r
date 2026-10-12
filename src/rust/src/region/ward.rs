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

use extendr_api::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

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

    // Cluster state, indexed by cluster id (initially one cluster per area).
    let mut size: Vec<f64> = vec![1.0; n];
    let mut sum: Vec<Vec<f64>> = attr_vecs.clone();
    let mut alive: Vec<bool> = vec![true; n];
    let mut version: Vec<u32> = vec![0; n];
    let mut adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for (&i, &j) in adj_i.iter().zip(adj_j.iter()) {
        let (i, j) = (i as usize, j as usize);
        if i != j && i < n && j < n {
            adjacency[i].insert(j);
            adjacency[j].insert(i);
        }
    }
    // Which cluster each area belongs to; updated through a union-find style
    // "absorbed by" chain to keep merges O(1).
    let mut owner: Vec<usize> = (0..n).collect();

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

    let mut n_clusters = n;
    while n_clusters > n_regions {
        let cand = match heap.pop() {
            Some(c) => c,
            None => break, // disconnected graph: no adjacent pairs left
        };
        let (c1, c2) = (cand.c1, cand.c2);
        if !alive[c1] || !alive[c2] || version[c1] != cand.v1 || version[c2] != cand.v2 {
            continue; // stale
        }

        // Merge c2 into c1.
        alive[c2] = false;
        owner[c2] = c1;
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
            }
        }
        adjacency[c1].remove(&c2);
        version[c1] += 1;
        n_clusters -= 1;

        // Fresh candidates between the merged cluster and its neighbours.
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

    // Resolve each area to its surviving cluster, then relabel 1..k in order
    // of first appearance.
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
    let actual_n_regions = next as usize;

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
