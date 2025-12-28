//! SKATER (Spatial 'K'luster Analysis by Tree Edge Removal)
//!
//! Implements the SKATER algorithm from:
//! Assuncao et al. (2006) "Efficient regionalization techniques for socio-economic
//! geographical units using minimum spanning trees"

use extendr_api::prelude::*;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::algo::min_spanning_tree;
use petgraph::data::FromElements;
use petgraph::visit::EdgeRef;
use std::collections::{HashMap, HashSet};

use crate::graph::cluster_ssd;

/// Solve SKATER regionalization
pub fn solve(
    attrs: RMatrix<f64>,
    adj_i: &[i32],
    adj_j: &[i32],
    n_regions: usize,
    floor_var: Option<Vec<f64>>,
    floor_value: f64,
    _seed: Option<u64>,
) -> Vec<i32> {
    let n = attrs.nrows();
    let n_attrs = attrs.ncols();

    // Convert R matrix to Vec<Vec<f64>> (row-major)
    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    // Build dissimilarity-weighted graph from adjacency
    let dissim = compute_edge_dissimilarity(&attr_vecs, adj_i, adj_j);

    // Compute MST
    let mut graph: UnGraph<usize, f64> = UnGraph::new_undirected();
    let nodes: Vec<NodeIndex> = (0..n).map(|i| graph.add_node(i)).collect();

    for (idx, (&row, &col)) in adj_i.iter().zip(adj_j.iter()).enumerate() {
        if row < col as i32 {
            graph.add_edge(nodes[row as usize], nodes[col as usize], dissim[idx]);
        }
    }

    let mst: UnGraph<usize, f64> = UnGraph::from_elements(min_spanning_tree(&graph));

    // Iteratively cut edges to create regions
    let mut current_labels: Vec<i32> = vec![0; n];
    let mut current_n_regions = 1;
    let mut mst_edges: Vec<(usize, usize, f64)> = mst
        .edge_references()
        .map(|e| (e.source().index(), e.target().index(), *e.weight()))
        .collect();

    while current_n_regions < n_regions && !mst_edges.is_empty() {
        // Find best edge to cut
        let best_cut = find_best_cut(
            &attr_vecs,
            &current_labels,
            &mst_edges,
            floor_var.as_deref(),
            floor_value,
        );

        if let Some((cut_idx, new_labels)) = best_cut {
            // Remove the cut edge
            mst_edges.remove(cut_idx);

            // Update labels
            current_labels = new_labels;
            current_n_regions += 1;
        } else {
            // No valid cut found (floor constraint prevents further cuts)
            break;
        }
    }

    // Convert to 1-based labels for R
    current_labels.iter().map(|&l| l + 1).collect()
}

/// Compute edge dissimilarities based on attribute differences
fn compute_edge_dissimilarity(
    attrs: &[Vec<f64>],
    adj_i: &[i32],
    adj_j: &[i32],
) -> Vec<f64> {
    adj_i
        .iter()
        .zip(adj_j.iter())
        .map(|(&i, &j)| {
            let i = i as usize;
            let j = j as usize;
            // Euclidean distance in attribute space
            attrs[i]
                .iter()
                .zip(attrs[j].iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect()
}

/// Find the best edge to cut that minimizes total SSD
fn find_best_cut(
    attrs: &[Vec<f64>],
    current_labels: &[i32],
    edges: &[(usize, usize, f64)],
    floor_var: Option<&[f64]>,
    floor_value: f64,
) -> Option<(usize, Vec<i32>)> {
    let n = attrs.len();
    let mut best_cut_idx: Option<usize> = None;
    let mut best_ssd = f64::MAX;
    let mut best_labels: Vec<i32> = Vec::new();

    // Build adjacency from current MST edges for component detection
    for (edge_idx, &(from, to, _)) in edges.iter().enumerate() {
        // Simulate cutting this edge
        let remaining_edges: Vec<_> = edges
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != edge_idx)
            .map(|(_, e)| e)
            .collect();

        // Find new components after cut
        let new_labels = compute_components_after_cut(n, &remaining_edges, current_labels);

        // Check floor constraint
        if let Some(fv) = floor_var {
            let floor_satisfied = check_floor_constraint(&new_labels, fv, floor_value);
            if !floor_satisfied {
                continue;
            }
        }

        // Compute SSD for this cut
        let ssd = compute_total_ssd(attrs, &new_labels);

        if ssd < best_ssd {
            best_ssd = ssd;
            best_cut_idx = Some(edge_idx);
            best_labels = new_labels;
        }
    }

    best_cut_idx.map(|idx| (idx, best_labels))
}

/// Compute connected components after removing an edge
fn compute_components_after_cut(
    n: usize,
    edges: &[&(usize, usize, f64)],
    current_labels: &[i32],
) -> Vec<i32> {
    // Union-find
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], x: usize, y: usize) {
        let rx = find(parent, x);
        let ry = find(parent, y);
        if rx != ry {
            parent[ry] = rx;
        }
    }

    for &&(from, to, _) in edges {
        union(&mut parent, from, to);
    }

    // Assign new labels
    let mut label_map: HashMap<usize, i32> = HashMap::new();
    let max_existing = current_labels.iter().max().copied().unwrap_or(0);
    let mut next_label = max_existing + 1;

    let mut result = vec![0i32; n];
    for i in 0..n {
        let root = find(&mut parent.clone(), i);
        let label = *label_map.entry(root).or_insert_with(|| {
            let l = next_label;
            next_label += 1;
            l
        });
        result[i] = label;
    }

    // Relabel to be contiguous from 0
    let unique_labels: HashSet<_> = result.iter().copied().collect();
    let label_remap: HashMap<i32, i32> = unique_labels
        .into_iter()
        .enumerate()
        .map(|(i, l)| (l, i as i32))
        .collect();

    result.iter().map(|l| label_remap[l]).collect()
}

/// Check if floor constraint is satisfied for all regions
fn check_floor_constraint(labels: &[i32], floor_var: &[f64], floor_value: f64) -> bool {
    let mut region_sums: HashMap<i32, f64> = HashMap::new();

    for (&label, &val) in labels.iter().zip(floor_var.iter()) {
        *region_sums.entry(label).or_insert(0.0) += val;
    }

    region_sums.values().all(|&sum| sum >= floor_value)
}

/// Compute total SSD across all regions
fn compute_total_ssd(attrs: &[Vec<f64>], labels: &[i32]) -> f64 {
    // Group indices by region
    let mut region_indices: HashMap<i32, Vec<usize>> = HashMap::new();
    for (i, &label) in labels.iter().enumerate() {
        region_indices.entry(label).or_default().push(i);
    }

    // Sum SSD across regions
    region_indices
        .values()
        .map(|indices| cluster_ssd(attrs, indices))
        .sum()
}
