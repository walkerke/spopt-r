use extendr_api::prelude::*;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::algo::min_spanning_tree;
use petgraph::data::FromElements;
use std::collections::HashMap;

/// Compute minimum spanning tree
pub fn minimum_spanning_tree(i: &[i32], j: &[i32], weights: &[f64], n: usize) -> List {
    // Build undirected weighted graph
    let mut graph: UnGraph<(), f64> = UnGraph::new_undirected();

    // Add nodes
    let nodes: Vec<NodeIndex> = (0..n).map(|_| graph.add_node(())).collect();

    // Add edges (only add each edge once for undirected graph)
    for (idx, (&row, &col)) in i.iter().zip(j.iter()).enumerate() {
        if row < col as i32 {
            graph.add_edge(nodes[row as usize], nodes[col as usize], weights[idx]);
        }
    }

    // Compute MST
    let mst: UnGraph<(), f64> = UnGraph::from_elements(min_spanning_tree(&graph));

    // Extract MST edges
    let mut mst_from: Vec<i32> = Vec::new();
    let mut mst_to: Vec<i32> = Vec::new();
    let mut mst_weights: Vec<f64> = Vec::new();

    for edge_idx in mst.edge_indices() {
        if let Some((src, tgt)) = mst.edge_endpoints(edge_idx) {
            if let Some(&weight) = mst.edge_weight(edge_idx) {
                mst_from.push(src.index() as i32);
                mst_to.push(tgt.index() as i32);
                mst_weights.push(weight);
            }
        }
    }

    list!(
        from = mst_from,
        to = mst_to,
        weight = mst_weights
    )
}

/// Check if graph is connected
pub fn is_connected(i: &[i32], j: &[i32], n: usize) -> bool {
    if n == 0 {
        return true;
    }
    if i.is_empty() {
        return n == 1;
    }

    let components = connected_components(i, j, n);
    let first = components[0];
    components.iter().all(|&c| c == first)
}

/// Find connected components using union-find
pub fn connected_components(i: &[i32], j: &[i32], n: usize) -> Vec<i32> {
    if n == 0 {
        return Vec::new();
    }

    // Union-find data structure
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank: Vec<usize> = vec![0; n];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]); // Path compression
        }
        parent[x]
    }

    fn union(parent: &mut [usize], rank: &mut [usize], x: usize, y: usize) {
        let root_x = find(parent, x);
        let root_y = find(parent, y);

        if root_x != root_y {
            // Union by rank
            if rank[root_x] < rank[root_y] {
                parent[root_x] = root_y;
            } else if rank[root_x] > rank[root_y] {
                parent[root_y] = root_x;
            } else {
                parent[root_y] = root_x;
                rank[root_x] += 1;
            }
        }
    }

    // Process edges
    for (&row, &col) in i.iter().zip(j.iter()) {
        union(&mut parent, &mut rank, row as usize, col as usize);
    }

    // Get final roots for all nodes
    let final_parents: Vec<usize> = (0..n).map(|idx| find(&mut parent, idx)).collect();

    // Assign component labels
    let mut label_map: HashMap<usize, i32> = HashMap::new();
    let mut next_label = 0i32;

    final_parents
        .iter()
        .map(|&root| {
            *label_map.entry(root).or_insert_with(|| {
                let label = next_label;
                next_label += 1;
                label
            })
        })
        .collect()
}

/// Compute sum of squared deviations for a cluster
pub fn cluster_ssd(attrs: &[Vec<f64>], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }

    let n_attrs = attrs[0].len();
    let n = indices.len() as f64;

    // Compute centroid
    let mut centroid = vec![0.0; n_attrs];
    for &idx in indices {
        for (k, val) in attrs[idx].iter().enumerate() {
            centroid[k] += val;
        }
    }
    for c in &mut centroid {
        *c /= n;
    }

    // Compute SSD
    let mut ssd = 0.0;
    for &idx in indices {
        for (k, val) in attrs[idx].iter().enumerate() {
            let diff = val - centroid[k];
            ssd += diff * diff;
        }
    }

    ssd
}
