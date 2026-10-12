//! SKATER (Spatial 'K'luster Analysis by Tree Edge Removal)
//!
//! Implements the SKATER algorithm from:
//! Assuncao et al. (2006) "Efficient regionalization techniques for socio-economic
//! geographical units using minimum spanning trees"
//!
//! The minimum spanning tree is cut k-1 times; each cut removes the tree edge
//! whose removal lowers the total within-region sum of squared deviations the
//! most, subject to an optional floor on a constraint variable in both halves.
//!
//! Evaluating a cut only needs the sums of the two halves it creates. Rooting
//! each region's subtree and accumulating per-subtree sums (count, attribute
//! sums, attribute sums of squares, floor sum) in one post-order pass gives
//! every candidate edge in that region in O(|region| * p). A region's best cut
//! is cached and only the two regions produced by a cut are re-evaluated, so
//! the whole procedure is O(k * n * p) instead of re-scanning every edge and
//! every region for each candidate.

use extendr_api::prelude::*;
use petgraph::algo::min_spanning_tree;
use petgraph::data::FromElements;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;

/// Best feasible cut inside one region.
#[derive(Clone, Copy)]
struct Cut {
    /// Decrease in total SSD if the cut is made.
    gain: f64,
    /// Index of the MST edge (its position in the original edge list).
    edge: usize,
    /// Node whose subtree (away from the region root) becomes the new region.
    child: usize,
}

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
    if n == 0 {
        return Vec::new();
    }

    // Convert R matrix to Vec<Vec<f64>> (row-major)
    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    let mst_edges = mst_edges(&attr_vecs, adj_i, adj_j);
    let labels = cut_tree(&attr_vecs, &mst_edges, n_regions, floor_var.as_deref(), floor_value);

    // 1-based labels for R
    labels.iter().map(|&l| l + 1).collect()
}

/// Minimum spanning tree of the contiguity graph weighted by attribute
/// Euclidean distance (Kruskal, via petgraph). Shared by SKATER and
/// REDCAP first-order single linkage, whose tree this is.
pub fn mst_edges(attr_vecs: &[Vec<f64>], adj_i: &[i32], adj_j: &[i32]) -> Vec<(usize, usize)> {
    let n = attr_vecs.len();
    let dissim = compute_edge_dissimilarity(attr_vecs, adj_i, adj_j);
    let mut graph: UnGraph<usize, f64> = UnGraph::new_undirected();
    let nodes: Vec<NodeIndex> = (0..n).map(|i| graph.add_node(i)).collect();
    for (idx, (&row, &col)) in adj_i.iter().zip(adj_j.iter()).enumerate() {
        if row < col {
            graph.add_edge(nodes[row as usize], nodes[col as usize], dissim[idx]);
        }
    }
    let mst: UnGraph<usize, f64> = UnGraph::from_elements(min_spanning_tree(&graph));
    mst.edge_references()
        .map(|e| (e.source().index(), e.target().index()))
        .collect()
}

/// Partition a spanning tree (or forest) over `n` areas into `n_regions`
/// regions by repeatedly removing the edge whose removal lowers the total
/// within-region SSD the most, subject to an optional floor on both halves.
/// Returns 0-based labels. Shared by SKATER (minimum spanning tree) and
/// REDCAP (tree from contiguity-constrained hierarchical clustering).
pub fn cut_tree(
    attr_vecs: &[Vec<f64>],
    tree_edges: &[(usize, usize)],
    n_regions: usize,
    floor: Option<&[f64]>,
    floor_value: f64,
) -> Vec<i32> {
    let n = attr_vecs.len();
    let n_attrs = attr_vecs.first().map(|a| a.len()).unwrap_or(0);
    let mst_edges = tree_edges;
    // Tree adjacency: (neighbour, edge index)
    let mut tree: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (idx, &(a, b)) in mst_edges.iter().enumerate() {
        tree[a].push((b, idx));
        tree[b].push((a, idx));
    }
    let mut cut_edge = vec![false; mst_edges.len()];

    let mut scratch = Scratch::new(n, n_attrs);

    // Initial regions: one per connected component of the tree (a connected
    // graph gives exactly one). Labels are 0-based region ids.
    let mut labels: Vec<i32> = vec![-1; n];
    let mut roots: Vec<usize> = Vec::new();
    for v in 0..n {
        if labels[v] < 0 {
            let r = roots.len() as i32;
            labels[v] = r;
            relabel_component(v, r, &tree, &cut_edge, &mut labels);
            roots.push(v);
        }
    }
    let mut best: Vec<Option<Cut>> = (0..roots.len())
        .map(|r| {
            evaluate_region(
                r, roots[r], attr_vecs, &tree, &cut_edge, &labels, floor, floor_value, &mut scratch,
            )
        })
        .collect();

    while roots.len() < n_regions {
        // Best cut over all regions: largest gain, ties to the earliest edge.
        let mut pick: Option<(usize, Cut)> = None;
        for (r, c) in best.iter().enumerate() {
            if let Some(c) = c {
                let better = match pick {
                    None => true,
                    Some((_, p)) => c.gain > p.gain || (c.gain == p.gain && c.edge < p.edge),
                };
                if better {
                    pick = Some((r, *c));
                }
            }
        }
        let (r, cut) = match pick {
            Some(p) => p,
            None => break, // no feasible cut anywhere
        };

        // Make the cut: the child's side becomes a new region.
        cut_edge[cut.edge] = true;
        let new_r = roots.len();
        relabel_component(cut.child, new_r as i32, &tree, &cut_edge, &mut labels);
        roots.push(cut.child);
        best[r] = evaluate_region(
            r, roots[r], attr_vecs, &tree, &cut_edge, &labels, floor, floor_value, &mut scratch,
        );
        best.push(evaluate_region(
            new_r, cut.child, attr_vecs, &tree, &cut_edge, &labels, floor, floor_value, &mut scratch,
        ));
    }

    labels
}


/// Compute edge dissimilarities based on attribute differences
pub fn compute_edge_dissimilarity(attrs: &[Vec<f64>], adj_i: &[i32], adj_j: &[i32]) -> Vec<f64> {
    adj_i
        .iter()
        .zip(adj_j.iter())
        .map(|(&i, &j)| {
            let (i, j) = (i as usize, j as usize);
            attrs[i]
                .iter()
                .zip(attrs[j].iter())
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt()
        })
        .collect()
}

/// Reusable buffers for region evaluation.
struct Scratch {
    order: Vec<usize>,
    parent: Vec<usize>,
    parent_edge: Vec<usize>,
    count: Vec<f64>,
    sum: Vec<Vec<f64>>,
    sumsq: Vec<Vec<f64>>,
    fsum: Vec<f64>,
    visited_stamp: Vec<u32>,
    stamp: u32,
    p: usize,
}

impl Scratch {
    fn new(n: usize, p: usize) -> Self {
        Self {
            order: Vec::with_capacity(n),
            parent: vec![usize::MAX; n],
            parent_edge: vec![usize::MAX; n],
            count: vec![0.0; n],
            sum: vec![vec![0.0; p]; n],
            sumsq: vec![vec![0.0; p]; n],
            fsum: vec![0.0; n],
            visited_stamp: vec![0; n],
            stamp: 0,
            p,
        }
    }
}

/// Relabel every node reachable from `start` over uncut tree edges.
fn relabel_component(
    start: usize,
    label: i32,
    tree: &[Vec<(usize, usize)>],
    cut_edge: &[bool],
    labels: &mut [i32],
) {
    labels[start] = label;
    let mut stack = vec![start];
    while let Some(v) = stack.pop() {
        for &(w, e) in &tree[v] {
            if !cut_edge[e] && labels[w] != label {
                labels[w] = label;
                stack.push(w);
            }
        }
    }
}

/// Within-region SSD from count, sums and sums of squares:
/// sum_j (Q_j - S_j^2 / n).
#[inline]
fn ssd_from_sums(count: f64, sum: &[f64], sumsq: &[f64]) -> f64 {
    if count <= 0.0 {
        return 0.0;
    }
    sum.iter()
        .zip(sumsq.iter())
        .map(|(s, q)| q - s * s / count)
        .sum()
}

/// Find the best feasible cut inside region `r`, rooted at `root`.
///
/// One traversal from the root records a parent array and a pre-order; a
/// reverse pass accumulates subtree sums. The two halves of cutting the edge
/// above node v are v's subtree and the remainder of the region.
#[allow(clippy::too_many_arguments)]
fn evaluate_region(
    r: usize,
    root: usize,
    attrs: &[Vec<f64>],
    tree: &[Vec<(usize, usize)>],
    cut_edge: &[bool],
    labels: &[i32],
    floor: Option<&[f64]>,
    floor_value: f64,
    s: &mut Scratch,
) -> Option<Cut> {
    let p = s.p;
    s.stamp = s.stamp.wrapping_add(1);
    if s.stamp == 0 {
        s.visited_stamp.iter_mut().for_each(|x| *x = 0);
        s.stamp = 1;
    }
    s.order.clear();

    // Pre-order traversal restricted to this region over uncut edges.
    let mut stack = vec![root];
    s.visited_stamp[root] = s.stamp;
    s.parent[root] = usize::MAX;
    s.parent_edge[root] = usize::MAX;
    while let Some(v) = stack.pop() {
        s.order.push(v);
        for &(w, e) in &tree[v] {
            if !cut_edge[e] && labels[w] == r as i32 && s.visited_stamp[w] != s.stamp {
                s.visited_stamp[w] = s.stamp;
                s.parent[w] = v;
                s.parent_edge[w] = e;
                stack.push(w);
            }
        }
    }

    // Initialise per-node sums, then accumulate children into parents in
    // reverse pre-order (children always come after their parent).
    for &v in &s.order {
        s.count[v] = 1.0;
        s.fsum[v] = floor.map_or(0.0, |f| f[v]);
        for j in 0..p {
            let x = attrs[v][j];
            s.sum[v][j] = x;
            s.sumsq[v][j] = x * x;
        }
    }
    for idx in (1..s.order.len()).rev() {
        let v = s.order[idx];
        let u = s.parent[v];
        s.count[u] += s.count[v];
        s.fsum[u] += s.fsum[v];
        for j in 0..p {
            s.sum[u][j] += s.sum[v][j];
            s.sumsq[u][j] += s.sumsq[v][j];
        }
    }

    let total_count = s.count[root];
    let total_f = s.fsum[root];
    let region_ssd = ssd_from_sums(total_count, &s.sum[root], &s.sumsq[root]);

    let mut best: Option<Cut> = None;
    let mut rest_sum = vec![0.0; p];
    let mut rest_sumsq = vec![0.0; p];
    for &v in s.order.iter().skip(1) {
        if let Some(_f) = floor {
            if s.fsum[v] < floor_value || total_f - s.fsum[v] < floor_value {
                continue;
            }
        }
        let rest_count = total_count - s.count[v];
        for j in 0..p {
            rest_sum[j] = s.sum[root][j] - s.sum[v][j];
            rest_sumsq[j] = s.sumsq[root][j] - s.sumsq[v][j];
        }
        let gain = region_ssd
            - ssd_from_sums(s.count[v], &s.sum[v], &s.sumsq[v])
            - ssd_from_sums(rest_count, &rest_sum, &rest_sumsq);
        let e = s.parent_edge[v];
        let better = match best {
            None => true,
            Some(b) => gain > b.gain || (gain == b.gain && e < b.edge),
        };
        if better {
            best = Some(Cut { gain, edge: e, child: v });
        }
    }
    best
}
