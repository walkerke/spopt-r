//! SPENC (Spatially-Encouraged Spectral Clustering)
//!
//! Implements spectral clustering with spatial constraints by combining
//! spatial weights with attribute affinity using kernel methods.
//!
//! Based on: Wolf, L.J. (2020) "Spatially-encouraged spectral clustering"

use extendr_api::prelude::*;
use nalgebra::{DMatrix, SymmetricEigen};
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;

/// Solve SPENC regionalization
///
/// Returns labels (1-based for R), n_regions, and objective value
pub fn solve(
    attrs: RMatrix<f64>,
    n_regions: usize,
    adj_i: &[i32],
    adj_j: &[i32],
    gamma: f64,
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

    // Build spatial weights matrix (binary adjacency)
    let spatial_weights = build_spatial_weights(n, adj_i, adj_j);

    // Compute attribute affinity using RBF kernel
    let attribute_affinity = compute_rbf_affinity(&attr_vecs, gamma);

    // Combine: affinity = spatial * attribute (element-wise)
    let affinity = combine_affinities(&spatial_weights, &attribute_affinity);

    // Compute normalized Laplacian
    let laplacian = compute_normalized_laplacian(&affinity);

    // Get k smallest eigenvectors (skip first which is trivial)
    let embedding = compute_spectral_embedding(&laplacian, n_regions);

    // Run k-means on embedding
    let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(42));
    let labels = kmeans_clustering(&embedding, n_regions, 100, &mut rng);

    // Compute objective (within-cluster sum of squared distances in embedding space)
    let objective = compute_embedding_objective(&embedding, &labels, n_regions);

    // Convert to 1-based labels
    let labels_1based: Vec<i32> = labels.iter().map(|&l| l as i32 + 1).collect();

    list!(
        labels = labels_1based,
        n_regions = n_regions as i32,
        objective = objective
    )
}

/// Build spatial weights matrix from adjacency indices
fn build_spatial_weights(n: usize, adj_i: &[i32], adj_j: &[i32]) -> DMatrix<f64> {
    let mut w = DMatrix::zeros(n, n);

    for (&i, &j) in adj_i.iter().zip(adj_j.iter()) {
        let i = i as usize;
        let j = j as usize;
        w[(i, j)] = 1.0;
        w[(j, i)] = 1.0;
    }

    w
}

/// Compute RBF (Gaussian) kernel affinity matrix
fn compute_rbf_affinity(attrs: &[Vec<f64>], gamma: f64) -> DMatrix<f64> {
    let n = attrs.len();
    let mut affinity = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in i..n {
            if i == j {
                affinity[(i, j)] = 1.0;
            } else {
                // Squared Euclidean distance
                let sq_dist: f64 = attrs[i]
                    .iter()
                    .zip(attrs[j].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum();

                // RBF kernel: exp(-gamma * ||x - y||^2)
                let similarity = (-gamma * sq_dist).exp();
                affinity[(i, j)] = similarity;
                affinity[(j, i)] = similarity;
            }
        }
    }

    affinity
}

/// Combine spatial and attribute affinities (element-wise multiplication)
fn combine_affinities(spatial: &DMatrix<f64>, attribute: &DMatrix<f64>) -> DMatrix<f64> {
    let n = spatial.nrows();
    let mut combined = DMatrix::zeros(n, n);

    for i in 0..n {
        for j in 0..n {
            // Only keep affinity where spatial connection exists
            combined[(i, j)] = spatial[(i, j)] * attribute[(i, j)];
        }
        // Add self-connection for stability
        combined[(i, i)] = 1.0;
    }

    combined
}

/// Compute normalized Laplacian: L = I - D^(-1/2) * W * D^(-1/2)
fn compute_normalized_laplacian(affinity: &DMatrix<f64>) -> DMatrix<f64> {
    let n = affinity.nrows();

    // Compute degree vector
    let degrees: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| affinity[(i, j)]).sum())
        .collect();

    // Compute D^(-1/2)
    let d_inv_sqrt: Vec<f64> = degrees
        .iter()
        .map(|&d| if d > 1e-10 { 1.0 / d.sqrt() } else { 0.0 })
        .collect();

    // Compute normalized affinity: D^(-1/2) * W * D^(-1/2)
    let mut normalized = DMatrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            normalized[(i, j)] = d_inv_sqrt[i] * affinity[(i, j)] * d_inv_sqrt[j];
        }
    }

    // Laplacian = I - normalized
    let mut laplacian = DMatrix::identity(n, n);
    laplacian -= &normalized;

    laplacian
}

/// Compute spectral embedding (k smallest eigenvectors of Laplacian)
fn compute_spectral_embedding(laplacian: &DMatrix<f64>, k: usize) -> Vec<Vec<f64>> {
    let n = laplacian.nrows();

    // Compute eigendecomposition
    let eigen = SymmetricEigen::new(laplacian.clone());

    // Get eigenvalues and eigenvectors
    let eigenvalues = eigen.eigenvalues;
    let eigenvectors = eigen.eigenvectors;

    // Sort by eigenvalue (ascending)
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        eigenvalues[a].partial_cmp(&eigenvalues[b]).unwrap()
    });

    // Take k smallest eigenvectors (skip the first trivial one if it's near-zero)
    let start_idx = if eigenvalues[indices[0]].abs() < 1e-8 { 1 } else { 0 };
    let selected_indices: Vec<usize> = indices[start_idx..].iter().take(k).copied().collect();

    // Build embedding matrix (n x k)
    let mut embedding: Vec<Vec<f64>> = vec![vec![0.0; selected_indices.len()]; n];

    for (embed_col, &eigen_idx) in selected_indices.iter().enumerate() {
        for i in 0..n {
            embedding[i][embed_col] = eigenvectors[(i, eigen_idx)];
        }
    }

    // Normalize rows to unit length
    for row in &mut embedding {
        let norm: f64 = row.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-10 {
            for x in row.iter_mut() {
                *x /= norm;
            }
        }
    }

    embedding
}

/// K-means clustering on embedding
fn kmeans_clustering(
    embedding: &[Vec<f64>],
    k: usize,
    max_iter: usize,
    rng: &mut ChaCha8Rng,
) -> Vec<usize> {
    let n = embedding.len();
    if n == 0 || k == 0 {
        return vec![];
    }

    let dim = embedding[0].len();

    // Initialize centroids using k-means++
    let mut centroids = kmeans_plusplus_init(embedding, k, rng);

    let mut labels = vec![0usize; n];
    let mut changed = true;
    let mut iter = 0;

    while changed && iter < max_iter {
        changed = false;
        iter += 1;

        // Assign points to nearest centroid
        for i in 0..n {
            let mut best_cluster = 0;
            let mut best_dist = f64::MAX;

            for c in 0..k {
                let dist = euclidean_distance(&embedding[i], &centroids[c]);
                if dist < best_dist {
                    best_dist = dist;
                    best_cluster = c;
                }
            }

            if labels[i] != best_cluster {
                labels[i] = best_cluster;
                changed = true;
            }
        }

        // Update centroids
        let mut new_centroids = vec![vec![0.0; dim]; k];
        let mut counts = vec![0usize; k];

        for i in 0..n {
            let c = labels[i];
            counts[c] += 1;
            for d in 0..dim {
                new_centroids[c][d] += embedding[i][d];
            }
        }

        for c in 0..k {
            if counts[c] > 0 {
                for d in 0..dim {
                    new_centroids[c][d] /= counts[c] as f64;
                }
            } else {
                // Empty cluster: reinitialize randomly
                let rand_idx = rng.gen_range(0..n);
                new_centroids[c] = embedding[rand_idx].clone();
            }
        }

        centroids = new_centroids;
    }

    labels
}

/// K-means++ initialization
fn kmeans_plusplus_init(
    embedding: &[Vec<f64>],
    k: usize,
    rng: &mut ChaCha8Rng,
) -> Vec<Vec<f64>> {
    let n = embedding.len();
    let mut centroids: Vec<Vec<f64>> = Vec::with_capacity(k);

    // First centroid: random point
    let first_idx = rng.gen_range(0..n);
    centroids.push(embedding[first_idx].clone());

    // Remaining centroids: weighted by squared distance to nearest centroid
    for _ in 1..k {
        let mut distances: Vec<f64> = vec![f64::MAX; n];

        for i in 0..n {
            for centroid in &centroids {
                let dist = euclidean_distance(&embedding[i], centroid);
                if dist < distances[i] {
                    distances[i] = dist;
                }
            }
        }

        // Square distances for probability weighting
        let sq_distances: Vec<f64> = distances.iter().map(|d| d * d).collect();
        let total: f64 = sq_distances.iter().sum();

        if total < 1e-10 {
            // All points are at existing centroids, pick random
            let rand_idx = rng.gen_range(0..n);
            centroids.push(embedding[rand_idx].clone());
        } else {
            // Weighted random selection
            let threshold = rng.gen::<f64>() * total;
            let mut cumsum = 0.0;
            let mut selected = 0;

            for (i, &sq_dist) in sq_distances.iter().enumerate() {
                cumsum += sq_dist;
                if cumsum >= threshold {
                    selected = i;
                    break;
                }
            }

            centroids.push(embedding[selected].clone());
        }
    }

    centroids
}

/// Euclidean distance between two vectors
fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Compute objective: total within-cluster sum of squared distances
fn compute_embedding_objective(
    embedding: &[Vec<f64>],
    labels: &[usize],
    k: usize,
) -> f64 {
    let n = embedding.len();
    let dim = if n > 0 { embedding[0].len() } else { 0 };

    // Compute centroids
    let mut centroids = vec![vec![0.0; dim]; k];
    let mut counts = vec![0usize; k];

    for i in 0..n {
        let c = labels[i];
        counts[c] += 1;
        for d in 0..dim {
            centroids[c][d] += embedding[i][d];
        }
    }

    for c in 0..k {
        if counts[c] > 0 {
            for d in 0..dim {
                centroids[c][d] /= counts[c] as f64;
            }
        }
    }

    // Compute total squared distance
    let mut total = 0.0;
    for i in 0..n {
        let c = labels[i];
        total += euclidean_distance(&embedding[i], &centroids[c]).powi(2);
    }

    total
}
