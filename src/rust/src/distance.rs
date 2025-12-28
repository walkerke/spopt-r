use extendr_api::prelude::*;
use rayon::prelude::*;

/// Compute Euclidean distance matrix
pub fn euclidean_matrix(x1: &[f64], y1: &[f64], x2: &[f64], y2: &[f64]) -> RMatrix<f64> {
    let n1 = x1.len();
    let n2 = x2.len();

    // Compute distances in parallel for large matrices
    let distances: Vec<f64> = if n1 * n2 > 10000 {
        (0..n1)
            .into_par_iter()
            .flat_map(|i| {
                (0..n2)
                    .map(|j| {
                        let dx = x1[i] - x2[j];
                        let dy = y1[i] - y2[j];
                        (dx * dx + dy * dy).sqrt()
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    } else {
        let mut distances = Vec::with_capacity(n1 * n2);
        for i in 0..n1 {
            for j in 0..n2 {
                let dx = x1[i] - x2[j];
                let dy = y1[i] - y2[j];
                distances.push((dx * dx + dy * dy).sqrt());
            }
        }
        distances
    };

    // Create R matrix (column-major order)
    // extendr expects column-major, so we need to transpose
    let mut col_major = vec![0.0; n1 * n2];
    for i in 0..n1 {
        for j in 0..n2 {
            col_major[j * n1 + i] = distances[i * n2 + j];
        }
    }

    RMatrix::new_matrix(n1, n2, |r, c| col_major[c * n1 + r])
}

/// Compute Manhattan distance matrix
pub fn manhattan_matrix(x1: &[f64], y1: &[f64], x2: &[f64], y2: &[f64]) -> RMatrix<f64> {
    let n1 = x1.len();
    let n2 = x2.len();

    let distances: Vec<f64> = if n1 * n2 > 10000 {
        (0..n1)
            .into_par_iter()
            .flat_map(|i| {
                (0..n2)
                    .map(|j| (x1[i] - x2[j]).abs() + (y1[i] - y2[j]).abs())
                    .collect::<Vec<_>>()
            })
            .collect()
    } else {
        let mut distances = Vec::with_capacity(n1 * n2);
        for i in 0..n1 {
            for j in 0..n2 {
                distances.push((x1[i] - x2[j]).abs() + (y1[i] - y2[j]).abs());
            }
        }
        distances
    };

    let mut col_major = vec![0.0; n1 * n2];
    for i in 0..n1 {
        for j in 0..n2 {
            col_major[j * n1 + i] = distances[i * n2 + j];
        }
    }

    RMatrix::new_matrix(n1, n2, |r, c| col_major[c * n1 + r])
}

/// Compute pairwise attribute dissimilarity matrix
pub fn attribute_dissimilarity(attrs: &[Vec<f64>], metric: &str) -> Vec<Vec<f64>> {
    let n = attrs.len();
    let mut dissim = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let d = match metric {
                "manhattan" => attrs[i]
                    .iter()
                    .zip(attrs[j].iter())
                    .map(|(a, b)| (a - b).abs())
                    .sum(),
                "euclidean" => {
                    let sum_sq: f64 = attrs[i]
                        .iter()
                        .zip(attrs[j].iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    sum_sq.sqrt()
                }
                _ => {
                    // Default to euclidean
                    let sum_sq: f64 = attrs[i]
                        .iter()
                        .zip(attrs[j].iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum();
                    sum_sq.sqrt()
                }
            };
            dissim[i][j] = d;
            dissim[j][i] = d;
        }
    }

    dissim
}
