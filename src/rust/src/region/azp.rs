//! AZP (Automatic Zoning Procedure)
//!
//! Implements the AZP algorithm with basic, tabu, and simulated annealing variants.
//! Based on: Openshaw (1977), Openshaw & Rao (1995)
//!
//! The search itself lives in `local_search`; this module builds the initial
//! partition (random growth or a user-supplied warm start), dispatches to the
//! requested variant, and polishes the result with a greedy pass so that every
//! returned partition is a local optimum under single-area moves.

use extendr_api::prelude::*;
use rand_chacha::ChaCha8Rng;
use rand::SeedableRng;

use super::local_search::{anneal, greedy, grow_initial, tabu, Adj, Partition};

/// Solve AZP regionalization
pub fn solve(
    attrs: RMatrix<f64>,
    n_regions: usize,
    adj_i: &[i32],
    adj_j: &[i32],
    method: &str,
    max_iterations: usize,
    tabu_length: usize,
    cooling_rate: f64,
    initial_temperature: f64,
    initial_regions: Option<Vec<i32>>,
    seed: Option<u64>,
) -> List {
    let n = attrs.nrows();
    let n_attrs = attrs.ncols();

    if n == 0 || n_regions == 0 || n_regions > n {
        return list!(
            labels = Vec::<i32>::new(),
            n_regions = 0i32,
            objective = 0.0f64,
            initial_objective = 0.0f64,
            iterations = 0i32
        );
    }

    let attr_vecs: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n_attrs).map(|j| attrs[[i, j]]).collect())
        .collect();

    let adj = Adj::from_indices(adj_i, adj_j, n);
    let mut rng = ChaCha8Rng::seed_from_u64(seed.unwrap_or(42));

    // Initial partition: 0-based labels, validated on the R side.
    let labels0 = match initial_regions {
        Some(l) if l.len() == n => l,
        _ => grow_initial(n, n_regions, &adj, &mut rng),
    };
    let mut part = Partition::new(&labels0, n_regions, &attr_vecs, None);
    let initial_objective = part.objective(&attr_vecs);

    let iterations: usize = match method {
        "tabu" => {
            let conv = (10 * tabu_length).max(100);
            part = tabu(part, &attr_vecs, &adj, None, max_iterations, tabu_length, conv);
            greedy(&mut part, &attr_vecs, &adj, None, 1000, &mut rng);
            max_iterations
        }
        "sa" => {
            part = anneal(
                part,
                &attr_vecs,
                &adj,
                None,
                max_iterations,
                cooling_rate,
                initial_temperature,
                0,
                &mut rng,
            );
            greedy(&mut part, &attr_vecs, &adj, None, 1000, &mut rng);
            max_iterations
        }
        _ => greedy(&mut part, &attr_vecs, &adj, None, max_iterations, &mut rng),
    };

    let objective = part.objective(&attr_vecs);
    let labels: Vec<i32> = part.labels.iter().map(|&l| l + 1).collect();

    list!(
        labels = labels,
        n_regions = part.k() as i32,
        objective = objective,
        initial_objective = initial_objective,
        iterations = iterations as i32
    )
}
