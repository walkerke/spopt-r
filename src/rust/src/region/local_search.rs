//! Shared local-search engine for contiguity-constrained regionalization.
//!
//! Implements the Automatic Zoning Procedure family of moves (Openshaw 1977;
//! Openshaw & Rao 1995): an area on the border between two regions is moved
//! from its donor region to the recipient when that is feasible, where
//! feasible means the donor keeps at least one member, stays connected, and
//! (optionally) stays above a floor on a constraint variable. Three searches
//! are provided over that move set:
//!
//! * `greedy`: Openshaw's original AZP. Regions are visited in random order;
//!   for each region, border areas are pulled in while any improving move
//!   exists. One iteration is a full pass over all regions. Converges to a
//!   local optimum.
//! * `tabu`: steepest descent / mildest ascent with a short tabu list on
//!   reversed moves and an aspiration rule.
//! * `anneal`: simulated annealing where one iteration is a sweep over the
//!   current border in random order, followed by one cooling step. The
//!   start temperature is calibrated from a sample of uphill moves so that
//!   a typical uphill move is accepted with probability 0.8.
//!
//! Objective deltas are computed in O(p) from running attribute sums
//! (the exact update of a within-region sum of squared deviations when one
//! point is added or removed), and the connectivity check (a depth-first
//! search within the donor region) runs only for moves that would otherwise
//! be accepted. That is what makes thousands of sweeps affordable.
//!
//! Used by `azp` (fixed number of regions) and `maxp` (regions with a
//! floor constraint, number of regions fixed by the construction phase).

use rand::prelude::*;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha8Rng;
use std::collections::VecDeque;

/// Adjacency list built from 0-based edge index vectors.
pub struct Adj {
    pub neighbors: Vec<Vec<usize>>,
}

impl Adj {
    pub fn from_indices(adj_i: &[i32], adj_j: &[i32], n: usize) -> Self {
        let mut neighbors = vec![Vec::new(); n];
        for (&i, &j) in adj_i.iter().zip(adj_j.iter()) {
            let (i, j) = (i as usize, j as usize);
            if i == j || i >= n || j >= n {
                continue;
            }
            if !neighbors[i].contains(&j) {
                neighbors[i].push(j);
            }
            if !neighbors[j].contains(&i) {
                neighbors[j].push(i);
            }
        }
        Self { neighbors }
    }

    #[inline]
    pub fn get(&self, i: usize) -> &[usize] {
        &self.neighbors[i]
    }

}

/// Optional floor constraint: every region must keep `sum(var) >= floor`.
#[derive(Clone, Copy)]
pub struct Floor<'a> {
    pub var: &'a [f64],
    pub min: f64,
}

/// Reusable search buffers: stamp arrays for visited marks, per-search frontiers.
pub struct Scratch {
    mark: Vec<u32>,
    stamp: u32,
    /// Which search (index into `targets`) first reached each marked node.
    comp: Vec<u32>,
    targets: Vec<usize>,
    frontiers: Vec<VecDeque<usize>>,
    uf_parent: Vec<usize>,
    uf_active: Vec<usize>,
    /// Per-region stamp used to deduplicate an area's neighbouring regions.
    region_mark: Vec<u32>,
    region_stamp: u32,
    /// Reused buffer for the current border moves.
    pub moves: Vec<(usize, usize, usize)>,
}

impl Scratch {
    pub fn new(n: usize, k: usize) -> Self {
        Self {
            mark: vec![0; n],
            stamp: 0,
            comp: vec![0; n],
            targets: Vec::with_capacity(16),
            frontiers: Vec::new(),
            uf_parent: Vec::with_capacity(16),
            uf_active: Vec::with_capacity(16),
            region_mark: vec![0; k],
            region_stamp: 0,
            moves: Vec::new(),
        }
    }

    fn next_stamp(&mut self) {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.mark.iter_mut().for_each(|m| *m = 0);
            self.stamp = 1;
        }
    }
}

/// A partition of areas into regions with running sums for O(p) deltas.
#[derive(Clone)]
pub struct Partition {
    /// Region id for each area (always >= 0 here).
    pub labels: Vec<i32>,
    /// Areas in each region.
    pub members: Vec<Vec<usize>>,
    /// Position of each area within its region's member list.
    pos: Vec<usize>,
    /// Attribute sums per region.
    pub sums: Vec<Vec<f64>>,
    /// Sum of the floor variable per region (zeros when no floor).
    pub floor_sums: Vec<f64>,
}

impl Partition {
    /// Build from 0-based labels in `0..k`.
    pub fn new(labels: &[i32], k: usize, attrs: &[Vec<f64>], floor: Option<Floor>) -> Self {
        let p = attrs.first().map(|a| a.len()).unwrap_or(0);
        let mut members = vec![Vec::new(); k];
        let mut pos = vec![0usize; labels.len()];
        let mut sums = vec![vec![0.0; p]; k];
        let mut floor_sums = vec![0.0; k];
        for (i, &l) in labels.iter().enumerate() {
            let r = l as usize;
            pos[i] = members[r].len();
            members[r].push(i);
            for (c, v) in attrs[i].iter().enumerate() {
                sums[r][c] += v;
            }
            if let Some(f) = floor {
                floor_sums[r] += f.var[i];
            }
        }
        Self {
            labels: labels.to_vec(),
            members,
            pos,
            sums,
            floor_sums,
        }
    }

    #[inline]
    pub fn k(&self) -> usize {
        self.members.len()
    }

    /// Total within-region sum of squared deviations from region centroids.
    pub fn objective(&self, attrs: &[Vec<f64>]) -> f64 {
        let mut total = 0.0;
        for (r, members) in self.members.iter().enumerate() {
            if members.is_empty() {
                continue;
            }
            let n = members.len() as f64;
            for &m in members {
                for (c, v) in attrs[m].iter().enumerate() {
                    let d = v - self.sums[r][c] / n;
                    total += d * d;
                }
            }
        }
        total
    }

    /// Exact change in the objective if `area` moves from `from` to `to`.
    ///
    /// Adding x to a region of size n with mean mu raises its SSD by
    /// n/(n+1) * |x - mu|^2; removing x from a region of size n (x included)
    /// lowers its SSD by n/(n-1) * |x - mu|^2.
    #[inline]
    pub fn move_delta(&self, attrs: &[Vec<f64>], area: usize, from: usize, to: usize) -> f64 {
        let x = &attrs[area];
        let n_from = self.members[from].len() as f64;
        let n_to = self.members[to].len() as f64;
        let mut d_to = 0.0;
        let mut d_from = 0.0;
        for (c, &xv) in x.iter().enumerate() {
            let a = xv - self.sums[to][c] / n_to;
            d_to += a * a;
            let b = xv - self.sums[from][c] / n_from;
            d_from += b * b;
        }
        let gain = n_to / (n_to + 1.0) * d_to;
        let loss = if n_from > 1.0 {
            n_from / (n_from - 1.0) * d_from
        } else {
            0.0
        };
        gain - loss
    }

    /// Apply a move. The caller has already checked feasibility.
    pub fn execute(&mut self, attrs: &[Vec<f64>], floor: Option<Floor>, area: usize, to: usize) {
        let from = self.labels[area] as usize;
        let p = self.pos[area];
        let last = *self.members[from].last().expect("non-empty donor region");
        self.members[from].swap_remove(p);
        if last != area {
            self.pos[last] = p;
        }
        self.pos[area] = self.members[to].len();
        self.members[to].push(area);
        for (c, v) in attrs[area].iter().enumerate() {
            self.sums[from][c] -= v;
            self.sums[to][c] += v;
        }
        if let Some(f) = floor {
            self.floor_sums[from] -= f.var[area];
            self.floor_sums[to] += f.var[area];
        }
        self.labels[area] = to as i32;
    }

    /// Can `area` leave its region without emptying it, disconnecting it,
    /// or taking it below the floor?
    pub fn can_leave(&self, area: usize, adj: &Adj, floor: Option<Floor>, scratch: &mut Scratch) -> bool {
        let r = self.labels[area] as usize;
        if self.members[r].len() <= 1 {
            return false;
        }
        if let Some(f) = floor {
            if self.floor_sums[r] - f.var[area] < f.min {
                return false;
            }
        }
        !self.is_articulation(area, adj, scratch)
    }

    /// Would removing `area` disconnect its region?
    ///
    /// Removing a node disconnects a connected region exactly when its
    /// in-region neighbours stop being mutually reachable without it. A node
    /// with at most one in-region neighbour is never an articulation point.
    /// Otherwise a breadth-first search grows from every in-region neighbour
    /// at once, merging searches when they touch. It stops with "connected" as
    /// soon as all searches have merged, and with "disconnected" as soon as
    /// one merged group has nothing left to explore, so the cost is bounded by
    /// the smallest piece that would be cut off rather than by the region.
    /// Visits are marked in a stamp array; no per-call allocation.
    fn is_articulation(&self, area: usize, adj: &Adj, scratch: &mut Scratch) -> bool {
        let r = self.labels[area];
        if self.members[r as usize].len() <= 2 {
            return false;
        }
        scratch.targets.clear();
        for &nb in adj.get(area) {
            if self.labels[nb] == r {
                scratch.targets.push(nb);
            }
        }
        let t = scratch.targets.len();
        if t <= 1 {
            return false;
        }

        scratch.next_stamp();
        let stamp = scratch.stamp;
        scratch.mark[area] = stamp;
        while scratch.frontiers.len() < t {
            scratch.frontiers.push(VecDeque::new());
        }
        // Union-find over the t searches, with the number of searches in each
        // group that still have work.
        let parent = &mut scratch.uf_parent;
        let active = &mut scratch.uf_active;
        parent.clear();
        active.clear();
        for i in 0..t {
            parent.push(i);
            active.push(1);
            let sdx = scratch.targets[i];
            scratch.mark[sdx] = stamp;
            scratch.comp[sdx] = i as u32;
            scratch.frontiers[i].clear();
            scratch.frontiers[i].push_back(sdx);
        }
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        let mut groups = t;
        loop {
            let mut progressed = false;
            for i in 0..t {
                let cur = match scratch.frontiers[i].pop_front() {
                    Some(c) => c,
                    None => continue,
                };
                progressed = true;
                for &nb in adj.get(cur) {
                    if nb == area || self.labels[nb] != r {
                        continue;
                    }
                    if scratch.mark[nb] != stamp {
                        scratch.mark[nb] = stamp;
                        scratch.comp[nb] = i as u32;
                        scratch.frontiers[i].push_back(nb);
                    } else {
                        let j = scratch.comp[nb] as usize;
                        let (ri, rj) = (find(parent, i), find(parent, j));
                        if ri != rj {
                            parent[rj] = ri;
                            active[ri] += active[rj];
                            groups -= 1;
                            if groups == 1 {
                                return false;
                            }
                        }
                    }
                }
                if scratch.frontiers[i].is_empty() {
                    let ri = find(parent, i);
                    active[ri] -= 1;
                    if active[ri] == 0 {
                        // This group is a complete component holding none of
                        // the other neighbours: the region would split.
                        return true;
                    }
                }
            }
            if !progressed {
                return true;
            }
        }
    }

    /// Every (area, from, to) with `area` adjacent to a different region,
    /// in area order then first-seen neighbour order, written into
    /// `scratch.moves` (reused across calls to avoid reallocation).
    pub fn border_moves(&self, adj: &Adj, scratch: &mut Scratch) {
        scratch.moves.clear();
        for (area, &l) in self.labels.iter().enumerate() {
            let from = l as usize;
            scratch.region_stamp = scratch.region_stamp.wrapping_add(1);
            if scratch.region_stamp == 0 {
                scratch.region_mark.iter_mut().for_each(|m| *m = 0);
                scratch.region_stamp = 1;
            }
            let st = scratch.region_stamp;
            scratch.region_mark[from] = st;
            for &nb in adj.get(area) {
                let to = self.labels[nb] as usize;
                if scratch.region_mark[to] != st {
                    scratch.region_mark[to] = st;
                    scratch.moves.push((area, from, to));
                }
            }
        }
    }
}

/// Openshaw's AZP: random region order, pull in improving border areas until
/// none remain for that region, repeat passes until a full pass makes no
/// move or `max_passes` is reached. Returns the number of passes used.
pub fn greedy(
    part: &mut Partition,
    attrs: &[Vec<f64>],
    adj: &Adj,
    floor: Option<Floor>,
    max_passes: usize,
    rng: &mut ChaCha8Rng,
) -> usize {
    let k = part.k();
    let mut scratch = Scratch::new(part.labels.len(), part.k());
    let mut passes = 0;
    while passes < max_passes {
        passes += 1;
        let mut changed = false;
        let mut order: Vec<usize> = (0..k).collect();
        order.shuffle(rng);
        for &r in &order {
            loop {
                // Areas in other regions that touch region r.
                let mut cands: Vec<usize> = Vec::new();
                for &a in &part.members[r] {
                    for &nb in adj.get(a) {
                        if part.labels[nb] as usize != r && !cands.contains(&nb) {
                            cands.push(nb);
                        }
                    }
                }
                cands.shuffle(rng);
                let mut moved = false;
                for a in cands {
                    let from = part.labels[a] as usize;
                    let delta = part.move_delta(attrs, a, from, r);
                    if delta < -1e-12 && part.can_leave(a, adj, floor, &mut scratch) {
                        part.execute(attrs, floor, a, r);
                        moved = true;
                        changed = true;
                        break;
                    }
                }
                if !moved {
                    break;
                }
            }
        }
        if !changed {
            break;
        }
    }
    passes
}

/// Tabu search: at each iteration take the best feasible non-tabu move
/// (or a tabu move that beats the best objective so far). Stops after
/// `max_iter` moves or `conv` consecutive non-improving moves. Returns the
/// best partition visited.
pub fn tabu(
    start: Partition,
    attrs: &[Vec<f64>],
    adj: &Adj,
    floor: Option<Floor>,
    max_iter: usize,
    tabu_length: usize,
    conv: usize,
) -> Partition {
    let mut part = start;
    let mut scratch = Scratch::new(part.labels.len(), part.k());
    let mut current = part.objective(attrs);
    let mut best = current;
    let mut best_labels = part.labels.clone();
    let mut tabu_list: VecDeque<(usize, usize, usize)> = VecDeque::new();
    let mut no_improve = 0usize;

    for _ in 0..max_iter {
        if no_improve >= conv {
            break;
        }
        part.border_moves(adj, &mut scratch);
        let mut moves: Vec<(f64, usize, usize, usize)> = scratch
            .moves
            .iter()
            .map(|&(a, from, to)| (part.move_delta(attrs, a, from, to), a, from, to))
            .collect();
        if moves.is_empty() {
            break;
        }
        moves.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

        let mut chosen: Option<(f64, usize, usize, usize)> = None;
        for &(delta, a, from, to) in &moves {
            let is_tabu = tabu_list.contains(&(a, to, from));
            if is_tabu && current + delta >= best {
                continue;
            }
            if part.can_leave(a, adj, floor, &mut scratch) {
                chosen = Some((delta, a, from, to));
                break;
            }
        }
        let (delta, a, from, to) = match chosen {
            Some(m) => m,
            None => break,
        };
        part.execute(attrs, floor, a, to);
        current += delta;
        tabu_list.push_back((a, from, to));
        while tabu_list.len() > tabu_length {
            tabu_list.pop_front();
        }
        if current < best - 1e-12 {
            best = current;
            best_labels.clone_from(&part.labels);
            no_improve = 0;
        } else {
            no_improve += 1;
        }
    }
    Partition::new(&best_labels, part.k(), attrs, floor)
}

/// Calibrate a start temperature: the mean uphill delta over a sample of
/// border moves, scaled so that it is accepted with probability `p_accept`.
pub fn calibrate_temperature(
    part: &Partition,
    attrs: &[Vec<f64>],
    adj: &Adj,
    p_accept: f64,
    rng: &mut ChaCha8Rng,
) -> f64 {
    let mut scratch = Scratch::new(part.labels.len(), part.k());
    part.border_moves(adj, &mut scratch);
    let mut moves = std::mem::take(&mut scratch.moves);
    moves.shuffle(rng);
    let uphill: Vec<f64> = moves
        .iter()
        .take(500)
        .map(|&(a, from, to)| part.move_delta(attrs, a, from, to))
        .filter(|d| *d > 0.0)
        .collect();
    if uphill.is_empty() {
        return 1.0;
    }
    let mean = uphill.iter().sum::<f64>() / uphill.len() as f64;
    -mean / p_accept.ln()
}

/// Simulated annealing. One iteration is a sweep over the current border in
/// random order followed by one cooling step. Stops after `sweeps` sweeps,
/// when three consecutive sweeps accept nothing, or when the temperature
/// underflows. An optional short tabu list blocks immediate reversals.
/// Returns the best partition visited.
pub fn anneal(
    start: Partition,
    attrs: &[Vec<f64>],
    adj: &Adj,
    floor: Option<Floor>,
    sweeps: usize,
    cooling_rate: f64,
    initial_temperature: f64,
    tabu_length: usize,
    rng: &mut ChaCha8Rng,
) -> Partition {
    let mut part = start;
    let mut scratch = Scratch::new(part.labels.len(), part.k());
    let mut temperature = if initial_temperature > 0.0 {
        initial_temperature
    } else {
        calibrate_temperature(&part, attrs, adj, 0.8, rng)
    };
    let mut current = part.objective(attrs);
    let mut best = current;
    let mut best_labels = part.labels.clone();
    let mut tabu_list: VecDeque<(usize, usize, usize)> = VecDeque::new();
    let mut idle_sweeps = 0usize;

    for _ in 0..sweeps {
        if temperature < 1e-12 {
            break;
        }
        part.border_moves(adj, &mut scratch);
        let mut moves = std::mem::take(&mut scratch.moves);
        moves.shuffle(rng);
        let mut accepted = 0usize;
        for &(a, from, to) in &moves {
            // Earlier moves in this sweep may have moved the area, or moved
            // away the neighbour that made `to` adjacent to it.
            if part.labels[a] as usize != from
                || !adj.get(a).iter().any(|&nb| part.labels[nb] as usize == to)
            {
                continue;
            }
            if tabu_length > 0 && tabu_list.contains(&(a, to, from)) {
                continue;
            }
            let delta = part.move_delta(attrs, a, from, to);
            let accept = delta < 0.0 || rng.gen::<f64>() < (-delta / temperature).exp();
            if accept && part.can_leave(a, adj, floor, &mut scratch) {
                part.execute(attrs, floor, a, to);
                current += delta;
                accepted += 1;
                if tabu_length > 0 {
                    tabu_list.push_back((a, from, to));
                    while tabu_list.len() > tabu_length {
                        tabu_list.pop_front();
                    }
                }
                if current < best - 1e-12 {
                    best = current;
                    best_labels.clone_from(&part.labels);
                }
            }
        }
        scratch.moves = moves;
        temperature *= cooling_rate;
        if accepted == 0 {
            idle_sweeps += 1;
            if idle_sweeps >= 3 {
                break;
            }
        } else {
            idle_sweeps = 0;
        }
    }
    Partition::new(&best_labels, part.k(), attrs, floor)
}

/// Grow `k` contiguous regions from spread-out random seeds by round-robin
/// breadth-first expansion. Assumes the graph is connected.
pub fn grow_initial(n: usize, k: usize, adj: &Adj, rng: &mut ChaCha8Rng) -> Vec<i32> {
    let mut labels = vec![-1i32; n];
    let mut order: Vec<usize> = (0..n).collect();
    order.shuffle(rng);

    let mut seeds: Vec<usize> = Vec::new();
    for &c in &order {
        if seeds.len() >= k {
            break;
        }
        if seeds.is_empty() || seeds.iter().all(|&s| !adj.get(s).contains(&c)) {
            seeds.push(c);
        }
    }
    for &c in &order {
        if seeds.len() >= k {
            break;
        }
        if !seeds.contains(&c) {
            seeds.push(c);
        }
    }
    for (r, &s) in seeds.iter().enumerate() {
        labels[s] = r as i32;
    }

    let mut frontiers: Vec<VecDeque<usize>> = seeds
        .iter()
        .map(|&s| adj.get(s).iter().copied().collect())
        .collect();
    let mut assigned = seeds.len();
    let mut stalled = 0usize;
    while assigned < n && stalled < k {
        stalled = 0;
        for r in 0..k {
            let mut took = false;
            while let Some(a) = frontiers[r].pop_front() {
                if labels[a] >= 0 {
                    continue;
                }
                labels[a] = r as i32;
                assigned += 1;
                for &nb in adj.get(a) {
                    if labels[nb] < 0 {
                        frontiers[r].push_back(nb);
                    }
                }
                took = true;
                break;
            }
            if !took {
                // Refill from any unassigned area touching this region.
                for a in 0..n {
                    if labels[a] < 0 && adj.get(a).iter().any(|&nb| labels[nb] == r as i32) {
                        frontiers[r].push_back(a);
                    }
                }
                if frontiers[r].is_empty() {
                    stalled += 1;
                }
            }
        }
    }
    // Any area still unreached (disconnected graph) joins region 0.
    for l in labels.iter_mut() {
        if *l < 0 {
            *l = 0;
        }
    }
    labels
}
