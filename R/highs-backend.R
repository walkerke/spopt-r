# Internal solver functions using the CRAN highs package.
# These replicate the MIP models from the Rust HiGHS backend.
# Each function takes the same arguments and returns the same list
# structure as the corresponding rust_* wrapper.
#
# All constraint matrices are built with vectorized triplet construction
# (no nested R loops). The pattern for assignment models (p-median,
# p-center, cflp) uses shared index-generation helpers.

# ---------------------------------------------------------------------------
# Status check: stop with informative error if solver did not find optimal
# ---------------------------------------------------------------------------
.highs_is_optimal <- function(res) {
  # HiGHS status 7 = kOptimal. Use integer for stability across versions.
  res$status == 7L
}

.check_highs_status <- function(res, solver_name) {
  if (!.highs_is_optimal(res)) {
    stop(sprintf("%s solver returned non-optimal status: %s", solver_name,
                 res$status_message), call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# LSCP: Location Set Covering Problem
# Minimize number of facilities to cover all demand within service_radius.
# ---------------------------------------------------------------------------
.solve_lscp <- function(cost_matrix, service_radius) {
  n_demand <- nrow(cost_matrix)
  n_fac <- ncol(cost_matrix)

  L <- rep(1.0, n_fac)
  lower <- rep(0, n_fac)
  upper <- rep(1, n_fac)
  types <- rep("I", n_fac)

  coverage <- cost_matrix <= service_radius
  coverable <- which(rowSums(coverage) > 0)
  uncoverable <- n_demand - length(coverable)

  if (length(coverable) == 0) {
    return(list(
      selected = integer(0), n_selected = 0L, objective = 0,
      covered_demand = 0L, total_demand = n_demand,
      coverage_pct = 0, status = "Optimal", uncoverable_demand = uncoverable
    ))
  }

  A <- coverage[coverable, , drop = FALSE]
  storage.mode(A) <- "double"
  A <- Matrix::Matrix(A, sparse = TRUE)
  lhs <- rep(1, nrow(A))
  rhs <- rep(Inf, nrow(A))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "LSCP")

  sol <- res$primal_solution
  selected <- which(sol > 0.5)

  covered <- rowSums(cost_matrix[, selected, drop = FALSE] <= service_radius) > 0
  covered_demand <- sum(covered)
  coverage_pct <- (covered_demand / n_demand) * 100

  list(
    selected = as.integer(selected),
    n_selected = length(selected),
    objective = res$info$objective_function_value,
    covered_demand = as.integer(covered_demand),
    total_demand = as.integer(n_demand),
    coverage_pct = coverage_pct,
    status = res$status_message,
    uncoverable_demand = as.integer(uncoverable)
  )
}

# ---------------------------------------------------------------------------
# P-Median: minimize total weighted distance with exactly p facilities.
# ---------------------------------------------------------------------------
.solve_p_median <- function(cost_matrix, weights, n_facilities, fixed_facilities) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)
  p <- n_facilities

  if (is.null(fixed_facilities)) fixed_facilities <- integer(0)

  # Variable layout: y[1..n_f], x[n_f+1 .. n_f+n_d*n_f] (row-major x[i,j])
  n_vars <- n_f + n_d * n_f

  # Objective: 0 for y, weights[i]*cost[i,j] for x (vectorized)
  L <- c(rep(0, n_f), as.vector(t(cost_matrix * weights)))

  # Bounds + fixed facilities
  lower <- rep(0, n_vars)
  upper <- rep(1, n_vars)
  lower[fixed_facilities] <- 1

  types <- c(rep("I", n_f), rep("C", n_d * n_f))

  # Sparse constraint matrix (fully vectorized triplet construction)
  A <- .build_assignment_constraints(n_d, n_f)

  lhs <- c(p, rep(1, n_d), rep(-Inf, n_d * n_f))
  rhs <- c(p, rep(1, n_d), rep(0, n_d * n_f))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "P-Median")

  sol <- res$primal_solution
  selected <- which(sol[seq_len(n_f)] > 0.5)

  # Extract assignments: reshape x portion into matrix, take column of max per row
  x_mat <- matrix(sol[(n_f + 1L):n_vars], nrow = n_d, ncol = n_f, byrow = TRUE)
  assignments <- max.col(x_mat, ties.method = "first")

  # Mean weighted distance (vectorized)
  assigned_costs <- cost_matrix[cbind(seq_len(n_d), assignments)]
  mean_distance <- sum(weights * assigned_costs) / sum(weights)

  list(
    selected = as.integer(selected),
    assignments = as.integer(assignments),
    n_selected = length(selected),
    objective = res$info$objective_function_value,
    mean_distance = mean_distance
  )
}

# ---------------------------------------------------------------------------
# Shared helper: build the standard assignment constraint matrix.
#
# Variable layout: y[1..n_f], x[n_f+1 .. n_f+n_d*n_f] (row-major)
# Constraints:
#   Row 1:              sum_j y[j] = p
#   Rows 2..n_d+1:      sum_j x[i,j] = 1 for each i
#   Rows n_d+2..end:     x[i,j] - y[j] <= 0 for each i,j
#
# Returns a sparse Matrix of dimensions (1 + n_d + n_d*n_f) x (n_f + n_d*n_f)
# ---------------------------------------------------------------------------
.build_assignment_constraints <- function(n_d, n_f) {
  n_vars <- n_f + n_d * n_f
  n_con <- 1L + n_d + n_d * n_f

  # Constraint 1: sum_j y[j] = p  (n_f nonzeros)
  c1_i <- rep(1L, n_f)
  c1_j <- seq_len(n_f)
  c1_x <- rep(1, n_f)

  # Constraint 2: sum_j x[i,j] = 1  (n_d * n_f nonzeros)
  # Row index: 1 + i, repeated n_f times per demand
  c2_i <- rep(1L + seq_len(n_d), each = n_f)
  # Column index: n_f + (i-1)*n_f + j
  c2_j <- n_f + seq_len(n_d * n_f)
  c2_x <- rep(1, n_d * n_f)

  # Constraint 3: x[i,j] - y[j] <= 0  (2 * n_d * n_f nonzeros)
  # Row index for constraint (i,j): 1 + n_d + (i-1)*n_f + j
  con3_rows <- 1L + n_d + seq_len(n_d * n_f)

  # x[i,j] coefficient (+1): same column as c2_j
  c3a_i <- con3_rows
  c3a_j <- n_f + seq_len(n_d * n_f)
  c3a_x <- rep(1, n_d * n_f)

  # -y[j] coefficient (-1): column j cycles 1..n_f for each demand
  c3b_i <- con3_rows
  c3b_j <- rep(seq_len(n_f), times = n_d)
  c3b_x <- rep(-1, n_d * n_f)

  Matrix::sparseMatrix(
    i = c(c1_i, c2_i, c3a_i, c3b_i),
    j = c(c1_j, c2_j, c3a_j, c3b_j),
    x = c(c1_x, c2_x, c3a_x, c3b_x),
    dims = c(n_con, n_vars)
  )
}

# ---------------------------------------------------------------------------
# MCLP: Maximum Coverage Location Problem
# Maximize weighted demand covered with exactly p facilities.
# ---------------------------------------------------------------------------
.solve_mclp <- function(cost_matrix, weights, service_radius, n_facilities,
                        fixed_facilities) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)
  p <- n_facilities

  if (is.null(fixed_facilities)) fixed_facilities <- integer(0)

  # Variables: y[1..n_f] binary, z[n_f+1..n_f+n_d] binary
  n_vars <- n_f + n_d

  # Objective: maximize sum(weights[i] * z[i])
  L <- c(rep(0, n_f), weights)

  lower <- rep(0, n_vars)
  upper <- rep(1, n_vars)
  lower[fixed_facilities] <- 1
  types <- rep("I", n_vars)

  # Coverage: a[i,j] = 1 if cost_matrix[i,j] <= service_radius
  coverage <- cost_matrix <= service_radius

  # Constraints:
  # 1: sum_j y[j] = p                                   (1 row)
  # 2: z[i] - sum_j(a[i,j]*y[j]) <= 0 for each i       (n_d rows)
  n_con <- 1L + n_d

  # Constraint 1: facility count
  c1_i <- rep(1L, n_f)
  c1_j <- seq_len(n_f)
  c1_x <- rep(1, n_f)

  # Constraint 2: z[i] - sum covering y[j] <= 0 (vectorized via coverage matrix)
  cov_which <- which(coverage, arr.ind = TRUE)  # rows = demand, cols = facility
  # z[i] coefficient: +1 for each demand
  c2a_i <- 1L + seq_len(n_d)
  c2a_j <- n_f + seq_len(n_d)
  c2a_x <- rep(1, n_d)

  # -a[i,j]*y[j] coefficients: one entry per TRUE in coverage matrix
  c2b_i <- 1L + cov_which[, 1]       # constraint row = 1 + demand index
  c2b_j <- cov_which[, 2]            # variable column = facility index
  c2b_x <- rep(-1, nrow(cov_which))

  A <- Matrix::sparseMatrix(
    i = c(c1_i, c2a_i, c2b_i),
    j = c(c1_j, c2a_j, c2b_j),
    x = c(c1_x, c2a_x, c2b_x),
    dims = c(n_con, n_vars)
  )

  lhs <- c(p, rep(-Inf, n_d))
  rhs <- c(p, rep(0, n_d))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = TRUE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "MCLP")

  sol <- res$primal_solution
  selected <- which(sol[seq_len(n_f)] > 0.5)

  # Coverage statistics
  z_sol <- sol[(n_f + 1L):n_vars]
  covered_weight <- sum(weights[z_sol > 0.5])
  total_weight <- sum(weights)
  coverage_pct <- (covered_weight / total_weight) * 100

  list(
    selected = as.integer(selected),
    n_selected = length(selected),
    objective = res$info$objective_function_value,
    covered_weight = covered_weight,
    total_weight = total_weight,
    coverage_pct = coverage_pct
  )
}

# ---------------------------------------------------------------------------
# P-Center: minimize maximum distance from demand to nearest facility.
# Dispatches to binary_search or MIP method.
# ---------------------------------------------------------------------------
.solve_p_center <- function(cost_matrix, n_facilities, method, fixed_facilities) {
  if (method == "mip") {
    .solve_p_center_mip(cost_matrix, n_facilities, fixed_facilities)
  } else {
    .solve_p_center_bs(cost_matrix, n_facilities, fixed_facilities)
  }
}

# Binary search: find minimum feasible distance threshold
.solve_p_center_bs <- function(cost_matrix, n_facilities, fixed_facilities) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)
  p <- n_facilities

  if (is.null(fixed_facilities)) fixed_facilities <- integer(0)

  # Unique sorted distances
  distances <- sort(unique(as.vector(cost_matrix[is.finite(cost_matrix)])))

  if (length(distances) == 0) {
    return(list(
      selected = integer(0), assignments = rep(0L, n_d),
      n_selected = 0L, objective = NA_real_, max_distance = NA_real_
    ))
  }

  # Binary search for minimum feasible threshold
  lo <- 1L
  hi <- length(distances)
  best_selected <- NULL
  best_threshold <- NA_real_

  while (lo <= hi) {
    mid <- lo + (hi - lo) %/% 2L
    threshold <- distances[mid]
    selected <- .can_cover_with_p(cost_matrix, threshold, p, fixed_facilities)

    if (!is.null(selected)) {
      best_selected <- selected
      best_threshold <- threshold
      if (mid == 1L) break
      hi <- mid - 1L
    } else {
      lo <- mid + 1L
    }
  }

  if (is.null(best_selected)) {
    return(list(
      selected = integer(0), assignments = rep(0L, n_d),
      n_selected = 0L, objective = NA_real_, max_distance = NA_real_
    ))
  }

  # Assignments: nearest selected facility per demand (vectorized)
  sel_dists <- cost_matrix[, best_selected, drop = FALSE]
  assignments <- best_selected[max.col(-sel_dists, ties.method = "first")]

  list(
    selected = as.integer(best_selected),
    assignments = as.integer(assignments),
    n_selected = length(best_selected),
    objective = best_threshold,
    max_distance = best_threshold
  )
}

# Inner feasibility check: can all demand be covered with exactly p facilities?
.can_cover_with_p <- function(cost_matrix, threshold, p, fixed_facilities) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)

  coverage <- cost_matrix <= threshold

  # If any demand has no covering facility, infeasible
  if (any(rowSums(coverage) == 0)) return(NULL)

  # ILP: y[j] binary, sum y = p, each demand covered
  lower <- rep(0, n_f)
  upper <- rep(1, n_f)
  lower[fixed_facilities] <- 1

  # Constraint 1: sum y = p
  # Constraint 2+: coverage for each demand
  A_cov <- coverage
  storage.mode(A_cov) <- "double"

  # Build full A: row 1 = facility count, rows 2+ = coverage
  c1_row <- Matrix::sparseMatrix(
    i = rep(1L, n_f), j = seq_len(n_f), x = rep(1, n_f),
    dims = c(1L, n_f)
  )
  A <- rbind(c1_row, Matrix::Matrix(A_cov, sparse = TRUE))

  lhs <- c(p, rep(1, n_d))
  rhs <- c(p, rep(Inf, n_d))

  res <- highs::highs_solve(
    L = rep(0, n_f), lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = rep("I", n_f), maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  if (!.highs_is_optimal(res)) return(NULL)

  selected <- which(res$primal_solution > 0.5)
  if (length(selected) == p) selected else NULL
}

# MIP: direct minimax formulation
.solve_p_center_mip <- function(cost_matrix, n_facilities, fixed_facilities) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)
  p <- n_facilities

  if (is.null(fixed_facilities)) fixed_facilities <- integer(0)

  # Variables: W (col 1), y[2..n_f+1], x[n_f+2..n_f+1+n_d*n_f]
  n_vars <- 1L + n_f + n_d * n_f
  max_dist <- max(cost_matrix)

  # Objective: minimize W
  L <- c(1, rep(0, n_f + n_d * n_f))

  lower <- rep(0, n_vars)
  upper <- c(max_dist * 2, rep(1, n_f + n_d * n_f))
  lower[1L + fixed_facilities] <- 1  # y offset by 1 for W

  types <- c("C", rep("I", n_f), rep("C", n_d * n_f))

  # Reuse assignment constraints for y and x (offset columns by 1 for W)
  A_assign <- .build_assignment_constraints(n_d, n_f)
  # Extract triplets and shift columns by 1 to account for W in column 1
  trip <- Matrix::summary(A_assign)
  A_assign_shifted <- Matrix::sparseMatrix(
    i = trip$i, j = trip$j + 1L, x = trip$x,
    dims = c(nrow(A_assign), n_vars)
  )

  # Constraint 4: sum_j d[i,j]*x[i,j] - W <= 0 for each demand i (n_d rows)
  # W coefficient: -1 in column 1
  c4_w_i <- seq_len(n_d)
  c4_w_j <- rep(1L, n_d)
  c4_w_x <- rep(-1, n_d)

  # d[i,j]*x[i,j] coefficients
  x_offset <- 1L + n_f  # x variables start at column n_f+2
  c4_x_i <- rep(seq_len(n_d), each = n_f)
  c4_x_j <- x_offset + seq_len(n_d * n_f)
  c4_x_x <- as.vector(t(cost_matrix))

  A_minimax <- Matrix::sparseMatrix(
    i = c(c4_w_i, c4_x_i),
    j = c(c4_w_j, c4_x_j),
    x = c(c4_w_x, c4_x_x),
    dims = c(n_d, n_vars)
  )

  A <- rbind(A_assign_shifted, A_minimax)

  n_assign_con <- nrow(A_assign)
  lhs <- c(p, rep(1, n_d), rep(-Inf, n_d * n_f), rep(-Inf, n_d))
  rhs <- c(p, rep(1, n_d), rep(0, n_d * n_f), rep(0, n_d))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "P-Center (MIP)")

  sol <- res$primal_solution
  w_val <- sol[1]
  selected <- which(sol[2:(n_f + 1)] > 0.5)

  # Assignments from x variables
  x_mat <- matrix(sol[(n_f + 2L):n_vars], nrow = n_d, ncol = n_f, byrow = TRUE)
  assignments <- max.col(x_mat, ties.method = "first")

  list(
    selected = as.integer(selected),
    assignments = as.integer(assignments),
    n_selected = length(selected),
    objective = w_val,
    max_distance = w_val
  )
}

# ---------------------------------------------------------------------------
# P-Dispersion: maximize minimum inter-facility distance.
# ---------------------------------------------------------------------------
.solve_p_dispersion <- function(distance_matrix, n_facilities) {
  n <- nrow(distance_matrix)
  p <- n_facilities
  max_dist <- max(distance_matrix)
  big_m <- max_dist + 1

  # Variables: D (col 1, continuous), y[2..n+1] (binary)
  n_vars <- 1L + n
  n_pairs <- n * (n - 1L) / 2L

  L <- c(1, rep(0, n))
  lower <- c(0, rep(0, n))
  upper <- c(max_dist, rep(1, n))
  types <- c("C", rep("I", n))

  # Constraints:
  # 1: sum_i y[i] = p                              (1 row)
  # 2: D + M*y[i] + M*y[j] <= d[i,j] + 2M         (n_pairs rows)
  n_con <- 1L + n_pairs

  # Constraint 1: facility count (vectorized)
  c1_i <- rep(1L, n)
  c1_j <- 1L + seq_len(n)
  c1_x <- rep(1, n)

  # Constraint 2: big-M constraints for all i < j (vectorized)
  pairs <- which(upper.tri(distance_matrix), arr.ind = TRUE)  # rows: i, cols: j
  pair_rows <- seq_len(nrow(pairs)) + 1L  # constraint rows 2..n_pairs+1

  # D coefficient (+1)
  c2a_i <- pair_rows
  c2a_j <- rep(1L, nrow(pairs))
  c2a_x <- rep(1, nrow(pairs))

  # M*y[i] coefficient
  c2b_i <- pair_rows
  c2b_j <- 1L + pairs[, 1]  # y variable for facility i
  c2b_x <- rep(big_m, nrow(pairs))

  # M*y[j] coefficient
  c2c_i <- pair_rows
  c2c_j <- 1L + pairs[, 2]  # y variable for facility j
  c2c_x <- rep(big_m, nrow(pairs))

  A <- Matrix::sparseMatrix(
    i = c(c1_i, c2a_i, c2b_i, c2c_i),
    j = c(c1_j, c2a_j, c2b_j, c2c_j),
    x = c(c1_x, c2a_x, c2b_x, c2c_x),
    dims = c(n_con, n_vars)
  )

  # RHS: constraint 1 = p, constraint 2 = d[i,j] + 2*M
  pair_rhs <- distance_matrix[pairs] + 2 * big_m
  lhs <- c(p, rep(-Inf, n_pairs))
  rhs <- c(p, pair_rhs)

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = TRUE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "P-Dispersion")

  sol <- res$primal_solution
  d_val <- sol[1]
  selected <- which(sol[2:(n + 1)] > 0.5)

  list(
    selected = as.integer(selected),
    n_selected = length(selected),
    objective = d_val,
    min_distance = d_val
  )
}

# ---------------------------------------------------------------------------
# CFLP: Capacitated Facility Location Problem.
# Minimize transport + facility costs with capacity constraints.
# ---------------------------------------------------------------------------
.solve_cflp <- function(cost_matrix, weights, capacities, n_facilities,
                        facility_costs, fixed_facilities = NULL) {
  n_d <- nrow(cost_matrix)
  n_f <- ncol(cost_matrix)

  if (is.null(facility_costs)) facility_costs <- rep(0, n_f)

  # Variable layout: y[1..n_f], x[n_f+1..n_f+n_d*n_f] (row-major)
  n_vars <- n_f + n_d * n_f

  # Objective: facility_costs[j]*y[j] + weights[i]*cost[i,j]*x[i,j]
  L <- c(facility_costs, as.vector(t(cost_matrix * weights)))

  lower <- rep(0, n_vars)
  upper <- rep(1, n_vars)
  lower[fixed_facilities] <- 1
  types <- c(rep("I", n_f), rep("C", n_d * n_f))

  # Build constraints using assignment helper, then append capacity rows
  A_assign <- .build_assignment_constraints(n_d, n_f)
  n_assign_con <- nrow(A_assign)

  # If n_facilities == 0, remove the first row (facility count constraint)
  if (n_facilities == 0L) {
    A_assign <- A_assign[-1, , drop = FALSE]
    n_assign_con <- nrow(A_assign)
    assign_lhs <- c(rep(1, n_d), rep(-Inf, n_d * n_f))
    assign_rhs <- c(rep(1, n_d), rep(0, n_d * n_f))
  } else {
    assign_lhs <- c(n_facilities, rep(1, n_d), rep(-Inf, n_d * n_f))
    assign_rhs <- c(n_facilities, rep(1, n_d), rep(0, n_d * n_f))
  }

  # Capacity constraints: sum_i w[i]*x[i,j] - cap[j]*y[j] <= 0 for each j
  # (n_f rows, each with n_d + 1 nonzeros)
  cap_rows <- seq_len(n_f)

  # -cap[j]*y[j] coefficient
  cap_y_i <- cap_rows
  cap_y_j <- seq_len(n_f)
  cap_y_x <- -capacities

  # w[i]*x[i,j] coefficients: for each facility j, n_d entries
  # x[i,j] is at column n_f + (i-1)*n_f + j
  cap_x_i <- rep(cap_rows, each = n_d)
  cap_x_j <- as.vector(sapply(seq_len(n_f), function(j) n_f + (seq_len(n_d) - 1L) * n_f + j))
  cap_x_x <- rep(weights, times = n_f)

  A_cap <- Matrix::sparseMatrix(
    i = c(cap_y_i, cap_x_i),
    j = c(cap_y_j, cap_x_j),
    x = c(cap_y_x, cap_x_x),
    dims = c(n_f, n_vars)
  )

  A <- rbind(A_assign, A_cap)
  lhs <- c(assign_lhs, rep(-Inf, n_f))
  rhs <- c(assign_rhs, rep(0, n_f))

  res <- highs::highs_solve(
    L = L, lower = lower, upper = upper,
    A = A, lhs = lhs, rhs = rhs,
    types = types, maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  .check_highs_status(res, "CFLP")

  sol <- res$primal_solution
  selected <- which(sol[seq_len(n_f)] > 0.5)

  # Allocation matrix (n_d x n_f, row-major)
  x_vec <- sol[(n_f + 1L):n_vars]
  alloc_matrix <- matrix(x_vec, nrow = n_d, ncol = n_f, byrow = TRUE)

  # Primary assignments
  assignments <- max.col(alloc_matrix, ties.method = "first")

  # Split demand: primary allocation < 0.999
  primary_alloc <- alloc_matrix[cbind(seq_len(n_d), assignments)]
  n_split <- sum(primary_alloc < 0.999)

  # Utilizations: sum(w[i]*x[i,j]) / cap[j] for selected facilities
  utilizations <- numeric(n_f)
  for (j in selected) {
    utilizations[j] <- sum(weights * alloc_matrix[, j]) / capacities[j]
  }

  # Mean weighted distance (from full allocation, not just primary)
  total_weighted_dist <- sum(weights * rowSums(alloc_matrix * cost_matrix))
  mean_distance <- total_weighted_dist / sum(weights)

  list(
    selected = as.integer(selected),
    assignments = as.integer(assignments),
    allocation_matrix = as.vector(t(alloc_matrix)),  # flattened row-major
    n_selected = length(selected),
    n_split_demand = as.integer(n_split),
    objective = res$info$objective_function_value,
    mean_distance = mean_distance,
    utilizations = utilizations,
    status = res$status_message
  )
}
