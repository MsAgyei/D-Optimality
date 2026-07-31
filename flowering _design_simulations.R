# ============================================================
# SIMULATION 1 — Two-point visit design comparison
# True model: logit(p) = alpha + beta * x   (alpha = 0, beta = 1)
# True Day50% occurs at x = 0 (where p = 0.5)
#
# For each grid point (x1, x2) with x2 > x1, computes:
#   1. det(M)   — D-optimality criterion
#   2. M.c      — C-optimality criterion (delta method for Day50%)
#   3. GLM RMSE — Day50% estimation error from the fitted GLM
#   4. Method 1 — quick Day50% estimate on the proportion scale
#   5. Method 2 — quick Day50% estimate on the logit scale
# Each criterion is averaged over n_sim simulations per grid point,
# repeated for N_bunches = 10, 20, 50.
# ============================================================

library(dplyr)

set.seed(42)

# ------------------------------------------------------------
# TRUE MODEL PARAMETERS
# ------------------------------------------------------------
alpha <- 0   # true intercept
beta  <- 1   # true slope

# ------------------------------------------------------------
# FIXED PARAMETERS
# ------------------------------------------------------------
N              <- 200               # flowers per bunch
n_sim          <- 1000              # simulations per grid point
n_bunches_vec  <- c(10, 20, 50)     # inflorescence counts to compare

# ------------------------------------------------------------
# FULL 2D TRIANGULAR GRID (x2 > x1)
# Both visits range from 5% to 95% flowering
# ------------------------------------------------------------
p_grid <- seq(0.05, 0.95, by = 0.005)
x_grid <- log(p_grid / (1 - p_grid))

heat_grid_base <- expand.grid(x1 = x_grid, x2 = x_grid) %>%
  filter(x2 > x1) %>%
  mutate(
    p1_pct = exp(x1) / (1 + exp(x1)) * 100,
    p2_pct = exp(x2) / (1 + exp(x2)) * 100
  )

# ============================================================
# SIMULATION FUNCTION
# Runs n_sim repeated GLM fits at a single grid point (x1, x2)
# and summarises each optimality/estimation criterion
# ============================================================
sim1_fn <- function(x1, x2, N_bunches, N, n_sim) {
  p1 <- exp(x1) / (1 + exp(x1))
  p2 <- exp(x2) / (1 + exp(x2))

  det_M_vec  <- numeric(n_sim)  # D-optimality: det(information matrix)
  Mc_vec     <- numeric(n_sim)  # C-optimality: delta-method variance of Day50%
  c_rmse_vec <- numeric(n_sim)  # GLM Day50% error (true Day50% = 0)
  m1_vec     <- numeric(n_sim)  # Method 1 Day50% error — proportion scale
  m2_vec     <- numeric(n_sim)  # Method 2 Day50% error — logit scale
  p1_obs     <- numeric(n_sim)  # observed proportion at early visit
  p2_obs     <- numeric(n_sim)  # observed proportion at late visit

  for (i in seq_len(n_sim)) {
    y1 <- rbinom(N_bunches, N, p1)  # open flowers at early visit
    y2 <- rbinom(N_bunches, N, p2)  # open flowers at late visit

    p1_obs[i] <- mean(y1 / N)
    p2_obs[i] <- mean(y2 / N)

    df <- data.frame(
      y = c(y1, y2),
      n = N,
      x = rep(c(x1, x2), each = N_bunches)
    )

    g <- glm(cbind(y, n - y) ~ x, data = df, family = binomial())

    a <- g$coef[1]
    b <- g$coef[2]

    # D-optimality: det of information matrix
    M            <- solve(vcov(g))
    det_M_vec[i] <- det(M)

    # C-optimality: delta-method variance of Day50% = -a/b
    c.vec     <- c(-1/b, a/b^2)
    Mc_vec[i] <- as.numeric(t(c.vec) %*% solve(M) %*% c.vec)

    # GLM Day50% error (true Day50% = 0)
    c_rmse_vec[i] <- -a / b

    # Method 1: proportion scale
    # Solve pbar = a + b*x at each visit, then Day50% where p = 0.5
    b_m1      <- (p2_obs[i] - p1_obs[i]) / (x2 - x1)
    a_m1      <- p1_obs[i] - b_m1 * x1
    m1_vec[i] <- (0.5 - a_m1) / b_m1

    # Method 2: logit scale
    # Solve logit(pbar) = a + b*x at each visit, then Day50% where logit(p) = 0
    # Guard against p = 0 or p = 1, which makes logit blow up
    if (p1_obs[i] > 0 & p1_obs[i] < 1 &
        p2_obs[i] > 0 & p2_obs[i] < 1) {
      lp1       <- log(p1_obs[i] / (1 - p1_obs[i]))
      lp2       <- log(p2_obs[i] / (1 - p2_obs[i]))
      b_m2      <- (lp2 - lp1) / (x2 - x1)
      a_m2      <- lp1 - b_m2 * x1
      m2_vec[i] <- -a_m2 / b_m2
    } else {
      m2_vec[i] <- NA_real_
    }
  }

  list(
    det_M   = mean(det_M_vec),
    Mc      = mean(Mc_vec),
    c_rmse  = sqrt(mean(c_rmse_vec^2)),
    m1_rmse = sqrt(mean(m1_vec^2, na.rm = TRUE)),
    m2_rmse = sqrt(mean(m2_vec^2, na.rm = TRUE)),
    p1_obs  = p1_obs,
    p2_obs  = p2_obs
  )
}

# ============================================================
# RUN SIMULATION 1 FOR EACH N_bunches VALUE
# Results are kept in memory (sim1_heat_grid_results /
# sim1_heat_res_results), keyed by N_bunches — nothing written to disk.
# ============================================================
sim1_heat_grid_results <- list()  # one heat_grid per N_bunches
sim1_heat_res_results  <- list()  # raw sim1_fn output per N_bunches

for (N_bunches in n_bunches_vec) {

  cat("Running Simulation 1 | N_bunches =", N_bunches,
      "| n_sim =", n_sim, "| grid points =", nrow(heat_grid_base), "\n")

  heat_grid <- heat_grid_base
  heat_res  <- vector("list", nrow(heat_grid))

  for (i in seq_len(nrow(heat_grid))) {
    if (i %% 50 == 0) cat("  point", i, "of", nrow(heat_grid), "\n")
    heat_res[[i]] <- sim1_fn(heat_grid$x1[i], heat_grid$x2[i],
                              N_bunches, N, n_sim)
  }

  # Extract results into the grid
  heat_grid$det_M   <- sapply(heat_res, function(r) r$det_M)
  heat_grid$Mc      <- sapply(heat_res, function(r) r$Mc)
  heat_grid$c_rmse  <- sapply(heat_res, function(r) r$c_rmse)
  heat_grid$m1_rmse <- sapply(heat_res, function(r) r$m1_rmse)
  heat_grid$m2_rmse <- sapply(heat_res, function(r) r$m2_rmse)

  # Identify optimal design points
  d_opt_row <- heat_grid[which.max(heat_grid$det_M), ]
  c_opt_row <- heat_grid[which.min(heat_grid$Mc),     ]
  r_opt_row <- heat_grid[which.min(heat_grid$c_rmse), ]

  # Summary table: criterion values at each key design point
  summary_tbl <- data.frame(
    Point     = c("D-optimal", "C-optimal", "RMSE-optimal"),
    Early_pct = round(c(d_opt_row$p1_pct,  c_opt_row$p1_pct,  r_opt_row$p1_pct),  1),
    Late_pct  = round(c(d_opt_row$p2_pct,  c_opt_row$p2_pct,  r_opt_row$p2_pct),  1),
    det_M     = round(c(d_opt_row$det_M,   c_opt_row$det_M,   r_opt_row$det_M),   2),
    Mc        = round(c(d_opt_row$Mc,      c_opt_row$Mc,      r_opt_row$Mc),      4),
    GLM_RMSE  = round(c(d_opt_row$c_rmse,  c_opt_row$c_rmse,  r_opt_row$c_rmse),  4),
    M1_RMSE   = round(c(d_opt_row$m1_rmse, c_opt_row$m1_rmse, r_opt_row$m1_rmse), 4),
    M2_RMSE   = round(c(d_opt_row$m2_rmse, c_opt_row$m2_rmse, r_opt_row$m2_rmse), 4)
  )

  cat("========================================\n")
  cat("SUMMARY TABLE | N_bunches =", N_bunches, "| n_sim =", n_sim, "\n")
  cat("========================================\n")
  print(summary_tbl, row.names = FALSE)
  cat("\n")

  sim1_heat_grid_results[[as.character(N_bunches)]] <- heat_grid
  sim1_heat_res_results[[as.character(N_bunches)]]  <- heat_res
}

cat("=== SIMULATION 1 COMPLETE ===\n")
cat("Results stored in sim1_heat_grid_results[[\"10\"]], [[\"20\"]], [[\"50\"]]\n\n")


# ============================================================
# SIMULATION 2 SETUP — Real-data design comparison
# Standard curve: logit(p) = scale_factor_t * t
# True Day50% at t = 0
# D-optimal design points: t = -1 (17.6%) and t = +1 (82.4%)
#
# Site parameters (alpha, beta) below are posterior means from
# fitting a Bayesian flowering model to real field data. The
# model-fitting step itself is not included here — only the
# fitted values are needed to reproduce this analysis. Any
# similar time-to-event dataset with a fitted intercept (alpha)
# and slope (beta) can be substituted in site_params below.
# ============================================================

# ------------------------------------------------------------
# FIXED PARAMETERS
# ------------------------------------------------------------
reference_date <- as.Date("2024-12-01")
scale_factor_t <- 1.5434

# D-optimal points on t_scaled
t_opt_early <- -1
t_opt_late  <-  1
p_opt_early <- plogis(t_opt_early * scale_factor_t)  # 17.6%
p_opt_late  <- plogis(t_opt_late  * scale_factor_t)  # 82.4%

# ------------------------------------------------------------
# SITE PARAMETERS
# alpha (beta0) and beta (beta_day): posterior means from the
# fitted flowering model (see paper Methods for model details)
# visits: actual field visit days (days from reference_date)
# ------------------------------------------------------------
sites <- c("WPN", "WSB", "LPN", "DJV")

site_params <- list(
  WPN = list(
    alpha  = 0.4173,
    beta   = 6.8788,
    visits = c(-18, -12, -9, -4, 2, 11, 16)
  ),
  WSB = list(
    alpha  = -4.0838,
    beta   = 6.3290,
    visits = c(-12, -4, 2, 5, 8, 11, 16)
  ),
  LPN = list(
    alpha  = -2.0372,
    beta   = 3.2705,
    visits = c(-18, -12, -5, 1, 8, 12, 16)
  ),
  DJV = list(
    alpha  = -3.3325,
    beta   = 5.5202,
    visits = c(-22, -5, 1, 5, 16)
  )
)

# day50 (days from reference_date) computed directly from alpha, beta
for (site in sites) {
  sp <- site_params[[site]]
  site_params[[site]]$day50 <- -10 * sp$alpha / sp$beta
}

# ------------------------------------------------------------
# CONVERSION FACTORS: t_scaled -> site-specific days
# days_per_tsc (dpsc) = 10 * scale_factor_t / beta, per site
# Multiplying an error/RMSE on the t_scaled axis by dpsc converts
# it to a real number of days for that site.
# ------------------------------------------------------------
dpsc_list <- setNames(
  lapply(sites, function(site)
    10 * scale_factor_t / site_params[[site]]$beta),
  sites
)

conv_df <- bind_rows(lapply(sites, function(site) {
  sp <- site_params[[site]]
  data.frame(
    Site         = site,
    alpha        = round(sp$alpha, 4),
    beta         = round(sp$beta,  4),
    day50_days   = round(sp$day50, 2),
    day50_date   = format(reference_date + round(sp$day50), "%d-%b-%Y"),
    days_per_tsc = round(dpsc_list[[site]], 3)
  )
}))

cat("========================================\n")
cat("SITE PARAMETERS AND days_per_tsc (dpsc) CONVERSION\n")
cat("dpsc = number of real days represented by one unit\n")
cat("of t_scaled, for each site (depends on beta)\n")
cat("========================================\n")
print(conv_df, row.names = FALSE)

# ------------------------------------------------------------
# FULL 2D TRIANGULAR GRID — on t_scaled
# Both visits range from t = -3 to t = +3
# filter(t2 > t1) gives the triangle
# ------------------------------------------------------------
t_grid <- seq(-3, 3, by = 0.1)
p_grid <- plogis(t_grid * scale_factor_t) * 100  # proportion %

heat_grid <- expand.grid(t1 = t_grid, t2 = t_grid) %>%
  filter(t2 > t1) %>%
  mutate(
    x1     = t1 * scale_factor_t,  # logit scale
    x2     = t2 * scale_factor_t,  # logit scale
    p1_pct = plogis(x1) * 100,
    p2_pct = plogis(x2) * 100
  )

cat("\nGrid points:", nrow(heat_grid), "\n")
cat("t_scaled range: ", min(t_grid), "to", max(t_grid), "\n")
cat("Proportion range:", round(min(p_grid), 1),
    "% to", round(max(p_grid), 1), "%\n")

# ------------------------------------------------------------
# ACTUAL VISITS IN t_SCALED — for reference/comparison
# ------------------------------------------------------------
site_visits_t <- lapply(sites, function(site) {
  sp      <- site_params[[site]]
  days    <- sp$visits
  p_visit <- plogis(sp$alpha + sp$beta * days / 10)
  t_visit <- logit(p_visit) / scale_factor_t
  data.frame(
    Site     = site,
    day      = days,
    p_pct    = round(p_visit * 100, 1),
    t_scaled = round(t_visit, 3)
  )
})
names(site_visits_t) <- sites

cat("\n========================================\n")
cat("ACTUAL VISITS IN t_SCALED\n")
cat("========================================\n")
for (site in sites) {
  cat("\n", site, "\n")
  print(site_visits_t[[site]], row.names = FALSE)
}

cat("\n=== SIMULATION 2 SETUP COMPLETE ===\n")
cat("Objects ready: sites, site_params, dpsc_list, heat_grid\n")
cat("  t_opt_early =", t_opt_early, "| t_opt_late =", t_opt_late, "\n")
cat("  p_opt_early =", round(p_opt_early * 100, 1),
    "% | p_opt_late =", round(p_opt_late * 100, 1), "%\n\n")


# ============================================================
# SIMULATION 2 — Real-data two-point visit design comparison
# Standard curve: logit(p) = scale_factor_t * t
# True Day50% at t = 0
# D-optimal at t = -1 (17.6%) and t = +1 (82.4%)
#
# For each grid point (t1, t2) on the full 2D triangular grid,
# computes:
#   1. log det(M) — D-optimality criterion (log scale for stability)
#   2. M.c        — C-optimality criterion (delta method for Day50%)
#   3. GLM RMSE   — Day50% estimation error
#   4. Method 1   — quick Day50% estimate on the proportion scale
#   5. Method 2   — quick Day50% estimate on the logit scale
# Repeated for N_bunches = 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100.
# Results are converted to site-specific days using dpsc_list.
#
# Uses: heat_grid, scale_factor_t, dpsc_list, sites, N
#       (all defined in the Simulation 2 setup above)
# ============================================================

set.seed(42)
n_sim2         <- 10000  # simulations per grid point
n_bunches_vec2 <- c(10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)

# ============================================================
# SIMULATION FUNCTION — standard curve
# logit(p) = scale_factor_t * t
# ============================================================
sim2_fn <- function(t1, t2, N_bunches, N, n_sim) {
  p1 <- plogis(t1 * scale_factor_t)
  p2 <- plogis(t2 * scale_factor_t)

  log_det_M_vec <- numeric(n_sim)
  Mc_vec        <- numeric(n_sim)
  c_rmse_vec    <- numeric(n_sim)
  m1_vec        <- numeric(n_sim)
  m2_vec        <- numeric(n_sim)
  p1_obs        <- numeric(n_sim)
  p2_obs        <- numeric(n_sim)

  x1 <- t1 * scale_factor_t
  x2 <- t2 * scale_factor_t

  for (i in seq_len(n_sim)) {
    y1 <- rbinom(N_bunches, N, p1)
    y2 <- rbinom(N_bunches, N, p2)

    p1_obs[i] <- mean(y1 / N)
    p2_obs[i] <- mean(y2 / N)

    df <- data.frame(
      y = c(y1, y2),
      n = N,
      x = rep(c(x1, x2), each = N_bunches)
    )

    g <- glm(cbind(y, n - y) ~ x, data = df, family = binomial())

    a <- g$coef[1]
    b <- g$coef[2]

    # D-optimality — log determinant for numerical stability
    M                <- solve(vcov(g))
    log_det_M_vec[i] <- determinant(M, logarithm = TRUE)$modulus

    # C-optimality: delta-method variance of Day50% = -a/b
    c.vec     <- c(-1/b, a/b^2)
    Mc_vec[i] <- as.numeric(t(c.vec) %*% solve(M) %*% c.vec)

    # GLM Day50% error — true Day50% at x = 0 (t = 0)
    c_rmse_vec[i] <- -a / b

    # Method 1: proportion scale
    b_m1      <- (p2_obs[i] - p1_obs[i]) / (x2 - x1)
    a_m1      <- p1_obs[i] - b_m1 * x1
    m1_vec[i] <- (0.5 - a_m1) / b_m1

    # Method 2: logit scale
    # Guard against p = 0 or p = 1, which makes logit blow up
    if (p1_obs[i] > 0 & p1_obs[i] < 1 &
        p2_obs[i] > 0 & p2_obs[i] < 1) {
      lp1       <- log(p1_obs[i] / (1 - p1_obs[i]))
      lp2       <- log(p2_obs[i] / (1 - p2_obs[i]))
      b_m2      <- (lp2 - lp1) / (x2 - x1)
      a_m2      <- lp1 - b_m2 * x1
      m2_vec[i] <- -a_m2 / b_m2
    } else {
      m2_vec[i] <- NA_real_
    }
  }

  list(
    log_det_M = mean(log_det_M_vec),
    Mc        = mean(Mc_vec),
    c_rmse    = sqrt(mean(c_rmse_vec^2)),
    m1_rmse   = sqrt(mean(m1_vec^2, na.rm = TRUE)),
    m2_rmse   = sqrt(mean(m2_vec^2, na.rm = TRUE)),
    p1_obs    = p1_obs,
    p2_obs    = p2_obs
  )
}

# ============================================================
# RUN SIMULATION 2 FOR EACH N_bunches VALUE
# Results kept in memory, keyed by N_bunches — nothing written to disk.
# ============================================================
sim2_heat_grid_results <- list()
sim2_heat_res_results  <- list()

for (N_bunches in n_bunches_vec2) {

  cat("Running Simulation 2 | N_bunches =", N_bunches,
      "| n_sim =", n_sim2, "| grid points =", nrow(heat_grid), "\n")

  hg  <- heat_grid
  res <- vector("list", nrow(hg))

  for (i in seq_len(nrow(hg))) {
    if (i %% 50 == 0) cat("  point", i, "of", nrow(hg), "\n")
    res[[i]] <- sim2_fn(hg$t1[i], hg$t2[i], N_bunches, N, n_sim2)
  }

  hg$log_det_M <- sapply(res, function(r) r$log_det_M)
  hg$Mc        <- sapply(res, function(r) r$Mc)
  hg$c_rmse    <- sapply(res, function(r) r$c_rmse)
  hg$m1_rmse   <- sapply(res, function(r) r$m1_rmse)
  hg$m2_rmse   <- sapply(res, function(r) r$m2_rmse)

  # Site-specific day conversion
  for (site in sites) {
    dpsc <- dpsc_list[[site]]
    hg[[paste0("c_rmse_days_",  site)]] <- hg$c_rmse  * dpsc
    hg[[paste0("m1_rmse_days_", site)]] <- hg$m1_rmse * dpsc
    hg[[paste0("m2_rmse_days_", site)]] <- hg$m2_rmse * dpsc
  }

  # Identify optimal design points
  d_opt_row <- hg[which.max(hg$log_det_M), ]
  c_opt_row <- hg[which.min(hg$Mc),        ]
  r_opt_row <- hg[which.min(hg$c_rmse),    ]

  summary_tbl <- data.frame(
    Point     = c("D-optimal", "C-optimal", "Empirical Day50%"),
    t1        = round(c(d_opt_row$t1,        c_opt_row$t1,        r_opt_row$t1),        2),
    t2        = round(c(d_opt_row$t2,        c_opt_row$t2,        r_opt_row$t2),        2),
    p1_pct    = round(c(d_opt_row$p1_pct,    c_opt_row$p1_pct,    r_opt_row$p1_pct),    1),
    p2_pct    = round(c(d_opt_row$p2_pct,    c_opt_row$p2_pct,    r_opt_row$p2_pct),    1),
    log_det_M = round(c(d_opt_row$log_det_M, c_opt_row$log_det_M, r_opt_row$log_det_M), 4),
    Mc        = round(c(d_opt_row$Mc,        c_opt_row$Mc,        r_opt_row$Mc),        4),
    GLM_RMSE  = round(c(d_opt_row$c_rmse,    c_opt_row$c_rmse,    r_opt_row$c_rmse),    4),
    M1_RMSE   = round(c(d_opt_row$m1_rmse,   c_opt_row$m1_rmse,   r_opt_row$m1_rmse),   4),
    M2_RMSE   = round(c(d_opt_row$m2_rmse,   c_opt_row$m2_rmse,   r_opt_row$m2_rmse),   4)
  )

  cat("========================================\n")
  cat("SUMMARY TABLE | N_bunches =", N_bunches, "| n_sim =", n_sim2, "\n")
  cat("========================================\n")
  print(summary_tbl, row.names = FALSE)

  cat("\nSite-specific days at optimal points:\n")
  for (site in sites) {
    dpsc <- dpsc_list[[site]]
    cat("\n", site, "| days_per_tsc =", round(dpsc, 3), "\n")
    site_tbl <- data.frame(
      Point        = c("D-optimal", "C-optimal", "Empirical Day50%"),
      c_rmse_days  = round(c(d_opt_row$c_rmse,  c_opt_row$c_rmse,  r_opt_row$c_rmse)  * dpsc, 3),
      m1_rmse_days = round(c(d_opt_row$m1_rmse, c_opt_row$m1_rmse, r_opt_row$m1_rmse) * dpsc, 3),
      m2_rmse_days = round(c(d_opt_row$m2_rmse, c_opt_row$m2_rmse, r_opt_row$m2_rmse) * dpsc, 3)
    )
    print(site_tbl, row.names = FALSE)
  }
  cat("\n")

  sim2_heat_grid_results[[as.character(N_bunches)]] <- hg
  sim2_heat_res_results[[as.character(N_bunches)]]  <- res
}

cat("=== SIMULATION 2 COMPLETE ===\n")
cat("Results stored in sim2_heat_grid_results[[\"10\"]], [[\"15\"]], ..., [[\"100\"]]\n")





# ============================================================
# ACTUAL-VISIT RMSE — Real data
# Computes RMSE using the real, actual visit days per site
# (as opposed to the optimal grid search in Simulation 2).
# Serves as a reference point: how well do the real visit
# schedules perform compared to the D-/C-optimal designs?
#
# N_bunches fixed at 20 (matches the actual field dataset)
# Evaluation grid matches the Simulation 2 grid range (t = -3 to 3)
#
# The GLM here uses x_vis (actual visit days, logit scale) as the
# predictor. On that scale, true Day50% = 0, and the true position
# of each evaluation proportion is simply x_eval. dpsc converts the
# RMSE from t_scaled units to site-specific days at the end.
#
# Uses: site_params, sites, dpsc_list, scale_factor_t, N
#       (all defined in the Simulation 2 setup)
# ============================================================

set.seed(42)
n_sim_actual     <- 10000  # simulations — matches Simulation 2
N_bunches_actual <- 20     # fixed — matches the actual field dataset

# ============================================================
# SIMULATION FUNCTION — actual visits per site
# ============================================================
actual_rmse_fn <- function(site, N_bunches, N, n_sim) {
  sp    <- site_params[[site]]
  days  <- sp$visits
  n_vis <- length(days)

  # True proportions at actual visit days
  p_vis <- plogis(sp$alpha + sp$beta * days / 10)
  # Convert actual visits to logit scale
  x_vis <- log(p_vis / (1 - p_vis))

  # Evaluation proportions — matches Simulation 2 grid range
  t_eval <- seq(-3, 3, by = 0.1)
  p_eval <- plogis(t_eval * scale_factor_t)
  x_eval <- log(p_eval / (1 - p_eval))

  # True position of each evaluation proportion on the fitted scale
  x_eval_true <- x_eval

  c_rmse_vec <- numeric(n_sim)
  d_rmse_vec <- numeric(n_sim)

  for (i in seq_len(n_sim)) {
    y_obs <- rbinom(n_vis * N_bunches, N,
                     rep(p_vis, each = N_bunches))
    df <- data.frame(
      y = y_obs,
      n = N,
      x = rep(x_vis, each = N_bunches)
    )
    g <- glm(cbind(y, n - y) ~ x, data = df, family = binomial())
    a <- g$coef[1]
    b <- g$coef[2]

    # C-RMSE: error in estimating the 50% point
    # Estimated = -a/b, true = 0 on the x_vis scale
    c_rmse_vec[i] <- -a / b

    # D-RMSE: timing error across evaluation proportions
    t_hat         <- (x_eval - a) / b
    d_rmse_vec[i] <- mean((t_hat - x_eval_true)^2)
  }

  dpsc <- dpsc_list[[site]]
  data.frame(
    Site        = site,
    n_visits    = n_vis,
    N_bunches   = N_bunches,
    N           = N,
    c_rmse      = round(sqrt(mean(c_rmse_vec^2)), 4),
    d_rmse      = round(sqrt(mean(d_rmse_vec)),   4),
    c_rmse_days = round(sqrt(mean(c_rmse_vec^2)) * dpsc, 3),
    d_rmse_days = round(sqrt(mean(d_rmse_vec))   * dpsc, 3)
  )
}

# ============================================================
# RUN FOR ALL SITES
# ============================================================
actual_rmse_real <- bind_rows(lapply(sites, function(site) {
  cat("Running actual-visit RMSE for", site, "...\n")
  actual_rmse_fn(site, N_bunches_actual, N, n_sim_actual)
}))

# ============================================================
# PRINT TABLE
# ============================================================
cat("\n========================================\n")
cat("ACTUAL-VISIT RMSE — REAL DATA\n")
cat("N =", N, "| N_bunches =", N_bunches_actual,
    "| n_sim =", n_sim_actual, "\n")
cat("========================================\n")
print(actual_rmse_real, row.names = FALSE)
