#' Maximum Score Estimation and 95% Confidence Interval Calculation
#' @param y Numeric vector of binary dependent variable (0 or 1).
#' @param X Numeric matrix of independent variables.
#' @param B Number of bootstrap replications (default 50).
#' @param plugin_bw Bandwidth for plugin estimator (default 2.0).
#' @param num_eps Epsilon for numerical derivative (default 0.5).
#' @param n_starts Number of DEoptim random restarts for global optimization (default 3).
#' @param x1_coef Coefficient for the first column of X (must be either 1 or -1; default is 1).
#' @return A list containing the point estimates and 95% confidence intervals.
#' @export
maxscore_ci <- function(y, X, B = 50, plugin_bw = 2.0, num_eps = 0.5, n_starts = 3, x1_coef = 1) {
  if (!(x1_coef %in% c(1, -1))) {
    stop("Parameter 'x1_coef' must be either 1 or -1.")
  }
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  k_dim <- ncol(X) - 1
  
  # 1. Point estimation using DEoptim
  bhat <- fit_maxscore_de(y, X, k = k_dim, n_starts = n_starts, x1_coef = x1_coef)
  
  # 2. Compute H matrices using C++ routines
  H.ms <- hat_H_ms(y, X, plugin_bw, bhat, x1_coef = x1_coef)
  H.nd <- hat_H_nd(y, X, bhat, num_eps, x1_coef = x1_coef)
  
  plugin.boot <- matrix(NA, nrow = B, ncol = k_dim)
  std.boot    <- matrix(NA, nrow = B, ncol = k_dim)
  
  for (b in 1:B) {
    cat("bootstrap iteration...:",b,"\n")
    boot.idx <- sample(n, replace = TRUE)
    boot.y <- y[boot.idx]
    boot.X <- X[boot.idx, , drop = FALSE]
    
    # Standard bootstrap estimation
    std.res <- fit_maxscore_de(boot.y, boot.X, k = k_dim, n_starts = n_starts, x1_coef = x1_coef)
    std.boot[b, ] <- std.res
    
    # Reshaped Bootstrap: Plug-in
    p.res <- fit_boot_de(y, X, boot.y, boot.X, H.ms, bhat, k = k_dim, n_starts = n_starts, x1_coef = x1_coef)
    plugin.boot[b, ] <- p.res
  }
  
  # Compute 95% Confidence Intervals using Plugin Reshaped Bootstrap
  ci_lower <- apply(plugin.boot - matrix(bhat, nrow = B, ncol = k_dim, byrow = TRUE), 2, function(x) quantile(x, 0.025, type = 1))
  ci_upper <- apply(plugin.boot - matrix(bhat, nrow = B, ncol = k_dim, byrow = TRUE), 2, function(x) quantile(x, 0.975, type = 1))
  
  ci_matrix <- cbind(
    Estimate = bhat,
    CI_Lower = bhat + ci_lower,
    CI_Upper = bhat + ci_upper
  )
  
  list(
    point_estimate = bhat,
    confidence_interval_95 = ci_matrix,
    plugin_boot_samples = plugin.boot,
    std_boot_samples = std.boot
  )
}

#' Maximum Score Point Estimator using DEoptim
#' @export
fit_maxscore_de <- function(y, X, k = ncol(X) - 1, lower = rep(-50.0, k), upper = rep(50.0, k), n_starts = 3, x1_coef = 1) {
  best_val <- -Inf
  best_par <- NULL
  
  obj_fnc <- function(beta) {
    -score_fnc_cpp(y, X, beta, x1_coef = x1_coef)
  }
  
  for (run in 1:n_starts) {
    ctrl <- DEoptim.control(
      NP = 1500, 
      itermax = 5000, 
      strategy = 1,
      trace = FALSE
    )
    res <- DEoptim(fn = obj_fnc, lower = lower, upper = upper, control = ctrl)
    
    val <- -res$optim$bestval
    if (val > best_val) {
      best_val <- val
      best_par <- res$optim$bestmem
    }
  }
  return(as.numeric(best_par))
}

#' Reshaped Bootstrap Estimator using DEoptim
#' @export
fit_boot_de <- function(y, X, bootY, bootX, H, beta_hat, k = length(beta_hat), lower = rep(-50.0, k), upper = rep(50.0, k), n_starts = 3, x1_coef = 1) {
  best_val <- -Inf
  best_par <- NULL
  
  obj_fnc <- function(beta) {
    -boot_score_cpp(y, X, bootY, bootX, beta, H, beta_hat, x1_coef = x1_coef)
  }
  
  for (run in 1:n_starts) {
    ctrl <- DEoptim.control(
      NP = 1500, 
      itermax = 5000,
      strategy = 1,
      trace = FALSE
    )
    res <- DEoptim(fn = obj_fnc, lower = lower, upper = upper, control = ctrl)
    
    val <- -res$optim$bestval
    if (val > best_val) {
      best_val <- val
      best_par <- res$optim$bestmem
    }
  }
  return(as.numeric(best_par))
}
