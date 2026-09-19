# this script implements regsensitivy in R
# based on Diegert, Masten, Poirier (2026)
# https://arxiv.org/abs/2206.02303

get_dgp <- function(y, x, w, weight) {
  vw <- if (is.null(weight)) cov(w) else cov.wt(w, wt = weight)$cov
  keep <- diag(vw) != 0
  w <- w[, keep, drop = FALSE]

  data <- cbind(y, x, w)
  
  v <- if (is.null(weight)) cov(data) else cov.wt(data, wt = weight)$cov
  
  s.var_y <- v[1, 1]
  s.var_x <- v[2, 2]
  s.var_w <- v[3:nrow(v), 3:ncol(v)]
  
  s.wt <- chol2inv(chol(s.var_w))
  
  covwx <- as.matrix(v[3:nrow(v), 2])
  covwy <- as.matrix(v[3:nrow(v), 1])
  covxy <- v[1, 2]
  
  s.k0 <- s.var_x - t(covwx) %*% s.wt %*% covwx
  s.k1 <- covxy - t(covwx) %*% s.wt %*% covwy
  s.k2 <- s.var_y - t(covwy) %*% s.wt %*% covwy
  
  s.covwx_norm_sq <- s.var_x - s.k0
  
  s.var_x_resid <- s.k0
  
  s.beta_short <- covxy / s.var_x
  s.beta_med <- s.k1 / s.k0
  
  s.gamma_med <- (covwy - c(s.beta_med) * covwx)
  s.pi_med <- covwx
  
  s.r_short <- s.beta_short ^ 2 * s.var_x / s.var_y
  s.r_med <- (
    s.beta_med ^ 2 * s.var_x
    + t(s.gamma_med) %*% s.wt %*% s.gamma_med
    + 2 * s.beta_med %*% t(s.gamma_med) %*% s.wt %*% covwx
  ) / s.var_y
  
  s <- list(
    var_y = s.var_y,
    var_x = s.var_x,
    var_w = s.var_w,
    wt = s.wt,
    k0 = s.k0,
    k1 = s.k1,
    k2 = s.k2,
    covwx_norm_sq = s.covwx_norm_sq,
    var_x_resid = s.var_x_resid,
    beta_short = s.beta_short,
    beta_med = s.beta_med,
    gamma_med = s.gamma_med,
    pi_med = s.pi_med,
    r_short = s.r_short,
    r_med = s.r_med
  )
  class(s) <- "dmp objects"
  
  return(s)
}

zmax <- function(c, rx, s) {
  cmax <- min(c(c, rx))
  
  z <- sqrt(s$covwx_norm_sq) * rx * sqrt(1 - cmax ^ 2)
  z <- z / (1 - rx * cmax)
  
  return(z)
}

max_beta_bound <- function(c, s) {
  if (c == 1) {
    return(sqrt(s$k0 / s$var_x))
  }
  
  A = c ^ 2 * (s$k0 + s$covwx_norm_sq) - s$covwx_norm_sq
  B = s$k0 * c
  C = s$k0
  
  root1 <- (B + sqrt(B ^ 2 - A * C)) / A
  root2 <- (B - sqrt(B ^ 2 - A * C)) / A
  
  if (0 <= root1 & root1 <= root2) {
    return(root1)
  } else {
    return(root2)
  }
}

beta_deviation <- function(z, s) {
  z_sq = z ^ 2
  z_sq = min(c(z_sq, s$k0 - .000001))
  deviation_sq = (z_sq * (s$k2 / s$k0 - (s$k1 / s$k0) ^ 2)) / (s$k0 - z_sq)
  deviation = sqrt(deviation_sq)
  return(deviation)
}

beta_bounds <- function(c, rx, s) {
  finite_threshold <- max_beta_bound(c, s)
  
  bounds <- matrix(, nrow = length(rx), ncol = 2)
  for (i in 1:length(rx)) {
    finite <- rx[i] < finite_threshold
    if (finite == TRUE) {
      z <- zmax(c, rx[i], s)
      dev <- beta_deviation(z, s)
      bounds[i, 1:2] <- c(s$beta_med - dev, s$beta_med + dev)
    } else {
      bounds[i, 1:2] <- c(-Inf, Inf)
    }
  }
  return(bounds)
}

regsensitivity <- function(formula, data, weights, subset, na.action) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mm <- eval(mf, parent.frame())

  rhs_vars <- all.vars(formula[[3]])
  x_var <- rhs_vars[1]
  w_vars <- rhs_vars[-1]

  y  <- model.response(mm)
  x  <- mm[[x_var]]
  w  <- mm[w_vars]
  wt <- model.weights(mm)

  s <- get_dgp(y, x, w, wt)
  rx <- seq(0, 1, by = 0.1)
  bounds <- beta_bounds(1, rx, s)
  cbind(rx, bounds)
}
