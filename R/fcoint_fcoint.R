#' Fourier Cointegration Tests for Time Series
#'
#' Residual-based and ADL-based cointegration tests in which smooth
#' structural change in the deterministic component is approximated by a
#' single-frequency Fourier function: the Fourier ADL test (FADL), the
#' Fourier Engle-Granger test (FEG), a covariate-augmented variant (FEG2),
#' and the Fourier test for the null of cointegration of Tsong et al.
#' (2016).
#'
#' @param y Numeric vector. Dependent variable.
#' @param x Numeric vector or matrix. Regressors (one column per variable).
#' @param test Character. \code{"fadl"}, \code{"feg"}, \code{"feg2"},
#'   \code{"tsong"} or \code{"all"}. Default \code{"fadl"}.
#' @param model Character. \code{"constant"} (default) or \code{"trend"}.
#' @param max_freq Integer. Largest Fourier frequency searched (1 to 5; at
#'   most 3 is used for the Tsong test). Default \code{5}.
#' @param max_lag Integer. Largest lag order. \code{0} (default) selects a
#'   default: 6 for FADL and \code{floor(12 * (T/100)^0.25)} for FEG and
#'   FEG2.
#' @param criterion Character. Information criterion used by FADL to choose
#'   the lag orders: \code{"aic"} (default) or \code{"bic"}.
#' @param adf_lags Character. Lag selection in the ADF regressions of FEG
#'   and FEG2: \code{"tsig"} (default; general-to-specific, dropping the
#'   last lag while its absolute t-ratio is below 1.645), \code{"aic"} or
#'   \code{"bic"}.
#' @param dols_lags Integer. Leads and lags of the differenced regressors in
#'   the DOLS regression of the Tsong test. \code{0} (default) uses
#'   \code{floor(4 * (T/100)^(2/9))}.
#' @param lrv Character. Long-run variance estimator in the Tsong test:
#'   \code{"iid"} (default) or \code{"bartlett"}.
#' @param bandwidth Integer. Bartlett bandwidth for \code{lrv = "bartlett"}.
#'   \code{NULL} (default) uses \code{round(4 * (T/100)^(2/9))}.
#'
#' @return An object of class \code{"fcoint"}: a list with \code{test},
#'   \code{results} (one list per test), \code{model}, \code{criterion} and
#'   \code{nobs}.
#'
#' @details
#' The Fourier terms are \eqn{\sin(2\pi k t/T)} and \eqn{\cos(2\pi k t/T)},
#' \eqn{t = 1, \dots, T}. In every test the frequency is the one that
#' minimises the sum of squared residuals.
#'
#' \strong{FADL} (Banerjee, Arcabic and Lee, 2017). The ADL regression
#' \deqn{\Delta y_t = d(t) + \delta y_{t-1} + \gamma' x_{t-1} +
#' \sum_{j=0}^{q} \phi_j' \Delta x_{t-j} + \sum_{i=1}^{p} \alpha_i
#' \Delta y_{t-i} + e_t}
#' is estimated on a common sample. For each frequency the lag orders
#' \eqn{p, q \in \{1, \dots, }\code{max_lag}\eqn{\}} minimise the chosen
#' information criterion (one \eqn{q} per regressor when there are at most
#' two regressors, a common \eqn{q} otherwise); the statistic is the
#' t-ratio of \eqn{\delta}. Critical values are those of the Stata module
#' \code{fcoint}, which attributes them to Tables 1a and 1b of the paper
#' (\eqn{T} = 100, 500, 2500; interpolated linearly in \eqn{T} between 100
#' and 500).
#'
#' \strong{FEG} (Yilanci, 2019). The regression of \eqn{y_t} on a constant
#' (and a trend), the Fourier terms and \eqn{x_t} gives residuals
#' \eqn{\hat u_t}; the statistic is the t-ratio of \eqn{\rho} in
#' \eqn{\Delta\hat u_t = \rho\hat u_{t-1} + \sum_{i=1}^{p}\gamma_i
#' \Delta\hat u_{t-i} + \varepsilon_t} (no deterministic terms). Critical
#' values are Table 1 of Yilanci (2019) for one to three regressors and
#' frequencies 1 to 5 (\eqn{T} = 100, 500, 1000; linear interpolation in
#' \eqn{T}). Two entries of that table carry evident typesetting errors and
#' are corrected here: the 5\% value for two regressors, \eqn{k = 3}, trend
#' model, \eqn{T = 500} is printed as 4.620 (sign lost), and the 5\% and
#' 10\% values for three regressors, \eqn{k = 4}, constant model,
#' \eqn{T = 1000} are printed in reverse order.
#'
#' \strong{FEG2}. As FEG, but the ADF regression includes a constant and
#' \eqn{\Delta x_t}, and the long-run squared correlation \eqn{\rho^2}
#' between the FEG and FEG2 errors is reported. The Stata module attributes
#' this test and its critical values to an unpublished working paper by
#' Banerjee and Lee, whose tables could not be verified (the Stata table
#' also fills several cells by copying other cells). No critical values are
#' therefore reported for FEG2.
#'
#' \strong{Tsong et al. (2016)}. Following the TSPDLIB implementation of
#' Nazlioglu, the regression of \eqn{y_t} on a constant (and a trend),
#' \eqn{x_t} and the Fourier terms is estimated by OLS and by DOLS, and
#' \eqn{CI_f = T^{-2}\sum_t S_t^2/\hat\omega^2} with \eqn{S_t} the partial
#' sums of the residuals. Large values reject the null of cointegration.
#' Critical values (DOLS statistic) are the TSPDLIB table for up to four
#' regressors and \eqn{k \le 3}; the entry for four regressors,
#' \eqn{k = 3}, constant model repeats the \eqn{k = 2} row in TSPDLIB and
#' differs from the Stata module, so it is reported as \code{NA}. The
#' F statistics for the Fourier terms are reported without critical values.
#'
#' @references
#' Banerjee, P., Arcabic, V. and Lee, H. (2017). Fourier ADL cointegration
#' test to approximate smooth breaks with new evidence from crude oil
#' market. \emph{Economic Modelling}, 67, 114-124.
#' \doi{10.1016/j.econmod.2016.11.004}
#'
#' Tsong, C.-C., Lee, C.-F., Tsai, L.-J. and Hu, T.-C. (2016). The Fourier
#' approximation and testing for the null of cointegration.
#' \emph{Empirical Economics}, 51(3), 1085-1113.
#' \doi{10.1007/s00181-015-1028-6}
#'
#' Yilanci, V. (2019). A residual-based cointegration test with a Fourier
#' approximation. MPRA Paper No. 95395, University Library of Munich.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' x <- cumsum(rnorm(n))
#' y <- 1 + 0.5 * x + sin(2 * pi * (1:n) / n) + rnorm(n, sd = 0.3)
#' fcoint(y, x, test = "feg", max_freq = 3)
#'
#' @export
fcoint <- function(y, x,
                   test      = c("fadl", "feg", "feg2", "tsong", "all"),
                   model     = c("constant", "trend"),
                   max_freq  = 5L,
                   max_lag   = 0L,
                   criterion = c("aic", "bic"),
                   adf_lags  = c("tsig", "aic", "bic"),
                   dols_lags = 0L,
                   lrv       = c("iid", "bartlett"),
                   bandwidth = NULL) {

  test      <- match.arg(test)
  model     <- match.arg(model)
  criterion <- match.arg(criterion)
  adf_lags  <- match.arg(adf_lags)
  lrv       <- match.arg(lrv)

  y <- as.numeric(y)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(y) != nrow(x))
    stop("'y' and 'x' must have the same number of observations.")
  if (anyNA(y) || anyNA(x))
    stop("'y' and 'x' must not contain missing values.")

  max_freq  <- as.integer(max_freq)
  max_lag   <- as.integer(max_lag)
  dols_lags <- as.integer(dols_lags)
  if (max_freq < 1L || max_freq > 5L)
    stop("'max_freq' must be between 1 and 5.")
  if (max_lag < 0L) stop("'max_lag' must be non-negative.")

  T_obs <- length(y)
  if (T_obs < 30L)
    stop("Insufficient observations (T = ", T_obs, "). Need at least 30.")

  results <- list()
  if (test %in% c("fadl", "all"))
    results[["fadl"]] <- .fcoint_fadl(y, x, model, max_freq, max_lag,
                                      criterion)
  if (test %in% c("feg", "all"))
    results[["feg"]] <- .fcoint_feg(y, x, model, max_freq, max_lag,
                                    adf_lags, augment = FALSE)
  if (test %in% c("feg2", "all"))
    results[["feg2"]] <- .fcoint_feg(y, x, model, max_freq, max_lag,
                                     adf_lags, augment = TRUE)
  if (test %in% c("tsong", "all"))
    results[["tsong"]] <- .fcoint_tsong(y, x, model, min(max_freq, 3L),
                                        dols_lags, lrv, bandwidth)

  structure(list(test = test, results = results, model = model,
                 criterion = criterion, nobs = T_obs),
            class = "fcoint")
}


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

# internal
.lagv <- function(v, j) {
  n <- length(v)
  if (j == 0L) return(v)
  if (j > 0L) c(rep(NA_real_, j), v[seq_len(n - j)])
  else c(v[(1 - j):n], rep(NA_real_, -j))
}

# internal
.fourier_terms <- function(T_obs, k) {
  tt <- seq_len(T_obs)
  cbind(sin(2 * pi * k * tt / T_obs), cos(2 * pi * k * tt / T_obs))
}

# internal
.det_terms <- function(T_obs, model) {
  if (model == "trend") cbind(1, seq_len(T_obs)) else matrix(1, T_obs, 1L)
}

# internal
.ols <- function(X, y) {
  qx <- qr(X)
  if (qx$rank < ncol(X)) return(NULL)
  b <- qr.coef(qx, y)
  e <- as.numeric(y - X %*% b)
  n <- length(y)
  s2 <- sum(e^2) / (n - ncol(X))
  XtXi <- chol2inv(qr.R(qx))
  XtXi[qx$pivot, qx$pivot] <- XtXi
  list(b = as.numeric(b), e = e, ssr = sum(e^2), n = n, k = ncol(X),
       se = sqrt(diag(XtXi) * s2))
}

# internal
.ic_value <- function(ssr, n, k, criterion) {
  if (criterion == "bic") log(ssr / n) + log(n) * k / n
  else log(ssr / n) + 2 * k / n
}


# ------------------------------------------------------------------
# FADL
# ------------------------------------------------------------------

# internal
.fcoint_fadl <- function(y, x, model, max_freq, max_lag, criterion) {
  T_obs <- length(y)
  n_x   <- ncol(x)
  nl    <- if (max_lag == 0L) 6L else min(max_lag, 6L)

  dy  <- c(NA_real_, diff(y))
  dx  <- rbind(NA_real_, apply(x, 2, diff))
  if (!is.matrix(dx)) dx <- matrix(dx, ncol = n_x)
  ylag <- .lagv(y, 1L)
  xlag <- apply(x, 2, .lagv, j = 1L)
  if (!is.matrix(xlag)) xlag <- matrix(xlag, ncol = n_x)
  dylags <- sapply(seq_len(nl), function(j) .lagv(dy, j))
  dxlags <- lapply(seq_len(n_x), function(i)
    sapply(0:nl, function(j) .lagv(dx[, i], j)))

  ## common estimation sample: all lags up to nl available
  ok <- seq(nl + 2L, T_obs)
  det <- .det_terms(T_obs, model)

  grid <- if (n_x <= 2L) {
    as.matrix(expand.grid(c(list(p = seq_len(nl)),
                            stats::setNames(rep(list(seq_len(nl)), n_x),
                                            paste0("q", seq_len(n_x))))))
  } else {
    g <- expand.grid(p = seq_len(nl), q = seq_len(nl))
    as.matrix(cbind(g$p, matrix(g$q, nrow(g), n_x)))
  }

  best <- list(ssr = Inf)
  for (k in seq_len(max_freq)) {
    fk <- .fourier_terms(T_obs, k)
    kbest <- list(ic = Inf)
    for (g in seq_len(nrow(grid))) {
      p <- grid[g, 1L]
      qs <- grid[g, -1L]
      X <- cbind(det, fk, ylag, xlag,
                 do.call(cbind, lapply(seq_len(n_x), function(i)
                   dxlags[[i]][, seq_len(qs[i] + 1L), drop = FALSE])),
                 dylags[, seq_len(p), drop = FALSE])
      fit <- .ols(X[ok, , drop = FALSE], dy[ok])
      if (is.null(fit)) next
      ic <- .ic_value(fit$ssr, fit$n, fit$k, criterion)
      if (ic < kbest$ic) {
        pos <- ncol(det) + 3L
        kbest <- list(ic = ic, ssr = fit$ssr, p = p, q = qs,
                      delta = fit$b[pos], se = fit$se[pos], n = fit$n)
      }
    }
    if (is.finite(kbest$ic) && kbest$ssr < best$ssr)
      best <- c(kbest, k = k)
  }
  if (!is.finite(best$ssr)) return(list(error = "FADL regression could not be estimated"))

  cv <- .fcoint_cv_fadl(n_x, best$k, T_obs, model)
  list(test = "fadl", tstat = best$delta / best$se, delta = best$delta,
       se_delta = best$se, frequency = best$k, lag = best$p,
       lag_dx = best$q, ssr = best$ssr, nobs = best$n,
       cv1 = cv[1L], cv5 = cv[2L], cv10 = cv[3L],
       model = model, criterion = criterion)
}


# ------------------------------------------------------------------
# FEG and FEG2
# ------------------------------------------------------------------

# internal
.fcoint_feg <- function(y, x, model, max_freq, max_lag, adf_lags, augment) {
  T_obs <- length(y)
  n_x   <- ncol(x)
  det   <- .det_terms(T_obs, model)
  pmax  <- if (max_lag == 0L) max(1L, floor(12 * (T_obs / 100)^0.25)) else max_lag

  ssr_k <- vapply(seq_len(max_freq), function(k) {
    f <- .ols(cbind(det, .fourier_terms(T_obs, k), x), y)
    if (is.null(f)) Inf else f$ssr
  }, numeric(1))
  k_star <- which.min(ssr_k)
  u <- .ols(cbind(det, .fourier_terms(T_obs, k_star), x), y)$e

  du   <- c(NA_real_, diff(u))
  ulag <- .lagv(u, 1L)
  dulags <- sapply(seq_len(pmax), function(j) .lagv(du, j))
  dulags <- matrix(dulags, nrow = T_obs)
  dx <- rbind(NA_real_, apply(x, 2, diff))
  if (!is.matrix(dx)) dx <- matrix(dx, ncol = n_x)

  build <- function(p) {
    X <- cbind(ulag)
    if (augment) X <- cbind(1, X, dx)
    if (p > 0L) X <- cbind(X, dulags[, seq_len(p), drop = FALSE])
    X
  }
  pos <- if (augment) 2L else 1L

  if (adf_lags == "tsig") {
    p_opt <- 0L
    for (p in pmax:1L) {
      rows <- seq(p + 2L, T_obs)
      f <- .ols(build(p)[rows, , drop = FALSE], du[rows])
      if (is.null(f)) next
      if (abs(f$b[ncol(build(p))] / f$se[ncol(build(p))]) >= 1.645) {
        p_opt <- p
        break
      }
    }
  } else {
    rows <- seq(pmax + 2L, T_obs)
    ics <- vapply(0:pmax, function(p) {
      f <- .ols(build(p)[rows, , drop = FALSE], du[rows])
      if (is.null(f)) Inf else .ic_value(f$ssr, f$n, f$k, adf_lags)
    }, numeric(1))
    p_opt <- which.min(ics) - 1L
  }

  rows <- seq(p_opt + 2L, T_obs)
  fit  <- .ols(build(p_opt)[rows, , drop = FALSE], du[rows])
  tstat <- fit$b[pos] / fit$se[pos]

  out <- list(test = if (augment) "feg2" else "feg", tstat = tstat,
              delta = fit$b[pos], se_delta = fit$se[pos],
              frequency = k_star, lag = p_opt, nobs = fit$n,
              model = model)

  if (!augment) {
    cv <- .fcoint_cv_feg(n_x, k_star, T_obs, model)
  } else {
    ## long-run squared correlation between FEG and FEG2 errors
    Xv <- cbind(ulag, if (p_opt > 0L) dulags[, seq_len(p_opt), drop = FALSE])
    v  <- .ols(Xv[rows, , drop = FALSE], du[rows])$e
    e  <- fit$e
    v  <- v - mean(v); e <- e - mean(e)
    n  <- length(e)
    bw <- floor(T_obs^(1/3))
    s_v <- sum(v^2) / n; s_e <- sum(e^2) / n; s_ve <- sum(v * e) / n
    for (j in seq_len(bw)) {
      w <- 1 - j / (bw + 1)
      s_v  <- s_v  + 2 * w * sum(v[(j + 1):n] * v[1:(n - j)]) / n
      s_e  <- s_e  + 2 * w * sum(e[(j + 1):n] * e[1:(n - j)]) / n
      s_ve <- s_ve + w * (sum(v[(j + 1):n] * e[1:(n - j)]) +
                          sum(e[(j + 1):n] * v[1:(n - j)])) / n
    }
    out$rho2 <- if (s_v > 0 && s_e > 0) min(1, s_ve^2 / (s_v * s_e)) else NA_real_
    cv <- rep(NA_real_, 3L)
  }
  out$cv1 <- cv[1L]; out$cv5 <- cv[2L]; out$cv10 <- cv[3L]
  out
}


# ------------------------------------------------------------------
# Tsong et al. (2016)
# ------------------------------------------------------------------

# internal
.lrv_est <- function(e, type, bw) {
  n <- length(e)
  s <- sum(e^2) / n
  if (type == "bartlett" && bw > 0) {
    for (j in seq_len(min(bw, n - 1L)))
      s <- s + 2 * (1 - j / (bw + 1)) * sum(e[1:(n - j)] * e[(1 + j):n]) / n
  }
  s
}

# internal
.fcoint_tsong <- function(y, x, model, kmax, dols_lags, lrv, bandwidth) {
  T_obs <- length(y)
  n_x   <- ncol(x)
  q  <- if (dols_lags == 0L) floor(4 * (T_obs / 100)^(2/9)) else dols_lags
  bw <- if (is.null(bandwidth)) round(4 * (T_obs / 100)^(2/9)) else as.integer(bandwidth)
  det <- .det_terms(T_obs, model)

  dx <- rbind(0, apply(x, 2, diff))
  if (!is.matrix(dx)) dx <- matrix(dx, ncol = n_x)
  leads <- do.call(cbind, lapply(seq_len(n_x), function(i)
    sapply(seq_len(q), function(j) .lagv(dx[, i], -j))))
  lags  <- do.call(cbind, lapply(seq_len(n_x), function(i)
    sapply(seq_len(q), function(j) .lagv(dx[, i], j))))
  keep <- seq(q + 2L, T_obs - q)

  res <- lapply(seq_len(kmax), function(k) {
    fk <- .fourier_terms(T_obs, k)
    z  <- cbind(det, x, fk)
    zr <- cbind(det, x)
    f1 <- .ols(z, y); r1 <- .ols(zr, y)
    zd  <- cbind(z, leads, dx, lags)[keep, , drop = FALSE]
    zrd <- cbind(zr, leads, dx, lags)[keep, , drop = FALSE]
    f2 <- .ols(zd, y[keep]); r2 <- .ols(zrd, y[keep])
    if (is.null(f1) || is.null(f2)) return(NULL)
    ci <- function(f) {
      S <- cumsum(f$e)
      sum(S^2) / (f$n^2 * .lrv_est(f$e, lrv, bw))
    }
    Fst <- function(fu, fr) ((fr$ssr - fu$ssr) / 2) / (fu$ssr / (fu$n - fu$k))
    list(k = k, ssr_ols = f1$ssr, ssr_dols = f2$ssr,
         CI_ols = ci(f1), CI_dols = ci(f2),
         F_ols = Fst(f1, r1), F_dols = Fst(f2, r2), n_dols = f2$n)
  })
  res <- Filter(Negate(is.null), res)
  if (!length(res)) return(list(error = "Tsong regressions could not be estimated"))
  i1 <- which.min(vapply(res, `[[`, numeric(1), "ssr_ols"))
  i2 <- which.min(vapply(res, `[[`, numeric(1), "ssr_dols"))
  cv <- .fcoint_cv_tsong(n_x, res[[i2]]$k, model)

  list(test = "tsong",
       CI_stat = res[[i2]]$CI_dols, F_stat = res[[i2]]$F_dols,
       frequency = res[[i2]]$k, nobs = res[[i2]]$n_dols,
       CI_ols = res[[i1]]$CI_ols, F_ols = res[[i1]]$F_ols,
       frequency_ols = res[[i1]]$k,
       dolslags = q, lrv = lrv, bandwidth = bw,
       ci_cv1 = cv[1L], ci_cv5 = cv[2L], ci_cv10 = cv[3L],
       model = model)
}


# ------------------------------------------------------------------
# Critical values
# ------------------------------------------------------------------

# internal
.interp_T <- function(T_obs, grid, vals) {
  ## vals: matrix length(grid) x 3 (1%, 5%, 10%)
  if (T_obs <= grid[1L]) return(vals[1L, ])
  if (T_obs >= grid[length(grid)]) return(vals[length(grid), ])
  j <- max(which(grid <= T_obs))
  w <- (T_obs - grid[j]) / (grid[j + 1L] - grid[j])
  vals[j, ] + w * (vals[j + 1L, ] - vals[j, ])
}

# internal
.fcoint_cv_fadl <- function(n_x, k, T_obs, model) {
  ## Stata module fcoint (_fcoint_cv_fadl_single), attributed to
  ## Banerjee, Arcabic and Lee (2017), Tables 1a and 1b.
  ## Columns: 1%,5%,10% at T = 100 | 500 | 2500.
  if (n_x > 3L || k > 5L) return(rep(NA_real_, 3L))
  tab <- list(
    constant = rbind(
      c(-4.73, -4.09, -3.76, -4.61, -4.03, -3.72, -4.60, -4.02, -3.72),
      c(-4.44, -3.75, -3.37, -4.33, -3.70, -3.36, -4.33, -3.70, -3.35),
      c(-4.21, -3.51, -3.14, -4.13, -3.48, -3.14, -4.15, -3.48, -3.14),
      c(-4.07, -3.38, -3.03, -4.01, -3.37, -3.04, -4.04, -3.38, -3.04),
      c(-4.00, -3.32, -2.97, -3.94, -3.32, -2.99, -3.96, -3.32, -3.00),
      c(-4.96, -4.32, -3.98, -4.85, -4.26, -3.95, -4.81, -4.25, -3.95),
      c(-4.79, -4.10, -3.73, -4.66, -4.03, -3.70, -4.64, -4.03, -3.70),
      c(-4.56, -3.87, -3.49, -4.48, -3.83, -3.49, -4.48, -3.83, -3.49),
      c(-4.43, -3.73, -3.36, -4.34, -3.70, -3.36, -4.36, -3.71, -3.38),
      c(-4.35, -3.65, -3.28, -4.27, -3.64, -3.30, -4.28, -3.65, -3.31),
      c(-5.17, -4.51, -4.17, -5.06, -4.46, -4.15, -5.02, -4.45, -4.14),
      c(-5.04, -4.36, -3.99, -4.94, -4.30, -3.96, -4.90, -4.29, -3.96),
      c(-4.90, -4.16, -3.79, -4.76, -4.13, -3.78, -4.75, -4.13, -3.78),
      c(-4.75, -4.03, -3.65, -4.64, -4.00, -3.65, -4.65, -4.01, -3.67),
      c(-4.66, -3.94, -3.57, -4.58, -3.94, -3.58, -4.56, -3.93, -3.60)),
    trend = rbind(
      c(-5.17, -4.55, -4.24, -5.04, -4.47, -4.19, -5.00, -4.46, -4.18),
      c(-5.01, -4.34, -4.00, -4.86, -4.27, -3.96, -4.86, -4.27, -3.95),
      c(-4.79, -4.11, -3.76, -4.70, -4.07, -3.75, -4.68, -4.07, -3.74),
      c(-4.64, -3.96, -3.61, -4.58, -3.94, -3.61, -4.55, -3.95, -3.61),
      c(-4.53, -3.87, -3.52, -4.48, -3.86, -3.54, -4.46, -3.86, -3.54),
      c(-5.36, -4.72, -4.40, -5.23, -4.66, -4.37, -5.20, -4.65, -4.37),
      c(-5.22, -4.57, -4.23, -5.10, -4.51, -4.20, -5.07, -4.49, -4.19),
      c(-5.07, -4.39, -4.02, -4.95, -4.34, -4.01, -4.91, -4.32, -4.00),
      c(-4.94, -4.23, -3.87, -4.84, -4.22, -3.88, -4.81, -4.20, -3.87),
      c(-4.86, -4.15, -3.79, -4.76, -4.13, -3.80, -4.74, -4.13, -3.80),
      c(-5.54, -4.89, -4.55, -5.39, -4.83, -4.53, -5.39, -4.83, -4.53),
      c(-5.47, -4.77, -4.42, -5.30, -4.73, -4.40, -5.28, -4.70, -4.39),
      c(-5.33, -4.63, -4.27, -5.18, -4.57, -4.25, -5.18, -4.55, -4.23),
      c(-5.19, -4.49, -4.11, -5.08, -4.46, -4.12, -5.06, -4.44, -4.11),
      c(-5.08, -4.38, -4.01, -5.01, -4.36, -4.03, -4.99, -4.37, -4.03)))
  r <- tab[[model]][(n_x - 1L) * 5L + k, ]
  .interp_T(T_obs, c(100, 500, 2500), matrix(r, 3L, byrow = TRUE))
}

# internal
.fcoint_cv_feg <- function(n_x, k, T_obs, model) {
  ## Yilanci (2019), Table 1. Columns: 1%,5%,10% at T = 100 | 500 | 1000.
  if (n_x > 3L || k > 5L) return(rep(NA_real_, 3L))
  tab <- list(
    constant = rbind(
      c(-4.906, -4.302, -3.988, -4.756, -4.198, -3.898, -4.738, -4.175, -3.886),
      c(-4.665, -3.995, -3.648, -4.517, -3.912, -3.589, -4.503, -3.898, -3.579),
      c(-4.437, -3.743, -3.380, -4.333, -3.685, -3.349, -4.314, -3.686, -3.342),
      c(-4.285, -3.599, -3.252, -4.183, -3.554, -3.231, -4.172, -3.546, -3.221),
      c(-4.190, -3.520, -3.187, -4.091, -3.478, -3.165, -4.081, -3.477, -3.165),
      c(-5.282, -4.655, -4.337, -5.067, -4.511, -4.220, -5.048, -4.487, -4.205),
      c(-5.168, -4.526, -4.189, -4.969, -4.394, -4.085, -4.949, -4.371, -4.065),
      c(-4.958, -4.283, -3.938, -4.804, -4.183, -3.870, -4.778, -4.172, -3.852),
      c(-4.805, -4.122, -3.767, -4.647, -4.048, -3.722, -4.657, -4.040, -3.716),
      c(-4.708, -4.033, -3.689, -4.587, -3.964, -3.633, -4.536, -3.935, -3.629),
      c(-5.596, -4.957, -4.640, -5.354, -4.796, -4.512, -5.315, -4.786, -4.497),
      c(-5.573, -4.918, -4.593, -5.330, -4.752, -4.460, -5.286, -4.727, -4.435),
      c(-5.393, -4.733, -4.394, -5.177, -4.597, -4.285, -5.150, -4.582, -4.277),
      ## T = 1000: printed as -5.035 -4.134 -4.455 (5% and 10% transposed)
      c(-5.271, -4.605, -4.252, -5.071, -4.468, -4.148, -5.035, -4.455, -4.134),
      c(-5.155, -4.478, -4.127, -4.976, -4.378, -4.056, -4.959, -4.352, -4.042)),
    trend = rbind(
      c(-5.354, -4.731, -4.423, -5.128, -4.576, -4.293, -5.074, -4.555, -4.274),
      c(-5.243, -4.582, -4.250, -4.995, -4.433, -4.136, -4.973, -4.410, -4.119),
      c(-5.002, -4.340, -3.997, -4.801, -4.230, -3.910, -4.804, -4.208, -3.901),
      c(-4.849, -4.175, -3.827, -4.697, -4.092, -3.767, -4.693, -4.088, -3.769),
      c(-4.774, -4.086, -3.739, -4.634, -3.997, -3.683, -4.593, -3.994, -3.677),
      c(-5.641, -5.026, -4.705, -5.404, -4.855, -4.571, -5.367, -4.826, -4.550),
      c(-5.598, -4.954, -4.633, -5.329, -4.772, -4.480, -5.295, -4.748, -4.460),
      ## T = 500, 5%: printed as 4.620 (sign lost)
      c(-5.450, -4.781, -4.436, -5.199, -4.620, -4.313, -5.167, -4.597, -4.292),
      c(-5.294, -4.622, -4.271, -5.089, -4.487, -4.183, -5.065, -4.469, -4.158),
      c(-5.203, -4.508, -4.164, -5.006, -4.404, -4.086, -4.945, -4.370, -4.063),
      c(-5.941, -5.294, -4.971, -5.638, -5.094, -4.814, -5.602, -5.070, -4.795),
      c(-5.926, -5.278, -4.961, -5.635, -5.078, -4.791, -5.590, -5.048, -4.762),
      c(-5.792, -5.141, -4.806, -5.515, -4.964, -4.659, -5.504, -4.940, -4.643),
      c(-5.698, -5.023, -4.681, -5.441, -4.843, -4.534, -5.404, -4.835, -4.529),
      c(-5.601, -4.905, -4.560, -5.361, -4.752, -4.436, -5.332, -4.743, -4.435)))
  r <- tab[[model]][(n_x - 1L) * 5L + k, ]
  .interp_T(T_obs, c(100, 500, 1000), matrix(r, 3L, byrow = TRUE))
}

# internal
.fcoint_cv_tsong <- function(n_x, k, model) {
  ## TSPDLIB (Nazlioglu), cv_coint_tsongetal: rows p = 1..4, columns
  ## 10%, 5%, 1%. Returned as 1%, 5%, 10%.
  if (n_x > 4L || k > 3L) return(rep(NA_real_, 3L))
  tab <- list(
    constant = list(
      rbind(c(0.095, 0.124, 0.198), c(0.070, 0.092, 0.155),
            c(0.059, 0.076, 0.130), c(0.050, 0.061, 0.096)),
      rbind(c(0.200, 0.276, 0.473), c(0.132, 0.182, 0.328),
            c(0.098, 0.132, 0.215), c(0.072, 0.097, 0.171)),
      rbind(c(0.224, 0.304, 0.507), c(0.148, 0.202, 0.383),
            c(0.112, 0.146, 0.250), c(NA, NA, NA))),
    trend = list(
      rbind(c(0.042, 0.048, 0.063), c(0.038, 0.045, 0.059),
            c(0.036, 0.042, 0.055), c(0.034, 0.038, 0.050)),
      rbind(c(0.078, 0.099, 0.163), c(0.063, 0.081, 0.127),
            c(0.051, 0.066, 0.103), c(0.044, 0.055, 0.086)),
      rbind(c(0.090, 0.114, 0.170), c(0.075, 0.094, 0.143),
            c(0.061, 0.075, 0.116), c(0.053, 0.065, 0.099))))
  rev(tab[[model]][[k]][n_x, ])
}


#' Print Method for fcoint Objects
#'
#' @param x An object of class \code{"fcoint"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.fcoint <- function(x, ...) {
  cat("Fourier Cointegration Tests\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Model        : %s\n", x$model))
  cat(sprintf("Observations : %d\n", x$nobs))
  cat("\n")
  fmt <- function(v) ifelse(is.na(v), "   NA  ", sprintf("%7.3f", v))

  for (nm in names(x$results)) {
    r <- x$results[[nm]]
    if (!is.null(r$error)) {
      cat(sprintf("  [%s] Error: %s\n", toupper(nm), r$error))
      next
    }
    cat(sprintf("--- %s ---\n", toupper(nm)))
    if (nm == "tsong") {
      cat("  H0: cointegration (reject for large values)\n")
      cat(sprintf("  CI (DOLS)    : %8.4f   k* = %d, leads/lags = %d\n",
                  r$CI_stat, r$frequency, r$dolslags))
      cat(sprintf("  CI (OLS)     : %8.4f   k* = %d\n", r$CI_ols, r$frequency_ols))
      cat(sprintf("  F (DOLS)     : %8.4f   F (OLS): %8.4f\n", r$F_stat, r$F_ols))
      cat(sprintf("  CV 1%%/5%%/10%% : %s / %s / %s\n",
                  fmt(r$ci_cv1), fmt(r$ci_cv5), fmt(r$ci_cv10)))
      dec <- if (anyNA(c(r$ci_cv1, r$ci_cv5, r$ci_cv10))) "critical values unavailable" else
        if (r$CI_stat > r$ci_cv1) "Reject H0 at 1%" else
        if (r$CI_stat > r$ci_cv5) "Reject H0 at 5%" else
        if (r$CI_stat > r$ci_cv10) "Reject H0 at 10%" else "Do not reject H0"
    } else {
      cat("  H0: no cointegration (reject for small values)\n")
      cat(sprintf("  t-statistic  : %8.4f   k* = %d, lags = %s\n",
                  r$tstat, r$frequency,
                  paste(c(r$lag, r$lag_dx), collapse = "/")))
      if (!is.null(r$rho2))
        cat(sprintf("  rho^2        : %8.4f\n", r$rho2))
      cat(sprintf("  CV 1%%/5%%/10%% : %s / %s / %s\n",
                  fmt(r$cv1), fmt(r$cv5), fmt(r$cv10)))
      dec <- if (anyNA(c(r$cv1, r$cv5, r$cv10))) "critical values unavailable" else
        if (r$tstat < r$cv1) "Reject H0 at 1%" else
        if (r$tstat < r$cv5) "Reject H0 at 5%" else
        if (r$tstat < r$cv10) "Reject H0 at 10%" else "Do not reject H0"
    }
    cat(sprintf("  Decision     : %s\n\n", dec))
  }
  invisible(x)
}
