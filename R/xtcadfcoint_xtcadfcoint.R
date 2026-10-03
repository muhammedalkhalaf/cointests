#' Panel CADF Cointegration Test with Structural Breaks
#'
#' Tests the null hypothesis of no cointegration in panel data using the
#' cross-sectionally augmented Dickey-Fuller (CADF) approach of Banerjee and
#' Carrion-i-Silvestre (2025). Accounts for cross-sectional dependence via the
#' Common Correlated Effects (CCE) estimator and allows for structural breaks.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying the
#'   cointegrating relationship to test.
#' @param data A data frame containing the panel data in long format.
#' @param index A character vector of length 2: \code{c("id_var", "time_var")}.
#' @param model Integer (0--5) specifying the deterministic component:
#'   \itemize{
#'     \item 0: No deterministic component
#'     \item 1: Constant (default)
#'     \item 2: Linear trend
#'     \item 3: Constant with level shifts (requires \code{breaks >= 1})
#'     \item 4: Linear trend with level shifts (requires \code{breaks >= 1})
#'     \item 5: Linear trend with level and slope shifts (requires \code{breaks >= 1})
#'   }
#' @param breaks Integer (0, 1, or 2). Number of structural breaks. Default is 0.
#' @param trimming Numeric trimming fraction for break date search. Default 0.15.
#' @param maxlags Maximum lag order for ADF augmentation. Default 4.
#' @param lagselect Lag selection method: \code{"bic"} (default), \code{"aic"},
#'   \code{"maic"}, \code{"mbic"}, or \code{"fixed"}.
#' @param nfactors Integer. Number of common factors for CCE. Default 1.
#' @param brk_slope Logical. If \code{TRUE}, allows breaks in the cointegrating
#'   vector slopes. Default \code{FALSE}.
#' @param brk_loadings Logical. If \code{TRUE}, allows breaks in factor loadings.
#'   Default \code{FALSE}.
#' @param cce Logical. If \code{TRUE} (default), applies CCE cross-sectional
#'   augmentation to account for common factors.
#' @param simulate Integer. Number of Monte Carlo replications used to
#'   simulate critical values from independent Gaussian random walks at the
#'   estimated break dates (see Details); this is not a bootstrap of the
#'   data. Use 0 (default) to skip simulation.
#' @param level Confidence level (in percent) for hypothesis test decisions.
#'   Default 95.
#'
#' @return An object of class \code{"xtcadfcoint"} with components:
#'   \describe{
#'     \item{panel_cips}{Panel CIPS statistic; with breaks, computed at the
#'       break dates \code{Tb_hat}.}
#'     \item{panel_cips_alt}{Panel statistic at the alternative break dates
#'       \code{Tb_tilde} (\code{NULL} without breaks).}
#'     \item{panel_cips_trim}{Panel statistic on the Kim and Perron (2009)
#'       trimmed data (model 5 only).}
#'     \item{t_individual, t_trim}{Individual CADF (or ADF) t-statistics.}
#'     \item{p_selected}{Selected lag orders per unit.}
#'     \item{beta_ccep, beta_ccep_alt}{Pooled CCE estimates of the
#'       cointegrating vector (regime-specific when \code{brk_slope = TRUE}).}
#'     \item{SSR}{Sums of squared residuals: levels and first differences of
#'       the defactored and of the detrended residuals.}
#'     \item{Tb_hat, Tb_tilde, Tb_trim}{Estimated break dates (positions in
#'       the sample; a break at \eqn{T_b} shifts the level from
#'       \eqn{T_b + 1}).}
#'     \item{cv}{Simulated critical values (if \code{simulate > 0}): rows
#'       \code{panel} and \code{individual}, columns 1\%, 2.5\%, 5\%,
#'       10\%.}
#'     \item{N, TT, k, model, breaks, ...}{Settings of the call.}
#'   }
#'
#' @details
#' The pooled CCE estimator (equation 13 of the paper) is computed with the
#' projection on the deterministic terms and the cross-section averages of
#' \eqn{y} and \eqn{x} (interacted with the break dummies when
#' \code{brk_loadings = TRUE}). With breaks, the break dates are chosen by
#' a grid search: \code{Tb_hat} minimises the sum of squared defactored
#' residuals (equation 14) and \code{Tb_tilde} the sum of squared residuals
#' after removing only the deterministic terms (equations 15 and 16). Each
#' unit's statistic is the t-ratio of the lagged residual in the
#' cross-section augmented ADF regression (equation 17), with impulse
#' dummies at \eqn{T_b + 1}; \code{nfactors} sets how many
#' cross-section averages enter. For model 5 the statistic is also computed
#' on data trimmed by three observations on each side of the break (Kim
#' and Perron, 2009), with \code{maxlags} fixed lags.
#'
#' The computations reproduce the Stata module \code{xtcadfcoint}, which
#' states that it translates the authors' GAUSS code; results were checked
#' against it. The information criteria and the residual variance use the
#' full \eqn{T} as in that code, and the MAIC and MBIC penalties follow Ng
#' and Perron (2001).
#'
#' Critical values depend on \eqn{N}, \eqn{T}, \eqn{k}, the model and
#' the break fractions (Tables B.1 to B.24 of the supplementary material
#' of the paper cover one break). \code{simulate} draws them from
#' independent random walks at the estimated break dates \code{Tb_hat},
#' without lag augmentation. The Stata module instead fixes the break
#' fractions at 0.5 (one break) or 0.3 and 0.7 (two breaks).
#'
#' @references
#' Banerjee, A. and Carrion-i-Silvestre, J.L. (2025).
#' Panel Data Cointegration Testing with Structural Instabilities.
#' \emph{Journal of Business & Economic Statistics}, 43(1), 122--133.
#' \doi{10.1080/07350015.2024.2327844}
#'
#' Kim, D. and Perron, P. (2009). Unit root tests allowing for a break in
#' the trend function at an unknown time under both the null and
#' alternative hypotheses. \emph{Journal of Econometrics}, 148(1), 1--13.
#' \doi{10.1016/j.jeconom.2008.08.019}
#'
#' Ng, S. and Perron, P. (2001). Lag length selection and the construction
#' of unit root tests with good size and power. \emph{Econometrica}, 69(6),
#' 1519--1554. \doi{10.1111/1468-0262.00256}
#'
#' Pesaran, M.H. (2006). Estimation and Inference in Large Heterogeneous Panels
#' with a Multifactor Error Structure. \emph{Econometrica}, 74(4), 967--1012.
#' \doi{10.1111/j.1468-0262.2006.00692.x}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 8; tt <- 30
#' dat <- data.frame(
#'   id   = rep(1:n, each = tt),
#'   time = rep(1:tt, times = n),
#'   y    = cumsum(rnorm(n * tt)),
#'   x1   = cumsum(rnorm(n * tt))
#' )
#' res <- xtcadfcoint(y ~ x1, data = dat, index = c("id", "time"),
#'                    model = 1, breaks = 0)
#' print(res)
#' summary(res)
#' }
#'
#' @export
xtcadfcoint <- function(formula, data, index,
                         model = 1L,
                         breaks = 0L,
                         trimming = 0.15,
                         maxlags = 4L,
                         lagselect = "bic",
                         nfactors = 1L,
                         brk_slope = FALSE,
                         brk_loadings = FALSE,
                         cce = TRUE,
                         simulate = 0L,
                         level = 95L) {

  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula object.", call. = FALSE)
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (!is.character(index) || length(index) != 2)
    stop("'index' must be a character vector of length 2.", call. = FALSE)
  if (!all(index %in% names(data)))
    stop("Variables in 'index' not found in 'data'.", call. = FALSE)
  model    <- as.integer(model)
  breaks   <- as.integer(breaks)
  maxlags  <- as.integer(maxlags)
  nfactors <- as.integer(nfactors)
  simulate <- as.integer(simulate)
  if (model < 0L || model > 5L)
    stop("'model' must be an integer between 0 and 5.", call. = FALSE)
  if (breaks < 0L || breaks > 2L)
    stop("'breaks' must be 0, 1, or 2.", call. = FALSE)
  if (breaks == 0L && model >= 3L)
    stop("Models 3-5 require breaks >= 1.", call. = FALSE)
  if (breaks > 0L && model < 3L)
    stop("breaks > 0 requires model >= 3.", call. = FALSE)
  if (maxlags < 0L) stop("'maxlags' must be non-negative.", call. = FALSE)
  if (nfactors < 1L) stop("'nfactors' must be at least 1.", call. = FALSE)
  lagselect <- match.arg(lagselect, c("bic", "aic", "maic", "mbic", "fixed"))

  ivar <- index[1]
  tvar <- index[2]
  mf      <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  depvar  <- names(mf)[1]
  indvars <- names(mf)[-1]
  k       <- length(indvars)
  if (k < 1) stop("At least one independent variable required.", call. = FALSE)

  keep   <- stats::complete.cases(data[, c(depvar, indvars, ivar, tvar), drop = FALSE])
  data_c <- data[keep, , drop = FALSE]
  data_c <- data_c[order(data_c[[ivar]], data_c[[tvar]]), , drop = FALSE]
  panels <- sort(unique(data_c[[ivar]]))
  times  <- sort(unique(data_c[[tvar]]))
  N  <- length(panels)
  TT <- length(times)
  if (N < 2) stop("At least 2 panel units are required.", call. = FALSE)
  if (TT < 10) stop("At least 10 time periods are required.", call. = FALSE)
  if (nrow(data_c) != N * TT)
    stop("Panel must be balanced for xtcadfcoint.", call. = FALSE)
  if (nfactors > k + 1L) {
    message("nfactors > k+1; setting nfactors = ", k + 1L)
    nfactors <- k + 1L
  }

  ## Y: TT x N; X: TT x (N*k), column i + (j-1)*N holds regressor j of unit i
  Y <- matrix(NA_real_, TT, N)
  X <- matrix(NA_real_, TT, N * k)
  pid <- match(data_c[[ivar]], panels)
  tid <- match(data_c[[tvar]], times)
  Y[cbind(tid, pid)] <- data_c[[depvar]]
  for (j in seq_len(k)) X[cbind(tid, pid + (j - 1L) * N)] <- data_c[[indvars[j]]]

  opt <- list(model = model, brk_slope = isTRUE(brk_slope),
              brk_loadings = isTRUE(brk_loadings), nf = nfactors,
              p_max = maxlags, auto = lagselect != "fixed",
              ic = switch(lagselect, aic = 0L, bic = 1L, maic = 2L,
                          mbic = 3L, fixed = 1L),
              cce = isTRUE(cce))

  if (breaks == 0L) {
    r <- .bcs_main(Y, X, 0L, opt)
    est <- list(main = r, alt = r, Tb = NULL, Tb_alt = NULL, trim = NULL)
  } else {
    est <- .bcs_endog(Y, X, breaks, trimming, opt)
  }

  cv_mat <- NULL
  if (simulate > 0L) {
    message("Simulating critical values (", simulate, " replications)...")
    cv_mat <- .bcs_simulate_cv(N, TT, k, if (breaks > 0L) est$Tb else 0L,
                               opt, simulate)
  }

  out <- list(
    panel_cips      = est$main$panel,
    panel_cips_alt  = if (breaks > 0L) est$alt$panel else NULL,
    panel_cips_trim = if (!is.null(est$trim)) est$trim$panel else NULL,
    t_individual    = est$main$t,
    t_trim          = if (!is.null(est$trim)) est$trim$t else NULL,
    p_selected      = est$main$p,
    beta_ccep       = est$main$beta,
    beta_ccep_alt   = if (breaks > 0L) est$alt$beta else NULL,
    SSR             = est$main$ssr,
    Tb_hat          = est$Tb,
    Tb_tilde        = est$Tb_alt,
    Tb_trim         = if (!is.null(est$trim)) est$trim$Tb else NULL,
    N = N, TT = TT, k = k, model = model, breaks = breaks,
    trimming = trimming, lagselect = lagselect, maxlags = maxlags,
    nfactors = nfactors, brk_slope = opt$brk_slope,
    brk_loadings = opt$brk_loadings, cce = opt$cce,
    depvar = depvar, indepvars = indvars, panels = panels, times = times,
    cv = cv_mat, level = level)
  class(out) <- "xtcadfcoint"
  out
}


## ------------------------------------------------------------------
## Engine: port of the Mata translation (xtcadfcoint_mata.ado) of the
## authors' GAUSS procedures CADFcoin_multiple, cadf_multiple,
## adf_multiple, CADFcoin_multiple_endog and KimPerron_trimdata.
## ------------------------------------------------------------------

# internal
.bcs_ginv <- function(A) {
  out <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(out)) {
    s <- svd(A)
    d <- ifelse(s$d > max(s$d) * 1e-12, 1 / s$d, 0)
    out <- s$v %*% (d * t(s$u))
  }
  out
}

# internal
.bcs_lag <- function(M, n) {
  M <- as.matrix(M)
  if (n >= nrow(M)) return(matrix(NA_real_, nrow(M), ncol(M)))
  rbind(matrix(NA_real_, n, ncol(M)), M[seq_len(nrow(M) - n), , drop = FALSE])
}

# internal
.bcs_dummies <- function(TT, Tb) {
  n_br <- length(Tb)
  DU <- DT <- DTb <- matrix(0, TT, n_br)
  for (i in seq_len(n_br)) {
    DU[, i]  <- c(rep(0, Tb[i]), rep(1, TT - Tb[i]))
    DT[, i]  <- c(rep(0, Tb[i]), seq_len(TT - Tb[i]))
    DTb[Tb[i] + 1L, i] <- 1
  }
  list(DU = DU, DT = DT, DTb = DTb)
}

# internal
.bcs_deter <- function(TT, model, dm, impulse) {
  tr <- seq_len(TT)
  switch(as.character(model),
    "0" = NULL,
    "1" = matrix(1, TT, 1),
    "2" = cbind(1, tr),
    "3" = cbind(1, dm$DU, if (impulse) dm$DTb),
    "4" = cbind(1, tr, dm$DU, if (impulse) dm$DTb),
    "5" = cbind(1, tr, dm$DU, dm$DT, if (impulse) dm$DTb))
}

# internal
.bcs_unit_x <- function(X, i, N, k, brk_slope, DU) {
  xi <- X[, i + (seq_len(k) - 1L) * N, drop = FALSE]
  if (brk_slope && !is.null(DU) && ncol(DU) > 0L)
    xi <- cbind(xi, do.call(cbind, lapply(seq_len(ncol(DU)),
                                          function(b) DU[, b] * xi)))
  xi
}

# internal
.bcs_main <- function(Y, X, Tb, opt, model = opt$model, auto = opt$auto,
                      stats = TRUE) {
  TT <- nrow(Y); N <- ncol(Y); k <- ncol(X) / N
  if (length(Tb) == 1L && Tb[1] == 0) Tb <- integer(0)
  n_br <- length(Tb)
  dm <- .bcs_dummies(TT, Tb)
  DU <- if (n_br > 0L && (model >= 3L || opt$brk_slope || opt$brk_loadings)) dm$DU else NULL
  x_deter <- .bcs_deter(TT, model, dm, impulse = FALSE)

  xbar <- sapply(seq_len(k), function(j) rowMeans(X[, (j - 1L) * N + seq_len(N), drop = FALSE]))
  xbar <- matrix(xbar, TT, k)
  ybar <- rowMeans(Y)
  if (!opt$cce) {
    H <- x_deter
  } else {
    cam <- cbind(ybar, xbar)
    H <- cam
    if (opt$brk_loadings && n_br > 0L)
      for (b in seq_len(n_br)) H <- cbind(H, dm$DU[, b] * cam)
    H <- cbind(x_deter, H)
  }
  M <- if (is.null(H)) diag(TT) else diag(TT) - H %*% .bcs_ginv(crossprod(H)) %*% t(H)

  kk  <- if (opt$brk_slope) k * (n_br + 1L) else k
  den <- matrix(0, kk, kk); num <- matrix(0, kk, 1)
  xs  <- vector("list", N)
  for (i in seq_len(N)) {
    xi <- .bcs_unit_x(X, i, N, k, opt$brk_slope, DU)
    xs[[i]] <- xi
    Mx <- M %*% xi
    den <- den + crossprod(xi, Mx)
    num <- num + crossprod(Mx, Y[, i])
  }
  beta <- as.numeric(.bcs_ginv(den) %*% num)

  EG <- sapply(seq_len(N), function(i) Y[, i] - xs[[i]] %*% beta)
  EG <- matrix(EG, TT, N)
  EGdt <- if (is.null(x_deter)) EG else
    EG - x_deter %*% (.bcs_ginv(crossprod(x_deter)) %*% crossprod(x_deter, EG))
  EGdf <- M %*% EG
  d1 <- function(A) A[-1L, , drop = FALSE] - A[-TT, , drop = FALSE]
  ssr <- c(EG_defactored = sum(EGdf^2), EG_detrended = sum(EGdt^2),
           DEG_defactored = sum(d1(EGdf)^2), DEG_detrended = sum(d1(EGdt)^2))

  out <- list(beta = beta, ssr = ssr)
  if (!stats) return(out)
  tp <- .bcs_cadf(EG, if (opt$cce) xbar else NULL, model, Tb, opt, auto)
  c(out, list(t = tp$t, p = tp$p, panel = mean(tp$t)))
}

# internal
.bcs_cadf <- function(EG, xbar, model, Tb, opt, auto) {
  TT <- nrow(EG); N <- ncol(EG)
  n_br <- length(Tb)
  cce <- !is.null(xbar)
  dm <- .bcs_dummies(TT, Tb)
  x_deter <- .bcs_deter(TT, model, dm, impulse = TRUE)
  p_max <- opt$p_max

  DEG  <- EG[-1L, , drop = FALSE] - EG[-TT, , drop = FALSE]
  EGl  <- EG[-TT, , drop = FALSE]                  # EG_resid_lag[2::T, ]
  if (cce) {
    cm  <- rowMeans(EG)
    lagterms <- matrix(cm[-TT], ncol = 1)
    dterms   <- matrix(rowMeans(DEG), ncol = 1)
    nf <- opt$nf
    if (nf > 1L) {
      xa <- xbar[, seq_len(min(nf - 1L, ncol(xbar))), drop = FALSE]
      lagterms <- cbind(lagterms, xa[-TT, , drop = FALSE])
      dterms   <- cbind(dterms, xa[-1L, , drop = FALSE] - xa[-TT, , drop = FALSE])
    }
    if (opt$brk_loadings && n_br > 0L) {
      for (b in seq_len(n_br)) {
        lagterms <- cbind(lagterms, dm$DU[-1L, b] * cm[-TT])
        dterms   <- cbind(dterms, dm$DU[-1L, b] * rowMeans(DEG))
      }
    }
  }

  t_out <- rep(NA_real_, N); p_out <- integer(N)
  fit_t <- function(y, x) {
    xtx <- .bcs_ginv(crossprod(x))
    b <- xtx %*% crossprod(x, y)
    e <- y - x %*% b
    s2 <- sum(e^2) / (TT - ncol(x))
    list(b = b, s2 = s2, t1 = b[1] / sqrt(s2 * xtx[1, 1]), x1 = x[, 1])
  }
  for (i in seq_len(N)) {
    y0 <- DEG[, i]
    x0 <- cbind(EGl[, i], if (cce) cbind(lagterms, dterms),
                if (!is.null(x_deter)) x_deter[-1L, , drop = FALSE])
    make <- function(p) {
      if (p == 0L) return(list(y = y0, x = x0))
      lt <- do.call(cbind, lapply(seq_len(p), function(ii)
        cbind(.bcs_lag(y0, ii), if (cce) .bcs_lag(dterms, ii))))
      rows <- (p_max + 1L):length(y0)
      list(y = y0[rows], x = cbind(x0, lt)[rows, , drop = FALSE])
    }
    if (!auto) {
      d <- make(p_max)
      t_out[i] <- fit_t(d$y, d$x)$t1
      p_out[i] <- p_max
    } else {
      best <- Inf
      for (ip in p_max:0L) {
        d <- make(ip)
        f <- fit_t(d$y, d$x)
        nc <- ncol(d$x)
        tau <- if (opt$ic >= 2L) (f$b[1]^2) * sum(f$x1^2) / TT^2 / f$s2 else 0
        ic <- switch(as.character(opt$ic),
          "0" = log(f$s2) + 2 * nc / TT,
          "1" = log(f$s2) + log(TT) * nc / TT,
          "2" = log(f$s2) + 2 * (nc + tau) / TT,
          "3" = log(f$s2) + log(TT) * (nc + tau) / TT)
        if (ic < best) {
          best <- ic
          p_out[i] <- ip
          t_out[i] <- f$t1
        }
      }
    }
  }
  list(t = t_out, p = p_out)
}

# internal
.bcs_endog <- function(Y, X, m, trimming, opt) {
  TT <- nrow(Y); N <- ncol(Y)
  h <- floor(trimming * TT)
  if (m == 1L) {
    cand <- matrix((h + 1L):(TT - h), ncol = 1)
  } else {
    cand <- do.call(rbind, lapply((h + 1L):(TT - 2L * h), function(i)
      cbind(i, (h + i):(TT - h))))
  }
  ssr <- t(apply(cand, 1, function(tb) .bcs_main(Y, X, tb, opt, stats = FALSE)$ssr))
  Tb  <- as.integer(cand[which.min(ssr[, "EG_defactored"]), ])
  Tb2 <- as.integer(cand[which.min(ssr[, "EG_detrended"]), ])
  main <- .bcs_main(Y, X, Tb, opt)
  alt  <- .bcs_main(Y, X, Tb2, opt)

  trim <- NULL
  if (opt$model == 5L) {
    kp <- .bcs_kimperron_trim(cbind(Y, X), Tb, 3L)
    if (nrow(kp$y) > 0L) {
      Yt <- kp$y[, seq_len(N), drop = FALSE]
      Xt <- kp$y[, -seq_len(N), drop = FALSE]
      rt <- if (kp$Tb[1] == 0) {
        o2 <- opt; o2$brk_slope <- FALSE; o2$brk_loadings <- FALSE
        .bcs_main(Yt, Xt, 0L, o2, model = 2L, auto = FALSE)
      } else {
        .bcs_main(Yt, Xt, kp$Tb, opt, auto = FALSE)
      }
      trim <- list(panel = rt$panel, t = rt$t, Tb = kp$Tb)
    }
  }
  list(main = main, alt = alt, Tb = Tb, Tb_alt = Tb2, trim = trim)
}

# internal
.bcs_kimperron_trim <- function(y, Tb, trm) {
  TT <- nrow(y); n_br <- length(Tb)
  lo <- Tb - trm; hi <- Tb + trm
  fin <- cbind(ifelse(lo > 1, lo, NA), ifelse(hi < TT, hi, NA))
  filt <- fin[stats::complete.cases(fin), , drop = FALSE]
  low_ind <- if (any(is.na(fin[, 1]))) max(which(is.na(fin[, 1]))) else 0L
  high_ind <- if (all(is.na(fin[, 2]))) 1L else 0L
  disc <- rep(0, TT)
  if (low_ind > 0L && !is.na(fin[low_ind, 2])) disc[seq_len(fin[low_ind, 2])] <- 1
  if (high_ind > 0L && !is.na(fin[high_ind, 1])) disc[(fin[high_ind, 1] + 1):TT] <- 1
  marker <- rep(0, TT)
  yh <- y
  if (nrow(filt) > 0L) {
    S <- matrix(0, TT, ncol(y))
    for (r in seq_len(nrow(filt))) {
      a <- filt[r, 1]; b <- filt[r, 2]
      disc[(a + 1):b] <- 1
      marker[a] <- a
      step <- c(rep(0, a + 1), rep(1, TT - a - 1))
      S <- S + outer(step, y[b, ] - y[a, ])
    }
    yh <- y - S
  }
  keepr <- which(disc == 0)
  if (!length(keepr)) return(list(y = y[0, , drop = FALSE], Tb = 0L))
  tb <- which(marker[keepr] > 0)
  list(y = yh[keepr, , drop = FALSE], Tb = if (length(tb)) tb else 0L)
}

# internal
.bcs_simulate_cv <- function(N, TT, k, Tb, opt, reps) {
  o <- opt; o$p_max <- 0L; o$ic <- 0L
  burn <- 50L
  panel <- ind <- numeric(reps)
  for (r in seq_len(reps)) {
    Ys <- rbind(0, apply(matrix(stats::rnorm((TT + burn - 1L) * N), ncol = N), 2, cumsum))
    Xs <- rbind(0, apply(matrix(stats::rnorm((TT + burn - 1L) * N * k), ncol = N * k), 2, cumsum))
    Ys <- Ys[(burn + 1L):(TT + burn), , drop = FALSE]
    Xs <- Xs[(burn + 1L):(TT + burn), , drop = FALSE]
    res <- .bcs_main(Ys, Xs, Tb, o, auto = FALSE)
    panel[r] <- res$panel
    ind[r]   <- res$t[1]
  }
  q <- c(0.01, 0.025, 0.05, 0.10)
  pick <- function(v) sort(v)[pmax(1, floor(q * reps))]
  cv <- rbind(panel = pick(panel), individual = pick(ind))
  colnames(cv) <- c("1%", "2.5%", "5%", "10%")
  cv
}


## ── S3 methods ───────────────────────────────────────────────────────────

#' Print method for xtcadfcoint objects
#'
#' @param x An object of class \code{"xtcadfcoint"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtcadfcoint <- function(x, ...) {
  model_desc <- c(
    "0" = "No deterministic component",
    "1" = "Constant",
    "2" = "Linear trend",
    "3" = "Constant with level shifts",
    "4" = "Linear trend with level shifts",
    "5" = "Linear trend with level and slope shifts"
  )
  lag_desc <- c(bic = "Automatic (BIC)", aic = "Automatic (AIC)",
                maic = "Automatic (MAIC)", mbic = "Automatic (MBIC)",
                fixed = "Fixed")

  cat("\n")
  cat(strrep("-", 78), "\n")
  cat("  Banerjee & Carrion-i-Silvestre (2025)\n")
  cat("  Panel CADF Cointegration Test with Structural Breaks\n")
  cat(strrep("-", 78), "\n")
  cat("  H0: No cointegration\n")
  cat("  H1: Cointegration (panel is cointegrated)\n")
  cat(strrep("-", 78), "\n")
  cat(sprintf("  Model specification    : %s\n",
              model_desc[as.character(x$model)]))
  cat(sprintf("  Number of breaks (m)   : %d\n", x$breaks))
  cat(sprintf("  Trimming fraction      : %.2f\n", x$trimming))
  cat(sprintf("  Cross-section dep (CCE): %s\n",
              if (x$cce) "Yes (CCE)" else "No"))
  cat(sprintf("  Number of factors      : %d\n", x$nfactors))
  cat(sprintf("  Lag selection          : %s\n",
              lag_desc[x$lagselect]))
  cat(sprintf("  Panel dimensions       : N = %d, T = %d\n", x$N, x$TT))
  cat(sprintf("  Dep. variable          : %s\n", x$depvar))
  cat(sprintf("  Regressors (k = %d)    : %s\n", x$k,
              paste(x$indepvars, collapse = ", ")))
  cat(strrep("-", 78), "\n\n")

  cat("  Panel CIPS Cointegration Test Results\n")
  cat(strrep("-", 78), "\n")
  cat(sprintf("  CIPS statistic (lambda_hat)  : %12.4f\n", x$panel_cips))
  if (x$breaks > 0) {
    cat(sprintf("  CIPS statistic (lambda_tilde): %12.4f\n", x$panel_cips_alt))
    if (!is.null(x$panel_cips_trim))
      cat(sprintf("  CIPS statistic (KP trimmed)  : %12.4f\n", x$panel_cips_trim))
  }
  cat(strrep("-", 78), "\n")

  ## Critical values note
  cat("\n  Note: Critical values depend on N, T, k, model, and break specification.\n")
  cat("  Refer to Tables B.13-B.24 in Banerjee & Carrion-i-Silvestre (2025)\n")
  cat("  or use simulate argument for bootstrap critical values.\n\n")

  ## Simulated critical values if available
  if (!is.null(x$cv)) {
    cv <- x$cv["panel", ]
    cat("  Simulated critical values of the panel statistic (known break dates,\n")
    cat("  independent random walks, no lag augmentation):\n")
    cat(strrep("-", 50), "\n")
    for (nm in names(cv)) {
      cat(sprintf("  %5s: %10.4f  | %s\n", nm, cv[[nm]],
                  if (x$panel_cips < cv[[nm]]) "Reject H0" else "Do not reject H0"))
    }
    cat(strrep("-", 50), "\n\n")
  }

  ## Break dates
  if (x$breaks > 0 && !is.null(x$Tb_hat)) {
    cat("  Estimated Break Dates:\n")
    cat(strrep("-", 50), "\n")
    for (j in seq_along(x$Tb_hat)) {
      cat(sprintf("  Break %d: Tb_hat = %d", j, x$Tb_hat[j]))
      if (!is.null(x$Tb_tilde) && length(x$Tb_tilde) >= j) {
        cat(sprintf("  | Tb_tilde = %d", x$Tb_tilde[j]))
      }
      cat("\n")
    }
    cat(strrep("-", 50), "\n\n")
  }

  ## Pooled CCE beta
  cat("  Pooled CCE Estimator (beta_hat):\n")
  cat(strrep("-", 46), "\n")
  for (j in seq_along(x$beta_ccep)) {
    reg <- (j - 1L) %/% x$k
    nm <- x$indepvars[(j - 1L) %% x$k + 1L]
    if (reg > 0L) nm <- paste0(nm, " x DU", reg)
    cat(sprintf("  %-20s: %12.6f\n", nm, x$beta_ccep[j]))
  }
  cat(strrep("-", 46), "\n\n")

  ## Individual statistics (first few)
  cat("  Individual CADF/ADF Statistics:\n")
  cat(strrep("-", 56), "\n")
  cat(sprintf("  %-15s  %12s  %10s\n", "Unit", "t-statistic", "Lag"))
  cat(strrep("-", 56), "\n")
  n_show <- min(x$N, 20L)
  for (pi in seq_len(n_show)) {
    cat(sprintf("  %-15s  %12.4f  %10d\n",
                as.character(x$panels[pi]),
                x$t_individual[pi],
                x$p_selected[pi]))
  }
  if (x$N > 20L) cat(sprintf("  ... (%d more units)\n", x$N - 20L))
  cat(strrep("-", 56), "\n")

  cat("\n  Reference: Banerjee & Carrion-i-Silvestre (2025, JBES)\n")
  cat(strrep("-", 78), "\n\n")

  invisible(x)
}

#' Summary method for xtcadfcoint objects
#'
#' @param object An object of class \code{"xtcadfcoint"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtcadfcoint <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
