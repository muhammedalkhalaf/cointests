## Reference values computed with the Stata modules fcoint and xtcadfcoint
## (SSC, 2026 releases) on the data generated below.

.ref_data <- function() {
  set.seed(20260928)
  N <- 10; T <- 50
  pd <- do.call(rbind, lapply(1:N, function(i) {
    x1 <- cumsum(rnorm(T)); x2 <- cumsum(rnorm(T))
    u <- as.numeric(arima.sim(list(ar = 0.5), T))
    y <- 1 + i / 10 + (0.5 + i / 50) * x1 - 0.3 * x2 + u
    data.frame(id = i, t = 1:T, y = y, x1 = x1, x2 = x2)
  }))
  T <- 80
  x1 <- cumsum(rnorm(T)); x2 <- cumsum(rnorm(T))
  brk <- c(rep(0, 40), rep(2, 40))
  y <- 1 + brk + 0.8 * x1 + 0.4 * x2 + as.numeric(arima.sim(list(ar = 0.3), T))
  list(panel = pd, ts = data.frame(y = y, x1 = x1, x2 = x2))
}

test_that("fcoint FADL, FEG and FEG2 reproduce Stata fcoint", {
  d <- .ref_data()$ts
  r <- fcoint(d$y, d$x1, test = "all", model = "constant")$results
  expect_equal(r$fadl$tstat, -5.058970, tolerance = 1e-6)
  expect_equal(r$feg$tstat,  -7.843630, tolerance = 1e-6)
  expect_equal(r$feg2$tstat, -7.958511, tolerance = 1e-6)
  r2 <- fcoint(d$y, cbind(d$x1, d$x2), test = "all", model = "trend")$results
  expect_equal(r2$fadl$tstat, -5.054618, tolerance = 1e-6)
  expect_equal(r2$feg$tstat,  -4.272204, tolerance = 1e-6)
  expect_equal(r2$feg$frequency, 2L)
  expect_equal(r2$feg2$tstat, -4.481756, tolerance = 1e-6)
})

test_that("fcoint critical values follow the published tables", {
  d <- .ref_data()$ts
  r <- fcoint(d$y, d$x1, test = "feg", max_freq = 1)$results$feg
  ## Yilanci (2019), Table 1, n = 1, k = 1, constant: T = 100 is the lower bound
  expect_equal(c(r$cv1, r$cv5, r$cv10), c(-4.906, -4.302, -3.988))
  r <- fcoint(d$y, d$x1, test = "feg2", max_freq = 1)$results$feg2
  expect_true(all(is.na(c(r$cv1, r$cv5, r$cv10))))
})

test_that("xtcadfcoint reproduces Stata xtcadfcoint", {
  p <- .ref_data()$panel
  r <- xtcadfcoint(y ~ x1 + x2, data = p, index = c("id", "t"), model = 1)
  expect_equal(r$panel_cips, -4.659360, tolerance = 1e-6)
  expect_equal(r$beta_ccep[1], 0.594208, tolerance = 1e-5)
  r <- xtcadfcoint(y ~ x1, data = p, index = c("id", "t"), model = 3,
                   breaks = 1, lagselect = "aic")
  expect_equal(r$panel_cips, -4.029347, tolerance = 1e-6)
  expect_equal(r$panel_cips_alt, -4.165257, tolerance = 1e-6)
  expect_equal(c(r$Tb_hat, r$Tb_tilde), c(36L, 14L))
  r <- xtcadfcoint(y ~ x1, data = p, index = c("id", "t"), model = 5,
                   breaks = 1, brk_slope = TRUE)
  expect_equal(r$panel_cips, -4.844233, tolerance = 1e-6)
  expect_equal(r$panel_cips_trim, -2.859926, tolerance = 1e-6)
  expect_equal(r$Tb_trim, 29L)
})
