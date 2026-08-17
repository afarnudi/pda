get_logL_D1 <- pda:::get_logL_D1
get_logL_D2 <- pda:::get_logL_D2

make_surv_data <- function(n = 200, seed = 1, p = 1) {
  set.seed(seed)
  df <- data.frame(
    time = rexp(n, rate = 0.1),
    status = rbinom(n, 1, 0.7)
  )
  for (i in seq_len(p)) {
    df[[paste0("X", i)]] <- rnorm(n)
  }
  df
}

fit_coxph <- function(p, ...) {
  df <- make_surv_data(p = p)
  covars <- paste(paste0("X", seq_len(p)), collapse = "+")
  form <- as.formula(paste("Surv(time, status) ~", covars))
  survival::coxph(form, data = df, ...)
}

test_that("get_logL_D1 returns a vector of the right length for a single covariate", {
  fit1 <- fit_coxph(p = 1)
  d1 <- get_logL_D1(fit1)

  expect_type(d1, "double")
  expect_length(d1, 1)
  # score residuals sum to ~0 at the MLE
  expect_equal(d1, 0, tolerance = 1e-6)
})

test_that("get_logL_D1 returns a vector of the right length for multiple covariates", {
  fit3 <- fit_coxph(p = 3)
  d1 <- get_logL_D1(fit3)

  expect_type(d1, "double")
  expect_length(d1, 3)
  expect_equal(d1, rep(0, 3), tolerance = 1e-6)
})

test_that("get_logL_D2 returns a 1x1 matrix for a single covariate and matches the model variance", {
  fit1 <- fit_coxph(p = 1)
  d2 <- get_logL_D2(fit1)

  expect_true(is.matrix(d2))
  expect_equal(dim(d2), c(1L, 1L))
  # observed information is the inverse of the naive variance at the MLE
  expect_equal(solve(-d2), fit1$var, tolerance = 1e-6)
})

test_that("get_logL_D2 returns a pxp matrix for multiple covariates and matches the model variance", {
  fit3 <- fit_coxph(p = 3)
  d2 <- get_logL_D2(fit3)

  expect_true(is.matrix(d2))
  expect_equal(dim(d2), c(3L, 3L))
  expect_equal(solve(-d2), fit3$var, tolerance = 1e-6)
})

test_that("get_logL_D1/D2 work off the MLE (as used in ODACH_CC.derive/estimate) for p = 1", {
  df <- make_surv_data(p = 1)
  fit0 <- survival::coxph(
    Surv(time, status) ~ X1, data = df,
    init = 0, method = "breslow",
    control = survival::coxph.control(iter.max = 0)
  )

  d1 <- expect_no_error(get_logL_D1(fit0))
  d2 <- expect_no_error(get_logL_D2(fit0))

  expect_length(d1, 1)
  expect_equal(dim(d2), c(1L, 1L))
})

test_that("get_logL_D1/D2 work off the MLE for multiple covariates", {
  df <- make_surv_data(p = 3)
  fit0 <- survival::coxph(
    Surv(time, status) ~ X1 + X2 + X3, data = df,
    init = c(0, 0, 0), method = "breslow",
    control = survival::coxph.control(iter.max = 0)
  )

  d1 <- expect_no_error(get_logL_D1(fit0))
  d2 <- expect_no_error(get_logL_D2(fit0))

  expect_length(d1, 3)
  expect_equal(dim(d2), c(3L, 3L))
})
