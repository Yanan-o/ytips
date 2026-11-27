test_that("piecewiseSEM formula and direct fit paths work", {
  skip_if_not_installed("piecewiseSEM")
  skip_if_not_installed("DiagrammeR")

  set.seed(1)
  n <- 60
  df <- data.frame(
    block = rep(1:6, each = 10),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
  df$y <- 0.4 * df$x1 + 0.2 * df$x2 + stats::rnorm(n)
  fit <- lmerTest::lmer(y ~ x1 + x2 + (1 | block), data = df)

  formula <- suppressWarnings(
    y_psem_formula(
      fit,
      sem.data = df,
      x.level1 = "x1",
      x.level2 = "x2"
    )
  )
  expect_match(formula$sem.formula, "piecewiseSEM::psem")

  sem.fit <- suppressWarnings(
    y_psem_boot(formula, nrep = 0, output = "sem.fit")
  )
  expect_s3_class(sem.fit, "psem")

  selected <- suppressWarnings(ytips:::yy_fit_step(fit, step = "step_aicc", output = "best_fit"))
  expect_true(inherits(selected, "merMod"))
})
