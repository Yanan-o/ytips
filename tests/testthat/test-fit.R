test_that("y_fit_summary can suppress diagnostics", {
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)

  expect_true("diagnostics" %in% names(formals(y_fit_summary)))
  result <- y_fit_summary(fit, output = "all", diagnostics = FALSE)

  expect_type(result, "list")
  expect_s3_class(result$plot, "ggplot")
})

test_that("y_fit_predict returns fitted values", {
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)
  result <- y_fit_predict(fit)

  expect_true(all(c("ypre", "ypre.Q2.5", "ypre.Q97.5") %in% names(result)))
  expect_equal(nrow(result), nrow(mtcars))
})

test_that("y_fit_summary supports plain lme4 models", {
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  result <- y_fit_summary(fit, output = "all", diagnostics = FALSE)

  expect_true("p" %in% names(result$Est.txt))
  expect_type(result, "list")
})

test_that("y_fit_classplot accepts a grouping data frame", {
  fit <- stats::lm(mpg ~ factor(cyl), data = mtcars)
  group_map <- data.frame(
    x = levels(factor(mtcars$cyl)),
    group = c("four", "six", "eight")
  )

  p <- y_fit_classplot(fit, group = group_map, hline = 0)
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})
