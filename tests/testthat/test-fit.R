test_that("y_fit_summary returns model information", {
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)

  expect_false("diagnostics" %in% names(formals(y_fit_summary)))
  expect_false("output" %in% names(formals(y_fit_summary)))
  result <- y_fit_summary(fit)

  expect_type(result, "list")
  expect_named(result, c("df", "name_y", "name_x", "Est.txt", "rr.txt"))
  expect_identical(result$name_y, "mpg")
  expect_identical(result$name_x, c("wt", "hp"))
  expect_equal(nrow(result$Est.txt), 2)
  expect_s3_class(y_fit_estplot(fit, diagnostics = FALSE), "ggplot")
})

test_that("y_fit_summary does not depend on model-data column positions", {
  df <- mtcars[, c("wt", "hp", "mpg")]
  result <- y_fit_summary(stats::lm(mpg ~ wt + hp, data = df))

  expect_identical(result$name_y, "mpg")
  expect_identical(result$name_x, c("wt", "hp"))
  expect_equal(nrow(result$Est.txt), 2)
})

test_that("y_fit_summary retains categorical and no-intercept coefficients", {
  factor_result <- y_fit_summary(stats::lm(mpg ~ factor(cyl), data = mtcars))
  no_intercept <- y_fit_summary(stats::lm(mpg ~ 0 + wt + hp, data = mtcars))

  expect_true(any(grepl("factor\\(cyl\\)", factor_result$Est.txt$x)))
  expect_equal(nrow(no_intercept$Est.txt), 2)
  expect_false(any(no_intercept$Est.txt$x == "(Intercept)"))
})

test_that("y_fit_predict returns fitted values", {
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)
  result <- y_fit_predict(fit)

  expect_true(all(c("ypre", "ypre.Q2.5", "ypre.Q97.5") %in% names(result)))
  expect_equal(nrow(result), nrow(mtcars))

  supplied <- y_fit_predict(fit, newdata = mtcars[1:5, ])
  expect_equal(nrow(supplied), 5)
})

test_that("y_fit_predict generates a prediction grid from the model", {
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)
  result <- y_fit_predict(fit, newdata = "new")

  expect_equal(nrow(result), 10000)
  expect_equal(range(result$wt), range(mtcars$wt))
  expect_equal(range(result$hp), range(mtcars$hp))
  expect_true(all(c("wt", "hp", "ypre", "ypre.Q2.5", "ypre.Q97.5") %in%
    names(result)))
})

test_that("prediction grids reduce resolution for three predictors", {
  fit <- stats::lm(mpg ~ wt + hp + drat, data = mtcars)
  result <- y_fit_predict(fit, newdata = "grid")

  expect_equal(nrow(result), 21^3)
  expect_equal(length(unique(result$wt)), 21)
  expect_equal(length(unique(result$hp)), 21)
  expect_equal(length(unique(result$drat)), 21)
})

test_that("automatic prediction grids preserve categorical levels", {
  fit <- stats::glm(
    am ~ wt + factor(cyl),
    data = mtcars,
    family = stats::binomial
  )
  result <- y_fit_predict(fit, newdata = "grid")

  expect_equal(sort(unique(result$cyl)), c(4, 6, 8))
  expect_true(all(result$ypre >= 0 & result$ypre <= 1))
  expect_false(anyNA(result$ypre))
})

test_that("fit helpers support generalized and mixed models", {
  glm_fit <- stats::glm(am ~ wt + hp, data = mtcars, family = stats::binomial)
  glm_result <- y_fit_summary(glm_fit)
  glm_pred <- y_fit_predict(glm_fit)
  expect_true(all(glm_pred$ypre >= 0 & glm_pred$ypre <= 1))
  expect_type(glm_result, "list")

  lmer_fit <- lmerTest::lmer(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy
  )
  lmer_pred <- y_fit_predict(lmer_fit)
  expect_equal(nrow(lmer_pred), nrow(lme4::sleepstudy))
  expect_true(all(c("ypre", "ypre.Q2.5", "ypre.Q97.5") %in% names(lmer_pred)))
})

test_that("y_fit_predict can retain random effects", {
  fit <- lmerTest::lmer(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy
  )
  newdata <- lme4::sleepstudy[c(1, 10, 20), ]

  population <- y_fit_predict(fit, newdata = newdata, re_formula = NA)
  conditional <- y_fit_predict(fit, newdata = newdata, re_formula = NULL)

  expect_false(isTRUE(all.equal(population$ypre, conditional$ypre)))
  expect_equal(
    conditional$ypre,
    unname(stats::predict(fit, newdata = newdata, re.form = NULL))
  )
})

test_that("fit helpers support glmmTMB models", {
  skip_if_not_installed("glmmTMB")
  fit <- suppressWarnings(glmmTMB::glmmTMB(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy
  ))

  summary_result <- y_fit_summary(fit)
  prediction <- y_fit_predict(fit, newdata = lme4::sleepstudy[1:5, ])

  expect_type(summary_result, "list")
  expect_equal(nrow(prediction), 5)
  expect_true(all(c("ypre", "ypre.Q2.5", "ypre.Q97.5") %in%
    names(prediction)))
})

test_that("y_fit_update orders fixed terms before random effects", {
  set.seed(1)
  df <- data.frame(
    block = factor(rep(1:8, each = 4)),
    z = stats::rnorm(32),
    x = stats::rnorm(32),
    q = stats::rnorm(32)
  )
  df$y <- df$z + df$x + stats::rnorm(32)

  fit <- lmerTest::lmer(y ~ x + z + (1 | block), data = df)
  updated <- y_fit_update(fit, str = "+ q")
  expect_equal(
    stats::formula(updated),
    y ~ z + x + q + (1 | block)
  )

  lm_fit <- stats::lm(y ~ x + z, data = df)
  lm_updated <- y_fit_update(lm_fit, str = "+ q", x.order = c("q", "x", "z"))
  expect_equal(stats::formula(lm_updated), y ~ q + x + z)

  interaction_fit <- lmerTest::lmer(y ~ x * z + (1 | block), data = df)
  interaction_updated <- y_fit_update(interaction_fit, str = "+ q")
  expect_equal(
    stats::formula(interaction_updated),
    y ~ z + x + q + x:z + (1 | block)
  )
})

test_that("y_fit_summary supports plain lme4 models", {
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  result <- y_fit_summary(fit)

  expect_true("p" %in% names(result$Est.txt))
  expect_type(result, "list")
})

test_that("y_fit_estplot draws the random-forest overlay", {
  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)
  local_mocked_bindings(
    yy_check_random_forest_dependencies = function() TRUE,
    y_random_forest_plot = function(...){
      list(RF.result = data.frame(
        name_x = c("wt", "hp"),
        X.IncMSE = c(20, 10)
      ))
    },
    .package = "ytips"
  )

  result <- y_fit_estplot(
    fit,
    diagnostics = FALSE,
    show.random.forest = TRUE
  )

  expect_s3_class(result, "ggplot")
  expect_silent(ggplot2::ggplot_build(result))
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

test_that("random-forest options report missing package installation commands", {
  random_forest_packages <- c("rfPermute", "randomForest")
  installed <- vapply(
    random_forest_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  skip_if(
    all(installed),
    "Random-forest dependencies are installed."
  )

  fit <- stats::lm(mpg ~ wt + hp, data = mtcars)

  expect_error(
    y_fit_estplot(
      fit,
      diagnostics = TRUE,
      show.random.forest = TRUE
    ),
    "Install with: install.packages"
  )
  expect_error(
    y_random_forest_plot(mtcars, x = c("wt", "hp"), y = "mpg"),
    "Install with: install.packages"
  )
})
