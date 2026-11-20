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

test_that("fit helpers support generalized and mixed models", {
  glm_fit <- stats::glm(am ~ wt + hp, data = mtcars, family = stats::binomial)
  glm_result <- y_fit_summary(glm_fit, output = "all", diagnostics = FALSE)
  glm_pred <- y_fit_predict(glm_fit)
  expect_s3_class(glm_result$plot, "ggplot")
  expect_true(all(glm_pred$ypre >= 0 & glm_pred$ypre <= 1))

  lmer_fit <- lmerTest::lmer(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy
  )
  lmer_pred <- y_fit_predict(lmer_fit)
  expect_equal(nrow(lmer_pred), nrow(lme4::sleepstudy))
  expect_true(all(c("ypre", "ypre.Q2.5", "ypre.Q97.5") %in% names(lmer_pred)))
})

test_that("fit helpers support glmmTMB models", {
  skip_if_not_installed("glmmTMB")
  fit <- suppressWarnings(glmmTMB::glmmTMB(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy
  ))

  summary_result <- y_fit_summary(
    fit,
    output = "all",
    diagnostics = FALSE
  )
  prediction <- y_fit_predict(fit, newdata = lme4::sleepstudy[1:5, ])

  expect_s3_class(summary_result$plot, "ggplot")
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
  result <- y_fit_summary(fit, output = "all", diagnostics = FALSE)

  expect_true("p" %in% names(result$Est.txt))
  expect_type(result, "list")
})

test_that("y_fit_summary draws the random-forest overlay", {
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

  result <- y_fit_summary(
    fit,
    output = "all",
    diagnostics = FALSE,
    show.random.forest = TRUE
  )

  expect_s3_class(result$plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(result$plot))
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
    y_fit_summary(
      fit,
      output = "all",
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
