test_that("plotting helpers return usable objects", {
  expect_length(color_y(5), 5)
  expect_s3_class(theme_y(), "theme")
  expect_s3_class(ytips:::no_legend, "theme")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()

  expect_s3_class(
    y_ggsave(
      p,
      save = FALSE,
      plot.background = ggplot2::element_rect(fill = "white", linewidth = 0)
    ),
    "ggplot"
  )
  expect_s3_class(
    y_classplot(mtcars[, c("mpg", "cyl", "am")], plot = "bar"),
    "ggplot"
  )
})

test_that("class plotting exposes only documented arguments", {
  expect_false("..." %in% names(formals(y_classplot)))
  expect_true("..." %in% names(formals(y_fit_classplot)))
  expect_true("..." %in% names(formals(y_fit_partial_plot)))
})

test_that("y_corrplot handles negative-only limits", {
  df <- data.frame(
    name_x = rep(c("a", "b"), each = 2),
    name_y = rep(c("a", "b"), 2),
    value = c(-1, -0.5, -0.5, -1),
    p.sig = "*"
  )

  p <- y_corrplot(df, fill.lim = c(-1, -0.2))
  built <- ggplot2::ggplot_build(p)
  expect_false(anyNA(built$data[[2]]$size))
})

test_that("partial plots support lm and lmer model options", {
  set.seed(1)
  df <- data.frame(
    block = rep(1:6, each = 4),
    x1 = stats::rnorm(24),
    x2 = stats::rnorm(24)
  )
  df$response <- df$x1 + 0.3 * df$x2 + stats::rnorm(24)

  lm_fit <- stats::lm(response ~ x1 + x2, data = df)
  lm_result <- y_fit_partial_plot(
    lm_fit,
    x.fixed = "x1"
  )
  expect_s3_class(lm_result$plot, "ggplot")

  lmer_fit <- lmerTest::lmer(
    response ~ x1 + x2 + (1 | block),
    data = df
  )
  lmer_result <- suppressWarnings(y_fit_partial_plot(
    lmer_fit,
    x.fixed = "x1"
  ))
  expect_s3_class(lmer_result$plot, "ggplot")
  expect_true(all(c("residual.y", "residual.x.fixed") %in%
    names(lmer_result$df_resid)))
})

test_that("decomposition plots accept custom response columns", {
  set.seed(2)
  df <- data.frame(
    block = rep(1:6, each = 4),
    pred1 = stats::rnorm(24),
    pred2 = stats::rnorm(24)
  )
  df$response <- df$pred1 + 0.5 * df$pred2 + stats::rnorm(24)
  df$between <- df$pred1 + stats::rnorm(24, sd = 0.1)
  df$within <- df$pred2 + stats::rnorm(24, sd = 0.1)
  lmer_fit <- lmerTest::lmer(
    response ~ pred1 + pred2 + (1 | block),
    data = df
  )

  result <- suppressWarnings(y_fit_decomposition_plot(
    lmer_fit,
    y1 = "between",
    y2 = "within",
    name.yy = c("component 1", "component 2")
  ))
  expect_s3_class(result$plot, "ggplot")
  expect_equal(nrow(result$df_p), 2)
  expect_equal(result$df_p$name, c("component 1", "component 2"))
  expect_equal(sum(result$df_p$explanatory.power), 1)

  lm_result <- y_fit_decomposition_plot(
    stats::lm(response ~ pred1 + pred2, data = df),
    y1 = "between",
    y2 = "within"
  )
  expect_s3_class(lm_result$plot, "ggplot")
  expect_equal(sum(lm_result$df_p$explanatory.power), 1)

  interaction_result <- y_fit_decomposition_plot(
    stats::lm(response ~ pred1 * pred2, data = df),
    y1 = "between",
    y2 = "within"
  )
  expect_s3_class(interaction_result$plot, "ggplot")
})
