test_that("plotting helpers return usable objects", {
  expect_length(color_y(5), 5)
  expect_s3_class(theme_y(), "theme")
  expect_s3_class(ytips:::no_legend, "theme")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()

  expect_s3_class(y_ggsave(p, save = FALSE), "ggplot")
  expect_s3_class(
    y_classplot(mtcars[, c("mpg", "cyl", "am")], plot = "bar"),
    "ggplot"
  )
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

  lm_result <- y_partial_plot(
    df,
    y = "response",
    x = c("x1", "x2"),
    x.fix = "x1",
    fit = "lm"
  )
  expect_s3_class(lm_result$plot, "ggplot")

  lmer_result <- suppressWarnings(y_partial_plot(
    df,
    y = "response",
    x = c("x1", "x2"),
    x.fix = "x1",
    fit = "lmer"
  ))
  expect_s3_class(lmer_result$plot, "ggplot")
})

test_that("decomposition plots use the documented column names", {
  set.seed(2)
  df <- data.frame(
    block = rep(1:6, each = 4),
    x.fixed = stats::rnorm(24),
    y = stats::rnorm(24),
    y1 = stats::rnorm(24),
    y2 = stats::rnorm(24)
  )

  result <- suppressWarnings(y_decomposition_plot(
    df,
    y1y2.name = c("component 1", "component 2"),
    output = "all"
  ))
  expect_s3_class(result[[2]], "ggplot")
  expect_equal(nrow(result[[1]]), 2)
})
