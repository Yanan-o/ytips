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
