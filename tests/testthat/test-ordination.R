test_that("ordination helpers run on vegan example data", {
  skip_if_not_installed("vegan")

  data(varespec, package = "vegan")
  data(varechem, package = "vegan")
  spe <- varespec[1:12, 1:20]
  env <- varechem[1:12, 1:3]
  trt <- factor(rep(c("CK", "A", "B"), each = 4))

  expect_s3_class(y_alpha_diversity(spe), "data.frame")
  expect_s3_class(y_pca_plot(spe, trt), "ggplot")
  expect_s3_class(y_pcoa_plot(spe, trt), "ggplot")
  expect_s3_class(y_rda_plot(spe, trt, env), "ggplot")
  expect_s3_class(y_cca_plot(spe, trt, env), "ggplot")
  expect_s3_class(y_db_rda_plot(spe, trt, env), "ggplot")
})

test_that("Price equation plotting uses default colors", {
  spe <- matrix(
    c(1, 2, 0, 1, 3, 1), nrow = 3,
    dimnames = list(c("CK1", "A1", "B1"), c("s1", "s2"))
  )
  price <- y_price_part(spe)

  expect_s3_class(y_price_part_plot(price), "ggplot")
})
