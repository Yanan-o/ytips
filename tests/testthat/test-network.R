test_that("network helpers handle one or more layout groups", {
  skip_if_not_installed("igraph")

  g <- igraph::make_ring(6)
  igraph::E(g)$weight <- rep(1, igraph::ecount(g))
  igraph::V(g)$name <- paste0("N", seq_len(igraph::vcount(g)))

  expect_type(ytips:::yy_network_modularity(g), "list")
  expect_s3_class(ytips:::yy_network_properties(g), "data.frame")
  expect_type(ytips:::yy_network_zipi(g), "list")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (n_groups in 1:4) {
    cg <- data.frame(
      row.names = igraph::V(g)$name,
      c.group = rep(letters[seq_len(n_groups)], length.out = igraph::vcount(g))
    )
    vg <- data.frame(
      row.names = igraph::V(g)$name,
      v.group = rep(c("a", "b"), length.out = igraph::vcount(g))
    )
    expect_type(
      ytips:::yy_network_plot(
        g,
        c.group = cg,
        v.group = vg,
        v.size = 3,
        show.module = 1
      ),
      "list"
    )
  }
})

test_that("network plotting handles graphs without edges", {
  skip_if_not_installed("igraph")

  g <- igraph::make_empty_graph(3, directed = FALSE)
  igraph::V(g)$name <- paste0("N", 1:3)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_type(ytips:::yy_network_plot(g, show.module = 1), "list")
})

test_that("y_ps_filter honours prune", {
  skip_if_not_installed("phyloseq")

  otu <- data.frame(
    taxon = paste0("OTU", 1:4),
    S1 = c(1, 2, 3, 3),
    S2 = c(2, 3, 4, 3),
    S3 = c(3, 4, 5, 3),
    S4 = c(4, 5, 6, 3),
    check.names = FALSE
  )
  sample <- data.frame(
    sample = paste0("S", 1:4),
    trt = c("A", "A", "B", "B")
  )
  ps <- y_ps(otu, sample)

  expect_equal(phyloseq::ntaxa(y_ps_filter(ps, prune = FALSE)), 4)
  expect_equal(phyloseq::ntaxa(y_ps_filter(ps, prune = TRUE)), 3)
})
