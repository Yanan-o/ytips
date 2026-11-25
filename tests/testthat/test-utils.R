test_that("onedrive returns paths or first-level folder contents", {
  root <- file.path(withr::local_tempdir(), "OneDrive-Test")
  dir.create(root)
  test_name <- paste0("ytips-test-", basename(root))
  file.create(file.path(root, paste0(test_name, ".R")))
  dir.create(file.path(root, test_name))
  file.create(file.path(root, test_name, "inside.csv"))
  withr::local_envvar(
    OneDrive = root,
    OneDriveConsumer = "",
    OneDriveCommercial = ""
  )

  expect_true(any(basename(onedrive(paste0(test_name, ".R"))) == paste0(test_name, ".R")))
  expect_true(any(basename(onedrive(test_name)) == test_name))
  expect_true(any(basename(onedrive(test_name, contents = TRUE)) == "inside.csv"))
})

test_that("onedrive supports parent/file paths and multiple OneDrive roots", {
  base <- withr::local_tempdir()
  root1 <- file.path(base, "OneDrive-Test1")
  root2 <- file.path(base, "OneDrive-Test2")
  dir.create(root1)
  dir.create(root2)
  dir.create(file.path(root1, "parent"))
  dir.create(file.path(root2, "parent"))
  file.create(file.path(root1, "parent", "same.R"))
  file.create(file.path(root2, "parent", "same.R"))
  withr::local_envvar(
    OneDrive = root1,
    OneDriveConsumer = root2,
    OneDriveCommercial = ""
  )

  expect_length(onedrive("same.R"), 2)
  expect_length(onedrive("parent/same.R"), 2)
  expect_silent(onedrive("parent/same.R"))
})

test_that("onedrive rejects missing names and non-OneDrive roots", {
  root <- withr::local_tempdir()
  test_name <- paste0("ytips-test-", basename(root))
  file.create(file.path(root, test_name))
  withr::local_envvar(
    OneDrive = root,
    OneDriveConsumer = "",
    OneDriveCommercial = ""
  )

  expect_error(onedrive(), "x must be one file or folder name")
  expect_false(test_name %in% basename(onedrive(test_name)))
})

test_that("dependency checks provide actionable installation commands", {
  local_mocked_bindings(
    yy_pkg_available = function(pkg) FALSE,
    .package = "ytips"
  )

  expect_error(ytips:::yy_check_phyloseq(), "bioc::phyloseq")
  expect_error(
    ytips:::yy_check_random_forest_dependencies(),
    "install.packages"
  )
  expect_error(ytips:::yy_check_sem_dependencies(), "piecewiseSEM")
})

test_that("dependency checks offer an interactive installation choice", {
  marker <- tempfile()
  local_mocked_bindings(
    yy_pkg_available = function(pkg) file.exists(marker),
    .package = "ytips"
  )

  with_mocked_bindings(
    menu = function(...) 1L,
    .package = "utils",
    expect_true(ytips:::yy_check_pkg(
      "example.package",
      install.hint = paste0("file.create(", shQuote(marker), ")"),
      ask = TRUE
    ))
  )
  expect_true(file.exists(marker))

  unlink(marker)
  with_mocked_bindings(
    menu = function(...) 2L,
    .package = "utils",
    expect_error(
      ytips:::yy_check_pkg(
        "example.package",
        install.hint = "install.packages('example.package')",
        ask = TRUE
      ),
      "Install with"
    )
  )
})

test_that("dependency checks use an RStudio dialog when available", {
  marker <- tempfile()
  local_mocked_bindings(
    yy_pkg_available = function(pkg) file.exists(marker),
    .package = "ytips"
  )

  with_mocked_bindings(
    isAvailable = function(...) TRUE,
    showQuestion = function(...) TRUE,
    .package = "rstudioapi",
    expect_true(ytips:::yy_check_pkg(
      "example.package",
      install.hint = paste0("file.create(", shQuote(marker), ")"),
      ask = TRUE
    ))
  )
  expect_true(file.exists(marker))
})

test_that("network dependency checks follow the selected method", {
  local_mocked_bindings(
    yy_pkg_available = function(pkg) pkg %in% c("phyloseq", "igraph"),
    .package = "ytips"
  )

  expect_error(
    ytips:::yy_check_network_dependencies(),
    "Hmisc"
  )
  expect_error(
    ytips:::yy_check_network_dependencies(spiec.easi = TRUE),
    "SpiecEasi"
  )
  expect_true(
    ytips:::yy_check_network_dependencies(taxa.shuffle = TRUE)
  )
})

test_that("public entry points fail early with installation instructions", {
  local_mocked_bindings(
    yy_pkg_available = function(pkg) FALSE,
    .package = "ytips"
  )

  expect_error(y_ps(NULL, NULL), "bioc::phyloseq")
  expect_error(y_ps_filter(NULL), "bioc::phyloseq")
  expect_error(y_ps_bar(NULL, pairwise = TRUE), "rcompanion")
  expect_error(y_ps_igraph(NULL), "bioc::phyloseq")
  expect_error(y_ps_network(NULL), "bioc::phyloseq")
  expect_error(y_psem_boot(NULL, NULL), "piecewiseSEM")
  expect_error(y_psem_step(NULL), "piecewiseSEM")
  expect_error(y_ps_to_spe(NULL), "bioc::phyloseq")
  expect_error(ytips:::yy_adonis_summary(NULL), "pairwiseAdonis")
  expect_error(y_beta_part(NULL), "betapart")
  expect_error(y_rdacca_hp(NULL, NULL), "rdacca.hp")
  expect_error(
    yy_ord_plot(
      data.frame(Dim1 = 1:2, Dim2 = 1:2, trt = c("a", "b")),
      show.density = TRUE
    ),
    "ggExtra"
  )
  expect_error(y_descdist(1:10), "fitdistrplus")
  expect_error(onedrive("example.csv"), "fs")
  expect_error(
    y_classplot(mtcars[, c("mpg", "cyl", "am")], plot = "vio"),
    "ggbeeswarm"
  )
  expect_error(
    y_classplot(mtcars[, c("mpg", "cyl", "am")], plot = "viobox"),
    "ggdist"
  )
  expect_error(
    y_fit_classplot(stats::lm(mpg ~ factor(cyl), data = mtcars)),
    "car"
  )
  expect_error(
    y_fit_estplot(
      stats::lm(mpg ~ wt + hp, data = mtcars),
      diagnostics = TRUE
    ),
    "DHARMa"
  )
})

test_that("fitted-class plots check their conditional dependencies", {
  fit <- stats::lm(mpg ~ factor(cyl), data = mtcars)

  local_mocked_bindings(
    yy_pkg_available = function(pkg) pkg == "car",
    .package = "ytips"
  )
  expect_error(
    y_fit_classplot(fit, pairwise = TRUE),
    "rcompanion"
  )
})

test_that("fitted-class brackets require ggpubr", {
  fit <- stats::lm(count ~ factor(spray), data = InsectSprays)

  local_mocked_bindings(
    yy_pkg_available = function(pkg) pkg != "ggpubr",
    .package = "ytips"
  )
  expect_error(
    y_fit_classplot(fit, pairwise = FALSE),
    "ggpubr"
  )
})
