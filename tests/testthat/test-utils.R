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
