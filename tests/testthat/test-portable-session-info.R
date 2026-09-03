test_that("session metadata retains versions without local file attributes", {
  helper <- new.env(parent = baseenv())
  sys.source(system.file("scripts", "portable-session-info.R", package = "grip"),
             envir = helper)
  description <- structure(list(Package = "example", Version = "1.2.3"),
                           file = "/home/example/library/example/DESCRIPTION")
  info <- structure(list(otherPkgs = list(example = description),
                         loadedOnly = list(example = description),
                         platform = "test-platform"), class = "sessionInfo")
  clean <- helper$portable_session_info(info)
  expect_null(attr(clean$otherPkgs$example, "file"))
  expect_null(attr(clean$loadedOnly$example, "file"))
  expect_identical(clean$otherPkgs$example$Version, "1.2.3")
  expect_identical(clean$platform, info$platform)
  expect_s3_class(clean, "sessionInfo")
  expect_identical(helper$portable_session_info(clean), clean)
  expect_identical(attr(info$otherPkgs$example, "file"),
                   "/home/example/library/example/DESCRIPTION")
})

test_that("bundled benchmark sessions do not retain DESCRIPTION paths", {
  for (artifact in c("hmp_u01_gc_coarse/vignette_results.rds",
                     "vs_alternatives/benchmark_results.rds")) {
    results <- readRDS(system.file("extdata", artifact, package = "grip"))
    expect_s3_class(results$session_info, "sessionInfo")
    for (section in c("otherPkgs", "loadedOnly")) {
      for (description in results$session_info[[section]]) {
        expect_null(attr(description, "file"))
      }
    }
  }
})
