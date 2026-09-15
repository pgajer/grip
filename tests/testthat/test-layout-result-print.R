test_that("saved metric results retain objective and stopping information when printed", {
  skip_if_not_installed("smacof")
  expect_warning(fit <- metric.mds(edges = edges.cycle(6), n = 6,
                                   init = "random", max_iter = 1,
                                   diagnostics = FALSE), "iteration_limit")
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(fit, saved)
  restored <- readRDS(saved)
  local_mocked_bindings(score.gmds = function(...) stop("print must not recompute"))
  output <- capture.output(returned <- print(restored))
  expect_identical(returned, restored)
  expect_true(any(grepl("objective: raw distance stress", output, fixed = TRUE)))
  expect_true(any(grepl("squared input-distance units", output, fixed = TRUE)))
  expect_true(any(grepl("converged: no", output, fixed = TRUE)))
  expect_true(any(grepl("termination: iteration_limit", output, fixed = TRUE)))
  expect_true(any(grepl("selected start: 1 of 1", output, fixed = TRUE)))
  expect_true(any(grepl(format(restored$metadata$raw_stress, digits = 6, trim = TRUE),
                        output, fixed = TRUE)))
  expect_true(any(grepl("graph diagnostics: not computed", output, fixed = TRUE)))
})

test_that("direct scaling and edge-only results do not invent convergence or path scores", {
  edges <- edges.path(4)
  classical <- classical.mds(edges = edges, n = 4)
  output <- capture.output(print(classical))
  expect_true(any(grepl("classical strain", output, fixed = TRUE)))
  expect_true(any(grepl("not applicable (direct scaling)", output, fixed = TRUE)))
  expect_true(any(grepl("scale policy: profiled", output, fixed = TRUE)))
  fit <- edge.kk(coords = cbind(0:3, 0), prepared = prepare.edge.kk(edges), max_iter = 1)
  output <- capture.output(print(fit))
  expect_true(any(grepl("convergence: not reported", output, fixed = TRUE)))
  expect_true(any(grepl("retained-path target-normalized RMSE: unavailable", output, fixed = TRUE)))
  expect_true(any(grepl("chord target-normalized RMSE: unavailable", output, fixed = TRUE)))
  expect_true(any(grepl("edge-KK energy", output, fixed = TRUE)))
})

test_that("reported convergence and legacy results are printable without mutation", {
  fit <- classical.mds(edges = edges.path(4), n = 4, diagnostics = FALSE)
  fit$method <- "metric_mds"
  fit$metadata <- list(objective = "raw_distance_stress", raw_stress = 0,
                       converged = TRUE, termination = "stress_tolerance")
  before <- fit
  output <- capture.output(print(fit))
  expect_identical(fit, before)
  expect_true(any(grepl("converged: yes", output, fixed = TRUE)))
  fit$method <- "older_method"
  fit$metadata <- NULL
  expect_output(print(fit), "objective: not reported")
})
