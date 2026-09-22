#!/usr/bin/env Rscript
# Standalone native probes under UndefinedBehaviorSanitizer (Clang/GCC).
# Run from the package root. Requires Rcpp and a compiler with UBSan support.
Sys.setenv(PKG_CXXFLAGS = "-fsanitize=undefined -fno-sanitize-recover=undefined",
           PKG_LIBS = "-fsanitize=undefined",
           UBSAN_OPTIONS = "halt_on_error=1:print_stacktrace=1")
code <- paste(c("// [[Rcpp::plugins(cpp17)]]",
                readLines("src/metric_mds_sgd.cpp")), collapse="\n")
Rcpp::sourceCpp(code=code, rebuild=TRUE)
set.seed(501)
for (dimension in 2:5) for (trial in 1:20) {
  n <- 5L + trial %% 10L
  x <- matrix(rnorm(n * dimension), n, dimension)
  if (trial %% 3 == 0) x[1, ] <- x[2, ]
  d <- runif(n * (n - 1) / 2)
  if (trial %% 5 == 0) d[1] <- 0
  fit <- grip_sgd_mds_cpp(x, d, c(.5, .2, .05), trial, 2L, 1e6)
  stopifnot(all(is.finite(fit$conf)), identical(fit$trace$epoch, c(0, 2, 3)))
  r <- as.double(dist(fit$conf))
  scale <- sum(r * d) / sum(r * r)
  stopifnot(abs(sum((scale * r - d)^2) - fit$stress) < 1e-10)
}
cat("80 native undefined-behavior-sanitizer probes passed\n")
