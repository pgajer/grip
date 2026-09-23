# Internal controls for the native stress backend.
grip.mds.sgd.control <- function(control, max_iter) {
  defaults <- list(scheduler = "hybrid", learning_rate = 0.5,
                   final_rate = 0.01, switch_ratio = 0.4,
                   checkpoint_every = 1L, max_workspace_bytes = 256 * 1024^2)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) || any(names(control) == "") ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults))))) {
    stop("sgd_control must be a uniquely named list of documented controls", call. = FALSE)
  }
  defaults[names(control)] <- control
  z <- defaults
  if (!is.character(z$scheduler) || length(z$scheduler) != 1L ||
      is.na(z$scheduler) || !z$scheduler %in% c("hybrid", "exponential")) {
    stop("sgd_control$scheduler must be 'hybrid' or 'exponential'", call. = FALSE)
  }
  for (key in c("learning_rate", "final_rate", "max_workspace_bytes")) {
    grip.validate.scalar(z[[key]], paste0("sgd_control$", key), lower = 0, open.lower = TRUE)
  }
  if (z$final_rate > z$learning_rate) stop("final_rate must not exceed learning_rate", call. = FALSE)
  grip.validate.scalar(z$switch_ratio, "sgd_control$switch_ratio", lower = 0, upper = 1)
  if (!is.numeric(z$checkpoint_every) || length(z$checkpoint_every) != 1L ||
      !is.finite(z$checkpoint_every) || z$checkpoint_every < 1 ||
      z$checkpoint_every != floor(z$checkpoint_every) ||
      z$checkpoint_every > .Machine$integer.max) {
    stop("checkpoint_every must be a positive integer", call. = FALSE)
  }
  z$checkpoint_every <- as.integer(z$checkpoint_every)
  if (max_iter * 8 > z$max_workspace_bytes) {
    stop("SGD rate vector alone exceeds max_workspace_bytes", call. = FALSE)
  }
  z
}

grip.mds.sgd.rates <- function(control, max_iter) {
  t <- seq_len(max_iter) - 1
  a <- control$learning_rate
  final <- control$final_rate
  if (control$scheduler == "exponential") {
    rates <- exp(log(a) + (log(final) - log(a)) * (t / max_iter))
  } else {
    tau <- floor(control$switch_ratio * max_iter)
    middle <- if (tau > 0) a / 10 else a
    rates <- rep(middle, max_iter)
    early <- t < tau
    rates[early] <- exp(log(a) - log(10) * t[early] / tau)
    if (any(!early) && middle > final) {
      # Algebraically equivalent harmonic decay without overflowing middle/final.
      u <- (t[!early] - tau) / (max_iter - tau)
      rates[!early] <- final / ((final / middle) * (1 - u) + u)
    }
  }
  if (any(!is.finite(rates)) || any(rates <= 0)) {
    stop("SGD schedule has nonfinite or underflowed rates", call. = FALSE)
  }
  rates
}

grip.mds.sgd.fit <- function(start, target, rates, control, native_seed, weights = NULL) {
  grip_sgd_mds_cpp(start, target, rates, native_seed,
                   control$checkpoint_every, control$max_workspace_bytes, TRUE, weights)
}
