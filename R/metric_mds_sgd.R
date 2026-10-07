# Internal controls for the native stress backend.
grip.mds.sgd.control <- function(control, max.iter) {
  defaults <- list(scheduler = "hybrid", learning.rate = 0.5,
                   final.rate = 0.01, switch.ratio = 0.4, schedule.epsilon = 0.1,
                   checkpoint.every = 1L, max.workspace.bytes = 256 * 1024^2)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) || any(names(control) == "") ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults))))) {
    stop("sgd.control must be a uniquely named list of documented controls", call. = FALSE)
  }
  defaults[names(control)] <- control
  z <- defaults
  if (!is.character(z$scheduler) || length(z$scheduler) != 1L ||
      is.na(z$scheduler) || !z$scheduler %in% c("hybrid", "exponential", "zheng")) {
    stop("sgd.control$scheduler must be 'hybrid', 'exponential' or 'zheng'", call. = FALSE)
  }
  if (z$scheduler == "zheng" && any(c("learning.rate","final.rate","switch.ratio") %in% names(control)))
    stop("zheng schedule derives rates from weights; omit learning.rate, final.rate and switch.ratio",call.=FALSE)
  grip.validate.scalar(z$schedule.epsilon,"sgd.control$schedule.epsilon",lower=0,upper=1,open.lower=TRUE)
  for (key in c("learning.rate", "final.rate", "max.workspace.bytes")) {
    grip.validate.scalar(z[[key]], paste0("sgd.control$", key), lower = 0, open.lower = TRUE)
  }
  if (z$final.rate > z$learning.rate) stop("final.rate must not exceed learning.rate", call. = FALSE)
  grip.validate.scalar(z$switch.ratio, "sgd.control$switch.ratio", lower = 0, upper = 1)
  if (!is.numeric(z$checkpoint.every) || length(z$checkpoint.every) != 1L ||
      !is.finite(z$checkpoint.every) || z$checkpoint.every < 1 ||
      z$checkpoint.every != floor(z$checkpoint.every) ||
      z$checkpoint.every > .Machine$integer.max) {
    stop("checkpoint.every must be a positive integer", call. = FALSE)
  }
  z$checkpoint.every <- as.integer(z$checkpoint.every)
  if (max.iter * 8 > z$max.workspace.bytes) {
    stop("SGD rate vector alone exceeds max.workspace.bytes", call. = FALSE)
  }
  z
}

grip.mds.sgd.rates <- function(control, max.iter, weights = NULL) {
  if (control$scheduler == "zheng") {
    positive <- weights[is.finite(weights) & weights > 0]
    if (!length(positive)) stop("zheng schedule requires positive weights",call.=FALSE)
    top <- -log(min(positive)); bottom <- log(control$schedule.epsilon)-log(max(positive))
    rates <- exp(if(max.iter==1L) top else seq(top,bottom,length.out=max.iter))
    if (any(!is.finite(rates)) || any(rates<=0)) stop("Published schedule rates are not representable",call.=FALSE)
    return(rates)
  }
  t <- seq_len(max.iter) - 1
  a <- control$learning.rate
  final <- control$final.rate
  if (control$scheduler == "exponential") {
    rates <- exp(log(a) + (log(final) - log(a)) * (t / max.iter))
  } else {
    tau <- floor(control$switch.ratio * max.iter)
    middle <- if (tau > 0) a / 10 else a
    rates <- rep(middle, max.iter)
    early <- t < tau
    rates[early] <- exp(log(a) - log(10) * t[early] / tau)
    if (any(!early) && middle > final) {
      # Algebraically equivalent harmonic decay without overflowing middle/final.
      u <- (t[!early] - tau) / (max.iter - tau)
      rates[!early] <- final / ((final / middle) * (1 - u) + u)
    }
  }
  if (any(!is.finite(rates)) || any(rates <= 0)) {
    stop("SGD schedule has nonfinite or underflowed rates", call. = FALSE)
  }
  rates
}

grip.mds.sgd.fit <- function(start, target, rates, control, native.seed, weights = NULL) {
  grip_sgd_mds_cpp(start, target, rates, native.seed,
                   control$checkpoint.every, control$max.workspace.bytes, TRUE, weights)
}
