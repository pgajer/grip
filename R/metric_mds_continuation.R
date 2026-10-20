# Internal SGD result/state handling. Public entry point: metric.mds(resume=...).
.metric.mds.state.version <- function() "grip-sgd-continuation-v1"

.metric.mds.state.hash <- function(x) {
  grip_sgd_state_hash_cpp(serialize(x,NULL,version=2))
}

.metric.mds.sgd.status <- function(output, fit, control, kind) {
  trace <- fit$trace
  if (kind=="sparse") {
    trace$profiled_stress <- NULL
    names(trace)[names(trace)=="raw_stress"] <- "sparse_proxy_stress"
  }
  if (!is.null(fit$update_trace)) {
    output$update_trace <- fit$update_trace
    output$update_trace$coordinates <- lapply(fit$update_trace$coordinates,
      function(x) x * output$metadata$input_rms_distance)
    output$update_trace$coordinate_units <- "input_distance"
    output$update_trace$alignment <- "none"
    output$update_trace$scope <- "current_call"
  }
  output$trace <- trace
  # Consumers that retain metadata but omit top-level trace still keep history.
  output$metadata$history <- trace
  output$metadata$converged <- FALSE
  rejected <- identical(output$metadata$termination,"rejected_increase")
  output$metadata$termination <- if(rejected) "rejected_increase" else if(fit$stabilized) "practical_stabilization" else "iteration_limit"
  output$metadata$stopping <- list(
    reason=if(rejected) "rejected_increase" else if(fit$stabilized) "practical_stabilization" else "epoch_budget",
    schedule_completed=fit$niter>=control$schedule.epochs,
    schedule_epochs=control$schedule.epochs, epochs_completed=fit$niter,
    assessed=fit$assessments>0L, assessments=fit$assessments,
    stabilized=fit$stabilized, certified=FALSE,
    rule=control$stopping,
    returned_iterate=if(rejected) "scaled_initial" else if(kind=="full") "best_profiled_checkpoint" else "terminal",
    best_diagnostic=fit$best_diagnostic,best_diagnostic_epoch=fit$best_diagnostic_epoch,
    terminal_diagnostic=utils::tail(if(kind=="full") fit$trace$profiled_stress else fit$trace$raw_stress,1L),
    diagnostic=if(kind=="full") "RMS-normalized scale-profiled stress" else "asymmetric endpoint-average proxy",
    limitation="Practical stabilization is not certified convergence or global optimality; sparse updates need not minimize the proxy.")
  output
}

.metric.mds.keep.state <- function(output, fit, control, arguments, context) {
  if (!control$keep.state) return(output)
  output$continuation <- NULL
  state <- list(version=.metric.mds.state.version(),platform=R.version$platform,
    native=fit$state,arguments=arguments,control=control,context=context,
    history=fit$trace,template=output,stabilized=fit$stabilized)
  output$continuation <- list(state=state,checksum=.metric.mds.state.hash(state))
  output
}

.metric.mds.resume <- function(object, max.iter) {
  had.seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old.seed <- if(had.seed) get(".Random.seed",envir=.GlobalEnv) else NULL
  on.exit({
    if(had.seed) assign(".Random.seed",old.seed,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  max.iter <- grip.validate.vertex.count(max.iter)
  if (!inherits(object,"grip_gmds_layout") || !identical(object$method,"metric_mds") ||
      is.null(object$continuation))
    stop("resume requires an SGD metric.mds fit created with sgd.control=list(keep.state=TRUE); historical fits cannot be reconstructed",call.=FALSE)
  capsule <- object$continuation; state <- capsule$state
  if (!is.list(state) || !identical(state$version,.metric.mds.state.version()) ||
      !identical(state$platform,R.version$platform))
    stop("Unsupported continuation state version or platform",call.=FALSE)
  if (!identical(capsule$checksum,.metric.mds.state.hash(state)))
    stop("Continuation state integrity check failed; inputs, controls or saved state changed",call.=FALSE)
  if (identical(state$template$metadata$termination,"rejected_increase"))
    stop("Cannot resume a run whose fitted result was rejected for increased stress",call.=FALSE)
  if (isTRUE(state$stabilized))
    stop("This run already met its practical stopping rule; start a separate fit to change that rule",call.=FALSE)
  offset <- state$native$epoch
  if (as.double(offset)+max.iter>.Machine$integer.max)
    stop("Continuation epoch range exceeds the integer limit",call.=FALSE)
  total <- as.integer(offset+max.iter)
  control <- state$control; context <- state$context; args <- state$arguments
  if (as.double(max(total,control$schedule.epochs))*8>control$max.workspace.bytes)
    stop("Continuation rate vector exceeds max.workspace.bytes",call.=FALSE)
  weights <- if (context$kind=="full") context$weights else c(args$weights,args$reverse_weights)
  rates <- grip.mds.sgd.rates(control,total,weights)
  args$start <- state$native$coords;args$rates <- rates[seq.int(offset+1L,total)]
  args$state <- state$native;args$stopping <- control$stopping
  args$schedule_epochs <- control$schedule.epochs;args$keep_state <- TRUE
  args$update_trace <- control$trace
  began <- proc.time()[["elapsed"]]
  fit <- do.call(grip_sgd_mds_cpp,args)
  elapsed <- proc.time()[["elapsed"]]-began
  fit$trace <- rbind(state$history,fit$trace)
  rownames(fit$trace) <- NULL
  output <- state$template
  if (context$kind=="full") {
    x <- sweep(fit$conf,2L,colMeans(fit$conf),"-")
    d <- as.double(stats::dist(x)); w <- context$weights; target <- context$targets
    denominator <- sum(w*d^2)
    if (!is.finite(denominator) || denominator<=0) stop("Continued configuration is collapsed",call.=FALSE)
    multiplier <- sum(w*d*target)/denominator
    loss <- sum(w*(multiplier*d-target)^2)
    normalized <- multiplier*x;d <- as.double(stats::dist(normalized))
    target.scale <- sum(w*d*target)/context$energy
    output$coords <- normalized*context$rms
    rownames(output$coords) <- rownames(state$template$coords)
    output$metadata$coordinate_scale <- multiplier*context$rms
    output$metadata$raw_stress <- loss*context$loss.units
    output$metadata$target_normalized_rmse <- sqrt(loss/context$energy)
    output$metadata$stress1_identity <- sqrt(loss/sum(w*d^2))
    output$metadata$stress1_profiled <- sqrt(sum(w*(d-target.scale*target)^2)/sum(w*d^2))
    runs <- output$metadata$starts
    runs$raw_stress <- loss*context$loss.units;runs$normalized_stress <- loss/context$energy
    runs$iterations <- fit$niter;runs$converged <- FALSE
    runs$termination <- if(fit$stabilized) "practical_stabilization" else "iteration_limit"
    runs$backend_stress <- fit$stress;runs$elapsed_seconds <- runs$elapsed_seconds+elapsed
    runs$pair_updates <- fit$pair_updates;runs$best_epoch <- fit$best_epoch
    output$metadata$starts <- runs
    output$metadata$sgd[[1L]] <- list(trace=fit$trace,workspace_bytes=fit$workspace_bytes)
    if (context$diagnostics$enabled) {
      params <- context$diagnostics;params$enabled <- NULL
      output$diagnostics <- do.call(score.gmds,c(list(coords=output$coords,prepared=output$prepared),params))
    }
  } else {
    output$coords <- sweep(fit$terminal_conf,2L,colMeans(fit$terminal_conf),"-")*context$rms
    proxy <- utils::tail(fit$trace$raw_stress,1L)
    output$metadata$sparse_proxy_stress <- proxy
    output$metadata$sparse_proxy_normalized_rmse <- sqrt(proxy/context$energy)
    output$metadata$pair_updates <- fit$pair_updates
    output$metadata$fitting_seconds <- output$metadata$fitting_seconds+elapsed
    output$metadata$fitting_workspace_bytes <- fit$workspace_bytes
  }
  output$metadata$settings$max_iter <- total
  output <- .metric.mds.sgd.status(output,fit,control,context$kind)
  .metric.mds.keep.state(output,fit,control,state$arguments,context)
}
