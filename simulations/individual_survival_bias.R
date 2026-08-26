# Utilities for the individual survival-curve bias analysis.
# A bias is defined as predicted survival minus true survival.

empty_method_result = function(metric_names, n_subject, pred_times) {
  list(
    metrics = setNames(rep(NA_real_, length(metric_names)), metric_names),
    individual_survival_bias = matrix(
      NA_real_, nrow = n_subject, ncol = length(pred_times)
    )
  )
}

valid_simulation_result = function(res, expected_bias_dim) {
  is.list(res) && is.matrix(res$metrics) &&
    identical(dim(res$individual_survival_bias), as.integer(expected_bias_dim))
}

combine_method_results = function(method_results, method_names,
                                  subject_ids, pred_times) {
  stopifnot(length(method_results) == length(method_names))
  n_subject = length(subject_ids)
  n_time = length(pred_times)
  n_method = length(method_names)
  
  metrics = do.call(rbind, lapply(method_results, function(x) x$metrics))
  bias = array(
    NA_real_,
    dim = c(n_subject, n_time, n_method)
  )
  for(m in seq_len(n_method)) {
    cur_bias = method_results[[m]]$individual_survival_bias
    if(!identical(dim(cur_bias), c(n_subject, n_time))) {
      stop("Unexpected dimensions for an individual survival-bias matrix.")
    }
    bias[,,m] = cur_bias
  }
  
  list(
    metrics = cbind(seq_len(n_method), metrics),
    individual_survival_bias = bias
  )
}

summarize_individual_survival_bias = function(sim_results, validation_data,
                                               pred_times, method_names,
                                               scenario = NULL) {
  if(length(sim_results)==0L) {
    stop("sim_results must contain at least one simulation replication.")
  }
  n_subject = nrow(validation_data)
  expected_dim = c(n_subject, length(pred_times), length(method_names))
  valid_rep = vapply(
    sim_results, valid_simulation_result, logical(1),
    expected_bias_dim = expected_dim
  )
  if(!any(valid_rep)) {
    stop("No simulation replication returned a valid survival-bias result.")
  }
  if(any(!valid_rep)) {
    warning(sum(!valid_rep), " simulation replication(s) failed and were excluded.")
  }
  bias_sum = array(0.0, dim = expected_dim)
  bias_n = array(0L, dim = expected_dim)
  
  for(res in sim_results[valid_rep]) {
    cur_bias = res$individual_survival_bias
    keep = is.finite(cur_bias)
    bias_sum[keep] = bias_sum[keep]+cur_bias[keep]
    bias_n[keep] = bias_n[keep]+1L
  }
  
  mean_bias = bias_sum/bias_n
  mean_bias[bias_n==0L] = NA_real_
  dimnames(mean_bias) = list(
    subject_id = as.character(validation_data[,"ID"]),
    prediction_time = as.character(pred_times),
    method = method_names
  )
  dimnames(bias_n) = dimnames(mean_bias)

  failed_rep = which(!valid_rep)
  failure_messages = vapply(
    sim_results[failed_rep],
    function(x) {
      if(inherits(x, "try-error")) {
        paste(as.character(x), collapse = " ")
      } else {
        paste("Invalid simulation result of class", paste(class(x), collapse = "/"))
      }
    },
    character(1)
  )
  
  bias_data = expand.grid(
    subject_id = validation_data[,"ID"],
    prediction_time = pred_times,
    method = method_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  bias_data$mean_bias = as.vector(mean_bias)
  bias_data$n_rep = as.vector(bias_n)
  x_names = setdiff(colnames(validation_data), c("ID", "Y", "Delta"))
  subject_row = match(bias_data$subject_id, validation_data[,"ID"])
  for(x_name in x_names) {
    bias_data[[x_name]] = validation_data[subject_row, x_name]
  }
  if(!is.null(scenario)) {
    bias_data$scenario = scenario
  }
  
  list(
    mean_bias = mean_bias,
    n_rep = bias_n,
    data = bias_data,
    prediction_times = pred_times,
    methods = method_names,
    validation_data = validation_data,
    requested_replications = length(sim_results),
    successful_replications = sum(valid_rep),
    failed_replications = failed_rep,
    failure_messages = failure_messages
  )
}

plot_individual_survival_bias = function(bias_summary, file) {
  mean_bias = bias_summary$mean_bias
  pred_times = bias_summary$prediction_times
  method_names = bias_summary$methods
  finite_bias = mean_bias[is.finite(mean_bias)]
  if(length(finite_bias)==0L) {
    warning("No finite individual survival-bias estimates to plot.")
    return(invisible(NULL))
  }
  
  y_lim = range(finite_bias)
  if(diff(y_lim)<1e-8) {
    y_lim = y_lim+c(-0.01, 0.01)
  } else {
    y_lim = y_lim+c(-0.05, 0.05)*diff(y_lim)
  }
  n_time = length(pred_times)
  n_col = min(2L, n_time)
  n_row = ceiling(n_time/n_col)
  method_col = grDevices::hcl.colors(length(method_names), "Dark 3")
  
  grDevices::pdf(file, width = 4.5*n_col, height = 4.5*n_row)
  old_par = graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(5, 4, 3, 1))
  
  for(j in seq_len(n_time)) {
    values = lapply(seq_along(method_names), function(m) mean_bias[,j,m])
    keep = vapply(values, function(x) any(is.finite(x)), logical(1))
    if(!any(keep)) {
      graphics::plot.new()
      graphics::title(main = sprintf("Prediction time = %g", pred_times[j]))
      next
    }
    graphics::boxplot(
      values[keep], at = which(keep), names = FALSE,
      xlim = c(0.5, length(method_names)+0.5), ylim = y_lim,
      col = method_col[keep], outline = FALSE,
      ylab = "Mean bias across replications",
      main = sprintf("Prediction time = %g", pred_times[j])
    )
    graphics::axis(
      1, at = seq_along(method_names), labels = method_names,
      las = 1, cex.axis = 0.8
    )
    graphics::abline(h = 0, lty = 2, col = "grey35")
    if(any(!keep)) {
      graphics::text(
        which(!keep), y_lim[1]+0.04*diff(y_lim), "NA",
        col = "grey45", cex = 0.75
      )
    }
  }
  if(n_time < n_row*n_col) {
    graphics::plot.new()
  }
  invisible(file)
}
