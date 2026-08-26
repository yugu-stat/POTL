# Create the two main-paper figures for the multi-source simulation.
#
# Each run of sim_multi_source.R produces one joint multi-source result for a
# pooling rule ("cv" or "ns"), rather than five scenario-specific results.
# This script combines the two runs: POTL-N is read from the ns files, while
# POTL-CV and all comparator methods are read from the cv files. It reads the
# four saved RData files from its own directory by default and writes:
#   pred_error_multi_source.pdf
#   indiv_bias_multi_source.pdf
#
# Usage:
#   Rscript plot_multi_source_figures.R
#   Rscript plot_multi_source_figures.R <data_dir>

suppressPackageStartupMessages(library(ggplot2))

full_args = commandArgs(trailingOnly = FALSE)
script_arg = grep("^--file=", full_args, value = TRUE)
if(length(script_arg) != 1L) {
  stop("Run this plotting program with Rscript.")
}
script_file = normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir = dirname(script_file)

args = commandArgs(trailingOnly = TRUE)
if(length(args) > 1L) {
  stop("Usage: Rscript plot_multi_source_figures.R [data_dir]")
}

data_dir = if(length(args) == 1L) args[1] else script_dir
data_dir = normalizePath(data_dir, mustWork = TRUE)

design = "multi_source"
scalar_paths = c(
  ns = file.path(data_dir, "res_multi_source_ns.RData"),
  cv = file.path(data_dir, "res_multi_source_cv.RData")
)
bias_paths = c(
  ns = file.path(data_dir, "individual_survival_bias_multi_source_ns.RData"),
  cv = file.path(data_dir, "individual_survival_bias_multi_source_cv.RData")
)

source_method_levels = c("POTL", "Target-only", "CoxTL", "Pooled")
method_levels = c("POTL-N", "POTL-CV", "Target-only", "CoxTL", "Pooled")
method_colors = c(
  "POTL-N" = "#0072B2",
  "POTL-CV" = "#D55E00",
  "Target-only" = "black",
  "CoxTL" = "green",
  "Pooled" = "purple"
)
method_shapes = c(
  "POTL-N" = 16,
  "POTL-CV" = 15,
  "Target-only" = 17,
  "CoxTL" = 18,
  "Pooled" = 8
)

load_saved_object = function(path, object_name) {
  if(!file.exists(path)) {
    stop("Saved simulation result not found: ", path)
  }
  env = new.env(parent = emptyenv())
  loaded = load(path, envir = env)
  if(!object_name %in% loaded) {
    stop("Expected object '", object_name, "' was not found in ", path)
  }
  env[[object_name]]
}

# ---------------------------------------------------------------------------
# Prediction errors: median across replications, with median +/- scaled MAD.
# ---------------------------------------------------------------------------

read_scalar_results = function(path, pool) {
  allres = load_saved_object(path, "allres")
  if(!is.list(allres) || length(allres) == 0L) {
    stop("allres must be a nonempty list in ", path)
  }

  method_indices = if(pool == "ns") 1L else seq_along(source_method_levels)
  method_names = if(pool == "ns") {
    "POTL-N"
  } else {
    c("POTL-CV", source_method_levels[-1])
  }
  data = do.call(rbind, lapply(seq_along(allres), function(replication) {
    value = allres[[replication]]
    if(!is.matrix(value) || !is.numeric(value) ||
       !identical(dim(value), c(length(source_method_levels), 6L))) {
      stop(
        "Expected a 4 x 6 numeric result matrix in ", path,
        ", replication ", replication, "."
      )
    }
    method_id = as.integer(value[,1])
    if(!identical(method_id, seq_along(source_method_levels)) ||
       any(value[,1] != method_id)) {
      stop(
        "Expected method IDs 1:4 in ", path,
        ", replication ", replication, "."
      )
    }
    errors = value[method_indices,c(2, 3), drop = FALSE]
    if(any(!is.finite(errors)) || any(errors < 0)) {
      stop(
        "Selected L2D and D_tau values must be finite and nonnegative in ",
        path, ", replication ", replication, "."
      )
    }
    data.frame(
      replication = replication,
      method = method_names,
      L2D = errors[,1],
      Dtau = errors[,2],
      stringsAsFactors = FALSE
    )
  }))
  list(data = data, n_replications = length(allres))
}

scalar_results = list(
  ns = read_scalar_results(scalar_paths[["ns"]], "ns"),
  cv = read_scalar_results(scalar_paths[["cv"]], "cv")
)
scalar_data = rbind(scalar_results$ns$data, scalar_results$cv$data)
scalar_data$method = factor(scalar_data$method, levels = method_levels)
n_scalar_replications = c(
  ns = scalar_results$ns$n_replications,
  cv = scalar_results$cv$n_replications
)
if(n_scalar_replications[["ns"]] != n_scalar_replications[["cv"]]) {
  stop("The ns and cv scalar files must contain the same number of replications.")
}

summarize_scalar_metric = function(variable, label) {
  rows = lapply(method_levels, function(method) {
    values = scalar_data[[variable]][scalar_data$method == method]
    expected = if(method == "POTL-N") {
      n_scalar_replications[["ns"]]
    } else {
      n_scalar_replications[["cv"]]
    }
    if(length(values) != expected) {
      stop("Incomplete scalar results for method ", method, ".")
    }
    center = stats::median(values)
    spread = stats::mad(values, center = center, constant = 1.4826)
    data.frame(
      method = method,
      metric = label,
      median = center,
      mad = spread,
      ymin = max(0, center-spread),
      ymax = center+spread,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

scalar_summary = rbind(
  summarize_scalar_metric("L2D", "L[2]*D"),
  summarize_scalar_metric("Dtau", "D[tau]")
)
scalar_summary$method = factor(scalar_summary$method, levels = method_levels)
scalar_summary$metric = factor(
  scalar_summary$metric, levels = c("L[2]*D", "D[tau]")
)

scalar_ymax = ceiling(100*max(scalar_summary$ymax))/100
prediction_error_plot = ggplot(
  scalar_summary,
  aes(x = method, y = median, color = method, shape = method)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax), width = 0.12, linewidth = 0.45
  ) +
  geom_point(size = 2.4, stroke = 0.45) +
  facet_wrap(~metric, nrow = 1, labeller = label_parsed) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, scalar_ymax),
    breaks = pretty(c(0, scalar_ymax), n = 5),
    expand = expansion(mult = c(0, 0.025))
  ) +
  labs(x = "Method", y = "Prediction error") +
  theme_bw(base_size = 9, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linewidth = 0.45),
    strip.background = element_rect(fill = "grey93", linewidth = 0.45),
    strip.text = element_text(face = "bold", size = 9),
    axis.title.x = element_text(margin = margin(t=6)),
    axis.title.y = element_text(margin = margin(r=6)),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(5, 6, 4, 5, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Individual survival bias. For subject i, method m, and prediction time t,
# mean_bias[i,t,m] = mean_r{S_hat_rim(t) - S_0(t | X_i)}. The figure shows
# the distribution of this signed mean bias across validation subjects.
# ---------------------------------------------------------------------------

read_individual_bias = function(path, pool, scalar_replications) {
  individual_bias = load_saved_object(path, "individual_bias")
  required = c(
    "mean_bias", "n_rep", "methods", "prediction_times", "validation_data",
    "requested_replications", "successful_replications"
  )
  if(!is.list(individual_bias) || !all(required %in% names(individual_bias))) {
    stop("Unexpected individual_bias structure in ", path)
  }
  if(!identical(as.character(individual_bias$methods), source_method_levels)) {
    stop("Unexpected method order in ", path)
  }

  prediction_times = as.numeric(individual_bias$prediction_times)
  if(length(prediction_times) != 4L || any(!is.finite(prediction_times)) ||
     any(diff(prediction_times) <= 0)) {
    stop("Expected four strictly increasing prediction times in ", path)
  }
  if(is.null(nrow(individual_bias$validation_data)) ||
     nrow(individual_bias$validation_data) < 1L) {
    stop("validation_data must contain at least one subject in ", path)
  }

  source_dim = as.integer(c(
    nrow(individual_bias$validation_data),
    length(prediction_times),
    length(source_method_levels)
  ))
  if(!identical(dim(individual_bias$mean_bias), source_dim) ||
     !identical(dim(individual_bias$n_rep), source_dim)) {
    stop("Unexpected individual-bias dimensions in ", path)
  }

  method_indices = if(pool == "ns") 1L else seq_along(source_method_levels)
  selected_bias = individual_bias$mean_bias[,,method_indices, drop = FALSE]
  selected_n_rep = individual_bias$n_rep[,,method_indices, drop = FALSE]
  if(!is.numeric(individual_bias$mean_bias) ||
     any(!is.finite(selected_bias))) {
    stop("Selected mean_bias values must be finite and numeric in ", path)
  }
  if(any(abs(selected_bias) > 1+sqrt(.Machine$double.eps))) {
    stop("Selected mean_bias values lie outside the survival-probability scale in ", path)
  }

  requested = as.integer(individual_bias$requested_replications)
  successful = as.integer(individual_bias$successful_replications)
  if(length(requested) != 1L || length(successful) != 1L ||
     is.na(requested) || is.na(successful) || requested < 1L ||
     successful != requested) {
    stop("Not all requested individual-bias replications succeeded in ", path)
  }
  if(anyNA(selected_n_rep) || any(selected_n_rep != successful)) {
    stop("Some selected individual-bias cells do not use all replications in ", path)
  }
  if(scalar_replications != successful) {
    stop("Scalar and individual-bias replication counts do not agree in ", path)
  }
  individual_bias
}

individual_bias_ns = read_individual_bias(
  bias_paths[["ns"]], "ns", n_scalar_replications[["ns"]]
)
individual_bias_cv = read_individual_bias(
  bias_paths[["cv"]], "cv", n_scalar_replications[["cv"]]
)
if(!identical(
  as.numeric(individual_bias_ns$prediction_times),
  as.numeric(individual_bias_cv$prediction_times)
)) {
  stop("The ns and cv prediction times do not agree.")
}
if(!isTRUE(all.equal(
  individual_bias_ns$validation_data,
  individual_bias_cv$validation_data,
  check.attributes = TRUE
))) {
  stop("The ns and cv validation data do not agree.")
}

prediction_times = as.numeric(individual_bias_cv$prediction_times)
expected_dim = as.integer(c(
  nrow(individual_bias_cv$validation_data),
  length(prediction_times),
  length(method_levels)
))
mean_bias = array(NA_real_, dim = expected_dim)
mean_bias[,,1] = individual_bias_ns$mean_bias[,,1]
mean_bias[,,2:5] = individual_bias_cv$mean_bias

safe_bandwidth = function(values) {
  if(length(unique(values)) < 2L) {
    return(NA_real_)
  }
  bandwidth = suppressWarnings(stats::bw.nrd0(values))
  if(is.finite(bandwidth) && bandwidth > 0) bandwidth else NA_real_
}

bias_pp = 100*mean_bias
index = expand.grid(
  subject = seq_len(expected_dim[1]),
  time_index = seq_len(expected_dim[2]),
  method_index = seq_len(expected_dim[3]),
  KEEP.OUT.ATTRS = FALSE
)
bias_data = data.frame(
  prediction_time = prediction_times[index$time_index],
  method = method_levels[index$method_index],
  bias_pp = as.vector(bias_pp),
  stringsAsFactors = FALSE
)

bias_summary_rows = vector(
  "list", length(prediction_times)*length(method_levels)
)
k = 0L
for(time_index in seq_along(prediction_times)) {
  for(method_index in seq_along(method_levels)) {
    k = k+1L
    values = bias_pp[,time_index,method_index]
    q = stats::quantile(
      values, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), names = FALSE
    )
    bias_summary_rows[[k]] = data.frame(
      prediction_time = prediction_times[time_index],
      method = method_levels[method_index],
      q05 = q[1], q25 = q[2], median = q[3], q75 = q[4], q95 = q[5],
      mean = mean(values),
      bandwidth = safe_bandwidth(values),
      stringsAsFactors = FALSE
    )
  }
}
bias_summary = do.call(rbind, bias_summary_rows)

time_labels = paste0(
  "t = ", format(prediction_times, trim = TRUE, scientific = FALSE)
)
bias_data$prediction_time = factor(
  bias_data$prediction_time, levels = prediction_times, labels = time_labels
)
bias_data$method = factor(bias_data$method, levels = method_levels)
bias_summary$prediction_time = factor(
  bias_summary$prediction_time, levels = prediction_times, labels = time_labels
)
bias_summary$method = factor(bias_summary$method, levels = method_levels)

nice_symmetric_limit = function(max_abs) {
  target = 1.05*max_abs
  if(!is.finite(target) || target <= 0) {
    return(0.1)
  }
  magnitude = 10^floor(log10(target))
  candidates = c(1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 5.5, 6, 8, 10)*magnitude
  candidates[which(candidates >= target)[1]]
}

# Absolute values are used only to make the signed-bias scale symmetric around
# zero; they are not averaged or reported as an error metric.
bias_limit = nice_symmetric_limit(max(abs(bias_data$bias_pp)))
bandwidth_candidates = bias_summary$bandwidth[
  is.finite(bias_summary$bandwidth) & bias_summary$bandwidth > 0
]
bias_bandwidth = if(length(bandwidth_candidates) > 0L) {
  stats::median(bandwidth_candidates)
} else {
  bias_limit/50
}

limit_data = expand.grid(
  prediction_time = time_labels,
  endpoint = c(-1, 1),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
limit_data$prediction_time = factor(
  limit_data$prediction_time, levels = time_labels
)
limit_data$method = factor("POTL-N", levels = method_levels)
limit_data$bias_pp = bias_limit*limit_data$endpoint

individual_bias_plot = ggplot(
  bias_data,
  aes(x = method, y = bias_pp, fill = method)
) +
  geom_blank(
    data = limit_data,
    aes(x = method, y = bias_pp),
    inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = 0, linetype = 2, linewidth = 0.4, color = "grey25"
  ) +
  geom_violin(
    bw = bias_bandwidth, scale = "width", trim = TRUE,
    linewidth = 0.25, alpha = 0.78, color = "grey20"
  ) +
  geom_linerange(
    data = bias_summary,
    inherit.aes = FALSE,
    aes(x = method, ymin = q05, ymax = q95),
    linewidth = 0.35, color = "black"
  ) +
  geom_linerange(
    data = bias_summary,
    inherit.aes = FALSE,
    aes(x = method, ymin = q25, ymax = q75),
    linewidth = 1.15, color = "black"
  ) +
  geom_point(
    data = bias_summary,
    inherit.aes = FALSE,
    aes(x = method, y = median),
    shape = 21, size = 1.45, stroke = 0.35,
    fill = "white", color = "black"
  ) +
  geom_point(
    data = bias_summary,
    inherit.aes = FALSE,
    aes(x = method, y = mean),
    shape = 23, size = 1.25, stroke = 0.30,
    fill = "black", color = "black"
  ) +
  facet_wrap(~prediction_time, nrow = 1) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Individual mean bias in predicted survival (\u00d7 100)",
    fill = NULL
  ) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
    panel.border = element_rect(linewidth = 0.4),
    strip.background = element_rect(fill = "grey93", linewidth = 0.4),
    strip.text = element_text(size = 8.5),
    axis.title.y = element_text(size = 9),
    axis.text.y = element_text(size = 7.2, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(12, "pt"),
    legend.spacing.x = grid::unit(3, "pt"),
    legend.margin = margin(t = -2, unit = "pt"),
    panel.spacing = grid::unit(4, "pt"),
    plot.margin = margin(4, 5, 4, 4, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

prediction_error_output = file.path(
  data_dir, paste0("pred_error_", design, ".pdf")
)
individual_bias_output = file.path(
  data_dir, paste0("indiv_bias_", design, ".pdf")
)

ggsave(
  prediction_error_output, prediction_error_plot, device = "pdf",
  width = 7.25, height = 3.25, units = "in", useDingbats = FALSE
)
ggsave(
  individual_bias_output, individual_bias_plot, device = "pdf",
  width = 7.5, height = 3.5, units = "in", useDingbats = FALSE
)

message("Created: ", prediction_error_output)
message("Created: ", individual_bias_output)
