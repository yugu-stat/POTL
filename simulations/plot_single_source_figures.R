# Create the two main-paper figures for a five-scenario simulation design.
#
# The script reads saved RData files from its own directory by default and
# writes exactly two figures to that directory:
#   pred_error_<design>.pdf  : L2D and D_tau
#   indiv_bias_<design>.pdf  : mean signed individual survival bias
#
# Built-in designs:
#   Rscript plot_simulation_figures.R single_source
#   Rscript plot_simulation_figures.R single_source_diffcovar
#
# Optional arguments:
#   Rscript plot_simulation_figures.R <design> <data_dir>
#   Rscript plot_simulation_figures.R <design> <data_dir> \
#     <scalar_pattern> <bias_pattern>
#
# Each pattern must contain one integer conversion (normally "%d") for the
# scenario number. The optional patterns make the program reusable when a
# design uses different input filenames. The design term is used unchanged in
# the two output filenames and may contain letters, numbers, underscores, and
# hyphens.

suppressPackageStartupMessages(library(ggplot2))

full_args = commandArgs(trailingOnly = FALSE)
script_arg = grep("^--file=", full_args, value = TRUE)
if(length(script_arg) != 1L) {
  stop("Run this plotting program with Rscript.")
}
script_file = normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir = dirname(script_file)

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 1L || length(args) > 4L || length(args) == 3L) {
  stop(
    "Usage: Rscript plot_simulation_figures.R <design> [data_dir] ",
    "[scalar_pattern bias_pattern]"
  )
}

design = args[1]
if(!grepl("^[A-Za-z0-9_-]+$", design)) {
  stop("design may contain only letters, numbers, underscores, and hyphens.")
}

data_dir = if(length(args) >= 2L) args[2] else script_dir
data_dir = normalizePath(data_dir, mustWork = TRUE)

design_files = list(
  single_source = list(
    scalar = "simres_sc%d.RData",
    bias = "individual_survival_bias_single_source_sc%d.RData"
  ),
  single_source_diffcovar = list(
    scalar = "simres_sc%d.RData",
    bias = "individual_survival_bias_single_source_diffcovar_sc%d.RData"
  )
)

if(length(args) == 4L) {
  scalar_pattern = args[3]
  bias_pattern = args[4]
} else {
  if(!design %in% names(design_files)) {
    stop(
      "Unknown built-in design '", design, "'. Supply scalar_pattern and ",
      "bias_pattern explicitly for a new design."
    )
  }
  scalar_pattern = design_files[[design]]$scalar
  bias_pattern = design_files[[design]]$bias
}

scenario_path = function(pattern, scenario) {
  filename = tryCatch(
    sprintf(pattern, scenario),
    error = function(e) {
      stop("Invalid filename pattern '", pattern, "': ", conditionMessage(e))
    }
  )
  file.path(data_dir, filename)
}

scalar_paths = vapply(1:5, function(sc) scenario_path(scalar_pattern, sc), "")
bias_paths = vapply(1:5, function(sc) scenario_path(bias_pattern, sc), "")
if(anyDuplicated(scalar_paths)) {
  stop("scalar_pattern must generate a different filename for each scenario.")
}
if(anyDuplicated(bias_paths)) {
  stop("bias_pattern must generate a different filename for each scenario.")
}

scenarios = 1:5
scenario_levels = paste0("SC", scenarios)
method_levels = c("POTL", "Target-only", "TransCox", "CoxTL", "Pooled")
method_colors = c(
  "POTL" = "blue",
  "Target-only" = "black",
  "TransCox" = "red",
  "CoxTL" = "green",
  "Pooled" = "purple"
)
method_shapes = c(
  "POTL" = 16,
  "Target-only" = 17,
  "TransCox" = 15,
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

read_scalar_scenario = function(sc) {
  path = scenario_path(scalar_pattern, sc)
  allres = load_saved_object(path, "allres")
  if(!is.list(allres) || length(allres) == 0L) {
    stop("allres must be a nonempty list in ", path)
  }

  pieces = lapply(seq_along(allres), function(replication) {
    value = allres[[replication]]
    if(!is.matrix(value) || !is.numeric(value) ||
       !identical(dim(value), c(length(method_levels), 6L))) {
      stop(
        "Expected a 5 x 6 numeric result matrix in ", path,
        ", replication ", replication, "."
      )
    }
    method_id = as.integer(value[,1])
    if(!identical(method_id, seq_along(method_levels)) ||
       any(value[,1] != method_id)) {
      stop(
        "Expected method IDs 1:5 in ", path,
        ", replication ", replication, "."
      )
    }
    errors = value[,c(2, 3), drop = FALSE]
    if(any(!is.finite(errors)) || any(errors < 0)) {
      stop(
        "L2D and D_tau must be finite and nonnegative in ", path,
        ", replication ", replication, "."
      )
    }
    data.frame(
      scenario = paste0("SC", sc),
      replication = replication,
      method = method_levels[method_id],
      L2D = errors[,1],
      Dtau = errors[,2],
      stringsAsFactors = FALSE
    )
  })

  list(data = do.call(rbind, pieces), n_rep = length(allres))
}

scalar_loaded = stats::setNames(
  lapply(scenarios, read_scalar_scenario), scenario_levels
)
scalar_replications = vapply(scalar_loaded, `[[`, integer(1), "n_rep")
if(length(unique(scalar_replications)) != 1L) {
  stop(
    "Scalar-result replication counts differ across scenarios: ",
    paste(names(scalar_replications), scalar_replications, sep = "=", collapse = ", ")
  )
}

scalar_data = do.call(rbind, lapply(scalar_loaded, `[[`, "data"))
scalar_data$scenario = factor(scalar_data$scenario, levels = scenario_levels)
scalar_data$method = factor(scalar_data$method, levels = method_levels)

summarize_scalar_metric = function(variable, label) {
  rows = vector("list", length(scenario_levels)*length(method_levels))
  k = 0L
  for(sc in scenario_levels) {
    for(method in method_levels) {
      k = k+1L
      keep = scalar_data$scenario == sc & scalar_data$method == method
      values = scalar_data[[variable]][keep]
      if(length(values) != scalar_replications[[sc]]) {
        stop("Incomplete scalar results for ", sc, ", ", method, ".")
      }
      center = stats::median(values)
      spread = stats::mad(values, center = center, constant = 1.4826)
      rows[[k]] = data.frame(
        scenario = sc,
        method = method,
        metric = label,
        median = center,
        mad = spread,
        ymin = max(0, center-spread),
        ymax = center+spread,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

scalar_summary = rbind(
  summarize_scalar_metric("L2D", "L[2]*D"),
  summarize_scalar_metric("Dtau", "D[tau]")
)
scalar_summary$scenario = factor(
  scalar_summary$scenario, levels = scenario_levels
)
scalar_summary$method = factor(scalar_summary$method, levels = method_levels)
scalar_summary$metric = factor(
  scalar_summary$metric, levels = c("L[2]*D", "D[tau]")
)

scalar_ymax = ceiling(100*max(scalar_summary$ymax))/100
scalar_dodge = position_dodge(width = 0.68)
prediction_error_plot = ggplot(
  scalar_summary,
  aes(x = scenario, y = median, color = method, shape = method)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = scalar_dodge, width = 0.12, linewidth = 0.45
  ) +
  geom_point(position = scalar_dodge, size = 2.25, stroke = 0.45) +
  facet_wrap(~metric, nrow = 1, labeller = label_parsed) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, scalar_ymax),
    breaks = pretty(c(0, scalar_ymax), n = 5),
    expand = expansion(mult = c(0, 0.025))
  ) +
  labs(
    x = "Simulation scenario", y = "Prediction error",
    color = NULL, shape = NULL
  ) +
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
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = -4, unit = "pt"),
    legend.key.width = grid::unit(12, "pt"),
    legend.spacing.x = grid::unit(3, "pt"),
    plot.margin = margin(5, 6, -2, 5, unit = "pt")
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  )

# ---------------------------------------------------------------------------
# Individual survival bias. For subject i, method m, and prediction time t,
# mean_bias[i,t,m] = mean_r{S_hat_rim(t) - S_0(t | X_i)}. The figure shows
# the distribution of this signed mean bias across validation subjects.
# ---------------------------------------------------------------------------

safe_bandwidth = function(values) {
  if(length(unique(values)) < 2L) {
    return(NA_real_)
  }
  bandwidth = suppressWarnings(stats::bw.nrd0(values))
  if(is.finite(bandwidth) && bandwidth > 0) bandwidth else NA_real_
}

read_bias_scenario = function(sc) {
  path = scenario_path(bias_pattern, sc)
  result = load_saved_object(path, "individual_bias")
  required = c(
    "mean_bias", "n_rep", "methods", "prediction_times", "validation_data",
    "requested_replications", "successful_replications"
  )
  if(!is.list(result) || !all(required %in% names(result))) {
    stop("Unexpected individual_bias structure in ", path)
  }
  if(!identical(as.character(result$methods), method_levels)) {
    stop("Unexpected method order in ", path)
  }

  prediction_times = as.numeric(result$prediction_times)
  if(length(prediction_times) != 4L || any(!is.finite(prediction_times)) ||
     any(diff(prediction_times) <= 0)) {
    stop("Expected four strictly increasing prediction times in ", path)
  }
  if(is.null(nrow(result$validation_data)) || nrow(result$validation_data) < 1L) {
    stop("validation_data must contain at least one subject in ", path)
  }

  expected_dim = as.integer(c(
    nrow(result$validation_data), length(prediction_times), length(method_levels)
  ))
  if(!identical(dim(result$mean_bias), expected_dim) ||
     !identical(dim(result$n_rep), expected_dim)) {
    stop("Unexpected individual-bias dimensions in ", path)
  }
  if(!is.numeric(result$mean_bias) || any(!is.finite(result$mean_bias))) {
    stop("mean_bias must contain only finite numeric values in ", path)
  }
  if(any(abs(result$mean_bias) > 1+sqrt(.Machine$double.eps))) {
    stop("mean_bias lies outside the survival-probability scale in ", path)
  }

  requested = as.integer(result$requested_replications)
  successful = as.integer(result$successful_replications)
  if(length(requested) != 1L || length(successful) != 1L ||
     is.na(requested) || is.na(successful) || requested < 1L ||
     successful != requested) {
    stop("Not all requested individual-bias replications succeeded in ", path)
  }
  if(anyNA(result$n_rep) || any(result$n_rep != successful)) {
    stop("Some individual-bias cells do not use all replications in ", path)
  }

  if("data" %in% names(result) && "scenario" %in% names(result$data)) {
    saved_scenario = unique(result$data$scenario[!is.na(result$data$scenario)])
    if(length(saved_scenario) != 1L || as.integer(saved_scenario) != sc) {
      stop("Saved scenario metadata does not match the filename in ", path)
    }
  }

  bias_pp = 100*result$mean_bias
  index = expand.grid(
    subject = seq_len(expected_dim[1]),
    time_index = seq_len(expected_dim[2]),
    method_index = seq_len(expected_dim[3]),
    KEEP.OUT.ATTRS = FALSE
  )
  plot_data = data.frame(
    scenario = paste0("SC", sc),
    prediction_time = prediction_times[index$time_index],
    method = method_levels[index$method_index],
    bias_pp = as.vector(bias_pp),
    stringsAsFactors = FALSE
  )

  summaries = vector(
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
      summaries[[k]] = data.frame(
        scenario = paste0("SC", sc),
        prediction_time = prediction_times[time_index],
        method = method_levels[method_index],
        q05 = q[1], q25 = q[2], median = q[3], q75 = q[4], q95 = q[5],
        mean = mean(values),
        bandwidth = safe_bandwidth(values),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    plot_data = plot_data,
    summary = do.call(rbind, summaries),
    validation_data = result$validation_data,
    prediction_times = prediction_times,
    requested_replications = requested,
    successful_replications = successful
  )
}

bias_loaded = stats::setNames(
  lapply(scenarios, read_bias_scenario), scenario_levels
)

reference_times = bias_loaded[[1]]$prediction_times
same_times = vapply(
  bias_loaded,
  function(x) identical(x$prediction_times, reference_times),
  logical(1)
)
if(!all(same_times)) {
  stop("Prediction times differ across scenarios.")
}

reference_validation = bias_loaded[[1]]$validation_data
same_validation = vapply(
  bias_loaded,
  function(x) identical(x$validation_data, reference_validation),
  logical(1)
)
if(!all(same_validation)) {
  stop("The five scenarios do not use an identical validation cohort.")
}

bias_successful = vapply(
  bias_loaded, `[[`, integer(1), "successful_replications"
)
if(!identical(unname(bias_successful), unname(scalar_replications))) {
  stop("Scalar and individual-bias replication counts do not agree.")
}

bias_data = do.call(rbind, lapply(bias_loaded, `[[`, "plot_data"))
bias_summary = do.call(rbind, lapply(bias_loaded, `[[`, "summary"))

time_labels = paste0(
  "t = ", format(reference_times, trim = TRUE, scientific = FALSE)
)
bias_data$scenario = factor(bias_data$scenario, levels = scenario_levels)
bias_data$prediction_time = factor(
  bias_data$prediction_time, levels = reference_times, labels = time_labels
)
bias_data$method = factor(bias_data$method, levels = method_levels)
bias_summary$scenario = factor(
  bias_summary$scenario, levels = scenario_levels
)
bias_summary$prediction_time = factor(
  bias_summary$prediction_time, levels = reference_times, labels = time_labels
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

# Absolute values are used here only to make each signed-bias scale symmetric
# around zero; they are not averaged or reported as an error metric.
observed_row_max = tapply(
  abs(bias_data$bias_pp), bias_data$scenario, max, na.rm = TRUE
)
row_limits = vapply(observed_row_max, nice_symmetric_limit, numeric(1))

row_bandwidths = vapply(scenario_levels, function(sc) {
  candidates = bias_summary$bandwidth[bias_summary$scenario == sc]
  candidates = candidates[is.finite(candidates) & candidates > 0]
  if(length(candidates) > 0L) {
    stats::median(candidates)
  } else {
    row_limits[[sc]]/50
  }
}, numeric(1))

limit_data = expand.grid(
  scenario = scenario_levels,
  prediction_time = time_labels,
  endpoint = c(-1, 1),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
limit_data$scenario = factor(limit_data$scenario, levels = scenario_levels)
limit_data$prediction_time = factor(
  limit_data$prediction_time, levels = time_labels
)
limit_data$method = factor("POTL", levels = method_levels)
limit_data$bias_pp =
  row_limits[as.character(limit_data$scenario)]*limit_data$endpoint

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
  )

for(sc in scenario_levels) {
  individual_bias_plot = individual_bias_plot + geom_violin(
    data = bias_data[bias_data$scenario == sc,],
    bw = row_bandwidths[[sc]], scale = "width", trim = TRUE,
    linewidth = 0.25, alpha = 0.78, color = "grey20"
  )
}

individual_bias_plot = individual_bias_plot +
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
  facet_grid(scenario ~ prediction_time, scales = "free_y") +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Individual-level mean bias in predicted survival (\u00d7 100)",
    fill = NULL
  ) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
    panel.border = element_rect(linewidth = 0.4),
    strip.background = element_rect(fill = "grey93", linewidth = 0.4),
    strip.text.x = element_text(size = 8.5),
    strip.text.y = element_text(size = 8.5),
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
  width = 7.5, height = 9.5, units = "in", useDingbats = FALSE
)

message("Created: ", prediction_error_output)
message("Created: ", individual_bias_output)
