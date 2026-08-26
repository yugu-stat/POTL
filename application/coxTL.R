# CoxTL implementation following YUYING-LU/CoxTL at commit
# cb2b95b80eda0a2fcc1ba3349a3e03158d42e11a.

# Exponential-tilt density-ratio objective used by CoxTL.
density_opt <- function(Xs, Xt, lam) {
  function(theta_0) {
    theta_0 <- as.matrix(theta_0, ncol = 1)
    mean(exp(Xs %*% theta_0)) - mean(Xt %*% theta_0) +
      lam * norm(theta_0, type = "2")^2
  }
}

.CoxTL_formula <- function(covariates) {
  # rms::cph recognizes strat() as its stratification special. Keeping the
  # alias in the formula environment avoids requiring rms to be attached.
  strat <- rms::strat
  stats::as.formula(
    paste(
      "survival::Surv(time, status) ~",
      paste(c(covariates, "strat(R)"), collapse = " + ")
    ),
    env = environment()
  )
}

# Fit CoxTL using source/target-specific risk sets and baseline hazards.
run_CoxTL <- function(data_t, data_s, weights = NULL, di, folds_num = 5,
                      lam_set = c(0.001, 0.005, 0.01, 0.05, 0.1, 1:10)) {
  data_t <- as.data.frame(data_t)
  data_s <- as.data.frame(data_s)
  num_t <- nrow(data_t)
  num_s <- nrow(data_s)

  if (di < 1 || di > ncol(data_t)) {
    stop("di must identify the leading CoxTL covariate columns")
  }
  covariates <- colnames(data_t)[seq_len(di)]
  if (!identical(covariates, colnames(data_s)[seq_len(di)])) {
    stop("target and source CoxTL covariates must have the same names and order")
  }
  required <- c("time", "status", "R")
  if (!all(required %in% colnames(data_t)) ||
      !all(required %in% colnames(data_s))) {
    stop("data_t and data_s must contain time, status, and R")
  }
  if (!is.factor(data_t$R) || !is.factor(data_s$R) ||
      !identical(levels(data_t$R), levels(data_s$R))) {
    stop("R must be a factor with common levels in data_t and data_s")
  }
  if (is.null(weights)) {
    weights <- rep(1, num_s)
  }
  weights <- as.numeric(weights)
  if (length(weights) != num_s || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    stop("source weights must be finite, positive, and match nrow(data_s)")
  }

  size_ratio <- num_t / num_s
  strat <- rms::strat
  cox_formula <- .CoxTL_formula(covariates)
  split_set <- sample(seq_len(folds_num), num_t, replace = TRUE)
  cv_record <- matrix(0, nrow = length(lam_set), ncol = folds_num)

  for (lam_idx in seq_along(lam_set)) {
    lambda <- size_ratio * lam_set[lam_idx]
    for (k in seq_len(folds_num)) {
      data_t1 <- data_t[split_set == k, , drop = FALSE]
      data_t2 <- data_t[split_set != k, , drop = FALSE]
      data_ts <- rbind(data_s, data_t2)
      cox_train <- rms::cph(
        cox_formula,
        data = data_ts,
        weights = c(lambda * weights, rep(1, nrow(data_t2))),
        x = TRUE,
        y = TRUE,
        surv = TRUE
      )
      cv_record[lam_idx, k] <- 1 - Hmisc::rcorr.cens(
        exp(as.matrix(data_t1[, seq_len(di), drop = FALSE]) %*%
              stats::coef(cox_train)),
        survival::Surv(data_t1$time, data_t1$status)
      )["C Index"]
    }
  }

  lam_best <- size_ratio * lam_set[which.max(rowMeans(cv_record))]
  rms::cph(
    cox_formula,
    data = rbind(data_s, data_t),
    weights = c(lam_best * weights, rep(1, num_t)),
    x = TRUE,
    y = TRUE,
    surv = TRUE
  )
}

# Estimate the source-to-target covariate density ratio as in the CoxTL
# tutorial (ridge penalty fixed at one).
CoxTL_density_weights <- function(data_t, data_s, di, lam = 1) {
  Xs_new <- as.matrix(cbind(
    1,
    data_s[, seq_len(di), drop = FALSE]
  ))
  Xt_new <- as.matrix(cbind(
    1,
    data_t[, seq_len(di), drop = FALSE]
  ))
  w_opt <- stats::nlm(
    density_opt(Xs_new, Xt_new, lam),
    rep(0, di + 1),
    iterlim = 500
  )
  weights <- as.numeric(exp(Xs_new %*% matrix(w_opt$estimate, ncol = 1)))
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("CoxTL density-ratio estimation produced invalid source weights")
  }
  weights
}

# A stratified CoxTL fit contains one baseline hazard per population. The
# manuscript metrics and predictions require the target-population baseline.
CoxTL_target_basehaz <- function(cox_tl, centered = FALSE) {
  cumhaz <- survival::basehaz(cox_tl, centered = centered)
  if (!"strata" %in% colnames(cumhaz)) {
    stop("the CoxTL fit is not stratified by R")
  }
  is_target <- as.character(cumhaz$strata) == "R=0"
  if (!any(is_target)) {
    stop("the CoxTL fit does not contain the target stratum R=0")
  }
  cumhaz[is_target, c("hazard", "time"), drop = FALSE]
}
