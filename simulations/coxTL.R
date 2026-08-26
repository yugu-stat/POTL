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

# Convert a simulation matrix to the data-frame layout expected by CoxTL.
# In particular, adding R with `$<-` directly to a matrix flattens the matrix
# into a list and makes the density-ratio step fail with "incorrect number of
# dimensions".
CoxTL_prepare_data <- function(data, di, population,
                               covariates = paste0("X", seq_len(di)),
                               strata_levels = c(0, 1)) {
  if (length(di) != 1L || !is.numeric(di) || !is.finite(di) ||
      di < 1 || di != as.integer(di)) {
    stop("di must be a positive integer")
  }
  di <- as.integer(di)
  data <- as.data.frame(data)
  if (nrow(data) == 0L) {
    stop("CoxTL data must contain at least one row")
  }
  if (ncol(data) != di + 2L) {
    stop("CoxTL data must contain di covariates followed by time and status")
  }
  if (length(covariates) != di || anyDuplicated(covariates)) {
    stop("covariates must contain di unique names")
  }
  if (!identical(as.character(strata_levels), c("0", "1"))) {
    stop("single-source CoxTL strata_levels must be c(0, 1)")
  }
  if (length(population) != 1L || !population %in% strata_levels) {
    stop("population must be one of strata_levels")
  }

  colnames(data) <- c(covariates, "time", "status")
  if (!is.numeric(data$time) || any(!is.finite(data$time)) ||
      any(data$time < 0)) {
    stop("CoxTL time must be finite and nonnegative")
  }
  if (!(is.numeric(data$status) || is.logical(data$status)) ||
      anyNA(data$status) || any(!data$status %in% c(0, 1))) {
    stop("CoxTL status must contain only 0 and 1")
  }
  data$status <- as.numeric(data$status)
  data$R <- factor(population, levels = strata_levels)
  data
}

.CoxTL_density_design <- function(data_t, data_s, di,
                                  standardize = FALSE) {
  if (length(di) != 1L || !is.numeric(di) || !is.finite(di) ||
      di < 1 || di != as.integer(di)) {
    stop("di must be a positive integer")
  }
  di <- as.integer(di)
  data_t <- as.data.frame(data_t)
  data_s <- as.data.frame(data_s)
  if (nrow(data_t) == 0L || nrow(data_s) == 0L) {
    stop("target and source density-ratio data must be nonempty")
  }
  if (di < 1L || di > ncol(data_t) || di > ncol(data_s)) {
    stop("di must identify the leading CoxTL covariate columns")
  }

  target_x <- data_t[, seq_len(di), drop = FALSE]
  source_x <- data_s[, seq_len(di), drop = FALSE]
  is_numeric <- vapply(c(target_x, source_x), is.numeric, logical(1))
  if (!all(is_numeric)) {
    stop("CoxTL density-ratio covariates must be numeric")
  }
  Xt <- as.matrix(target_x)
  Xs <- as.matrix(source_x)
  if (any(!is.finite(Xt)) || any(!is.finite(Xs))) {
    stop("CoxTL density-ratio covariates must be finite")
  }

  center <- rep(0, di)
  scale <- rep(1, di)
  if (standardize) {
    combined_x <- rbind(Xs, Xt)
    center <- colMeans(combined_x)
    scale <- apply(combined_x, 2L, stats::sd)
    scale[!is.finite(scale) | scale < sqrt(.Machine$double.eps)] <- 1
    Xs <- sweep(sweep(Xs, 2L, center, "-"), 2L, scale, "/")
    Xt <- sweep(sweep(Xt, 2L, center, "-"), 2L, scale, "/")
  }

  list(
    Xs = cbind(`(Intercept)` = 1, Xs),
    Xt = cbind(`(Intercept)` = 1, Xt),
    center = center,
    scale = scale
  )
}

.CoxTL_fit_density <- function(Xs, Xt, lam) {
  fit <- suppressWarnings(stats::nlm(
    density_opt(Xs, Xt, lam),
    rep(0, ncol(Xs)),
    iterlim = 500
  ))
  if (!fit$code %in% c(1L, 2L) || any(!is.finite(fit$estimate))) {
    stop("density-ratio optimizer did not converge")
  }
  fit$estimate
}

.CoxTL_density_validation_loss <- function(Xs, Xt, theta) {
  source_lp <- as.numeric(Xs %*% theta)
  target_lp <- as.numeric(Xt %*% theta)
  max_safe_lp <- log(.Machine$double.xmax) - log(max(1L, nrow(Xs)))
  if (any(!is.finite(source_lp)) || any(!is.finite(target_lp)) ||
      max(source_lp) > max_safe_lp) {
    return(NA_real_)
  }
  loss <- mean(exp(source_lp)) - mean(target_lp)
  if (is.finite(loss)) loss else NA_real_
}

.CoxTL_normalize_density_weights <- function(Xs, theta) {
  log_weights <- as.numeric(Xs %*% theta)
  if (any(!is.finite(log_weights))) {
    stop("density-ratio estimation produced non-finite log weights")
  }
  # Centering before exponentiation avoids overflow. The subsequent mean-one
  # normalization also separates density-ratio shape from CoxTL's source-scale
  # tuning parameter.
  weights <- exp(log_weights - max(log_weights))
  if (any(!is.finite(weights)) || any(weights <= 0) || mean(weights) <= 0) {
    stop("density-ratio estimation produced invalid source weights")
  }
  weights / mean(weights)
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
  if (num_t == 0L || num_s == 0L) {
    stop("target and source CoxTL data must be nonempty")
  }

  if (length(di) != 1L || !is.numeric(di) || !is.finite(di) ||
      di < 1 || di != as.integer(di) || di > ncol(data_t)) {
    stop("di must identify the leading CoxTL covariate columns")
  }
  di <- as.integer(di)
  covariates <- colnames(data_t)[seq_len(di)]
  if (!identical(covariates, colnames(data_s)[seq_len(di)])) {
    stop("target and source CoxTL covariates must have the same names and order")
  }
  required <- c("time", "status", "R")
  if (!all(required %in% colnames(data_t)) ||
      !all(required %in% colnames(data_s))) {
    stop("data_t and data_s must contain time, status, and R")
  }
  target_covariates <- data_t[, seq_len(di), drop = FALSE]
  source_covariates <- data_s[, seq_len(di), drop = FALSE]
  cox_covariates <- c(target_covariates, source_covariates)
  if (!all(vapply(cox_covariates, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(target_covariates))) ||
      any(!is.finite(as.matrix(source_covariates)))) {
    stop("CoxTL covariates must be numeric and finite")
  }
  valid_outcome <- function(data) {
    is.numeric(data$time) && all(is.finite(data$time)) &&
      all(data$time >= 0) &&
      (is.numeric(data$status) || is.logical(data$status)) &&
      !anyNA(data$status) && all(data$status %in% c(0, 1))
  }
  if (!valid_outcome(data_t) || !valid_outcome(data_s)) {
    stop("CoxTL time must be finite/nonnegative and status must be binary")
  }
  if (!is.factor(data_t$R) || !is.factor(data_s$R) ||
      !identical(levels(data_t$R), c("0", "1")) ||
      !identical(levels(data_t$R), levels(data_s$R)) ||
      anyNA(data_t$R) || anyNA(data_s$R) ||
      any(as.character(data_t$R) != "0") ||
      any(as.character(data_s$R) != "1")) {
    stop(
      "R must use common levels c('0', '1'), with target R=0 and source R=1"
    )
  }
  if (is.null(weights)) {
    weights <- rep(1, num_s)
  }
  weights <- as.numeric(weights)
  if (length(weights) != num_s || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    stop("source weights must be finite, positive, and match nrow(data_s)")
  }
  lam_set <- sort(unique(as.numeric(lam_set)))
  if (!length(lam_set) || any(!is.finite(lam_set)) || any(lam_set <= 0)) {
    stop("lam_set must contain finite, positive values")
  }
  if (length(folds_num) != 1L || !is.numeric(folds_num) ||
      !is.finite(folds_num) || folds_num != as.integer(folds_num) ||
      folds_num < 2L || folds_num > num_t) {
    stop("folds_num must be between 2 and the target sample size")
  }
  folds_num <- as.integer(folds_num)

  size_ratio <- num_t / num_s
  strat <- rms::strat
  cox_formula <- .CoxTL_formula(covariates)
  # Balanced folds ensure every validation fold is nonempty. An individual
  # numerical failure is recorded as NA rather than discarding the replicate.
  split_set <- sample(rep(seq_len(folds_num), length.out = num_t))
  cv_record <- matrix(
    NA_real_, nrow = length(lam_set), ncol = folds_num,
    dimnames = list(as.character(lam_set), paste0("fold", seq_len(folds_num)))
  )

  for (lam_idx in seq_along(lam_set)) {
    lambda <- size_ratio * lam_set[lam_idx]
    for (k in seq_len(folds_num)) {
      data_t1 <- data_t[split_set == k, , drop = FALSE]
      data_t2 <- data_t[split_set != k, , drop = FALSE]
      data_ts <- rbind(data_s, data_t2)
      cv_record[lam_idx, k] <- tryCatch({
        cox_train <- rms::cph(
          cox_formula,
          data = data_ts,
          weights = c(lambda * weights, rep(1, nrow(data_t2))),
          x = TRUE,
          y = TRUE,
          surv = TRUE
        )
        beta <- stats::coef(cox_train)
        if (length(beta) != di || any(!is.finite(beta))) {
          stop("CoxTL cross-validation fit has invalid coefficients")
        }
        linear_predictor <- as.numeric(
          as.matrix(data_t1[, seq_len(di), drop = FALSE]) %*% beta
        )
        # rcorr.cens treats larger predictions as longer survival. Negating
        # the Cox risk score gives the desired orientation without exp(lp),
        # which can overflow for an otherwise usable fit.
        score <- Hmisc::rcorr.cens(
          -linear_predictor,
          survival::Surv(data_t1$time, data_t1$status)
        )["C Index"]
        if (length(score) != 1L || !is.finite(score)) NA_real_ else score
      }, error = function(e) {
        NA_real_
      })
    }
  }

  valid_folds <- rowSums(is.finite(cv_record))
  min_valid_folds <- max(2L, ceiling(folds_num / 2))
  cv_score <- vapply(seq_len(nrow(cv_record)), function(i) {
    keep <- is.finite(cv_record[i, ])
    if (sum(keep) < min_valid_folds) NA_real_ else mean(cv_record[i, keep])
  }, numeric(1))
  if (!any(is.finite(cv_score))) {
    warning(
      "CoxTL source-scale cross-validation had no finite candidate; ",
      "falling back to the candidate closest to one"
    )
    first_idx <- which.min(abs(log(lam_set)))
    fit_order <- first_idx
  } else {
    fit_order <- order(cv_score, decreasing = TRUE, na.last = NA)
  }
  fit_order <- unique(c(
    fit_order,
    order(abs(log(lam_set)), na.last = NA),
    seq_along(lam_set)
  ))

  final_data <- rbind(data_s, data_t)
  final_fit <- NULL
  selected_idx <- NA_integer_
  for (lam_idx in fit_order) {
    lambda <- size_ratio * lam_set[lam_idx]
    candidate_fit <- tryCatch(
      rms::cph(
        cox_formula,
        data = final_data,
        weights = c(lambda * weights, rep(1, num_t)),
        x = TRUE,
        y = TRUE,
        surv = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(candidate_fit)) {
      beta <- stats::coef(candidate_fit)
      if (length(beta) == di && all(is.finite(beta))) {
        final_fit <- candidate_fit
        selected_idx <- lam_idx
        break
      }
    }
  }
  if (is.null(final_fit)) {
    stop("CoxTL failed for every source-scale tuning candidate")
  }
  final_fit$coxtl_tuning <- list(
    source_scale = size_ratio * lam_set[selected_idx],
    source_scale_multiplier = lam_set[selected_idx],
    cv_score = stats::setNames(cv_score, lam_set),
    valid_folds = stats::setNames(valid_folds, lam_set)
  )
  final_fit
}

# Estimate the source-to-target covariate density ratio. Supplying lam_set
# enables balanced K-fold tuning by held-out, unpenalized density loss; the
# scalar lam = 1 route retains the CoxTL tutorial's penalty default while also
# applying the numerical checks and mean-one normalization below.
CoxTL_density_weights <- function(data_t, data_s, di, lam = 1,
                                  lam_set = NULL, folds_num = 5,
                                  one_se = TRUE,
                                  standardize = !is.null(lam_set),
                                  min_ess_ratio = 0.1) {
  if (length(standardize) != 1L || !is.logical(standardize) ||
      is.na(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }
  if (length(one_se) != 1L || !is.logical(one_se) || is.na(one_se)) {
    stop("one_se must be TRUE or FALSE")
  }
  data_t <- as.data.frame(data_t)
  data_s <- as.data.frame(data_s)
  design <- .CoxTL_density_design(
    data_t, data_s, di, standardize = standardize
  )
  Xs <- design$Xs
  Xt <- design$Xt
  candidates <- if (is.null(lam_set)) lam else lam_set
  candidates <- sort(unique(as.numeric(candidates)))
  if (!length(candidates) || any(!is.finite(candidates)) ||
      any(candidates <= 0)) {
    stop("density-ratio tuning values must be finite and positive")
  }
  if (length(min_ess_ratio) != 1L || !is.finite(min_ess_ratio) ||
      min_ess_ratio <= 0 || min_ess_ratio > 1) {
    stop("min_ess_ratio must be in (0, 1]")
  }

  cv_loss <- cv_se <- rep(NA_real_, length(candidates))
  selected_idx <- 1L
  if (length(candidates) > 1L) {
    if (length(folds_num) != 1L || !is.numeric(folds_num) ||
        !is.finite(folds_num) || folds_num != as.integer(folds_num) ||
        folds_num < 2L) {
      stop("folds_num must be an integer allowing at least two density-ratio folds")
    }
    folds_num <- min(as.integer(folds_num), nrow(Xs), nrow(Xt))
    if (folds_num < 2L) {
      stop("target and source data must each allow at least two density-ratio folds")
    }
    source_fold <- sample(rep(seq_len(folds_num), length.out = nrow(Xs)))
    target_fold <- sample(rep(seq_len(folds_num), length.out = nrow(Xt)))
    fold_loss <- matrix(
      NA_real_, nrow = length(candidates), ncol = folds_num,
      dimnames = list(as.character(candidates), paste0("fold", seq_len(folds_num)))
    )

    for (lam_idx in seq_along(candidates)) {
      for (k in seq_len(folds_num)) {
        fold_loss[lam_idx, k] <- tryCatch({
          theta <- .CoxTL_fit_density(
            Xs[source_fold != k, , drop = FALSE],
            Xt[target_fold != k, , drop = FALSE],
            candidates[lam_idx]
          )
          .CoxTL_density_validation_loss(
            Xs[source_fold == k, , drop = FALSE],
            Xt[target_fold == k, , drop = FALSE],
            theta
          )
        }, error = function(e) {
          NA_real_
        })
      }
    }

    cv_loss <- vapply(seq_len(nrow(fold_loss)), function(i) {
      keep <- is.finite(fold_loss[i, ])
      if (sum(keep) < 2L) NA_real_ else mean(fold_loss[i, keep])
    }, numeric(1))
    cv_se <- vapply(seq_len(nrow(fold_loss)), function(i) {
      keep <- is.finite(fold_loss[i, ])
      if (sum(keep) < 2L) NA_real_ else
        stats::sd(fold_loss[i, keep]) / sqrt(sum(keep))
    }, numeric(1))

    if (any(is.finite(cv_loss))) {
      min_idx <- which.min(replace(cv_loss, !is.finite(cv_loss), Inf))
      selected_idx <- min_idx
      if (isTRUE(one_se) && is.finite(cv_se[min_idx])) {
        eligible <- which(
          is.finite(cv_loss) & cv_loss <= cv_loss[min_idx] + cv_se[min_idx]
        )
        selected_idx <- eligible[which.max(candidates[eligible])]
      }
    } else {
      # The most regularized candidate is the safest starting point when CV
      # contains no finite validation loss.
      selected_idx <- which.max(candidates)
    }
  }

  # If the selected fit collapses the effective source sample, move toward
  # stronger regularization before considering less stable candidates.
  fit_order <- unique(c(
    selected_idx,
    which(candidates > candidates[selected_idx]),
    rev(which(candidates < candidates[selected_idx]))
  ))
  weights <- NULL
  used_idx <- NA_integer_
  effective_sample_size <- NA_real_
  failure_messages <- character(0)
  for (lam_idx in fit_order) {
    candidate_weights <- tryCatch({
      theta <- .CoxTL_fit_density(Xs, Xt, candidates[lam_idx])
      .CoxTL_normalize_density_weights(Xs, theta)
    }, error = function(e) {
      failure_messages <<- c(failure_messages, conditionMessage(e))
      NULL
    })
    if (is.null(candidate_weights)) next
    candidate_ess <- sum(candidate_weights)^2 / sum(candidate_weights^2)
    if (!is.finite(candidate_ess) ||
        candidate_ess < min_ess_ratio * nrow(Xs)) {
      failure_messages <- c(
        failure_messages,
        sprintf(
          "lambda %g produced ESS %.1f below the minimum %.1f",
          candidates[lam_idx], candidate_ess, min_ess_ratio * nrow(Xs)
        )
      )
      next
    }
    weights <- candidate_weights
    effective_sample_size <- candidate_ess
    used_idx <- lam_idx
    break
  }

  fallback <- is.null(weights)
  if (fallback) {
    warning(
      "All CoxTL density-ratio candidates were unstable; using uniform ",
      "source weights. ", paste(unique(failure_messages), collapse = "; ")
    )
    weights <- rep(1, nrow(Xs))
    effective_sample_size <- as.double(nrow(Xs))
  }
  attr(weights, "lambda") <- if (fallback) NA_real_ else candidates[used_idx]
  attr(weights, "cv_loss") <- stats::setNames(cv_loss, candidates)
  attr(weights, "cv_se") <- stats::setNames(cv_se, candidates)
  attr(weights, "effective_sample_size") <- effective_sample_size
  attr(weights, "fallback") <- fallback
  attr(weights, "standardize_center") <- design$center
  attr(weights, "standardize_scale") <- design$scale
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

## Local multi-source extension of CoxTL.
run_CoxTL_ss<-function(sc, data_t,data_s_list,p,folds_num=5,
                       lam_set = 10^(seq(-4,1,1)),
                       nu_set=10^(seq(-4,1,1))){
  data_t = as.data.frame(data_t)
  data_s = as.data.frame(data_s_list[[sc]])
  num_t<-dim(data_t)[1]
  num_s<-dim(data_s)[1]
  data_t$R = factor(0, levels = c(0, 1))
  data_s$R = factor(1, levels = c(0, 1))
  strat <- rms::strat
  cox_formula = .CoxTL_formula(colnames(data_t)[seq_len(p)])
  
  cv_record = matrix(NA, nrow = length(lam_set), ncol = length(nu_set))
  for(lam_idx in 1:length(lam_set)) {
    lambda = lam_set[lam_idx]
    for(nu_idx in 1:length(nu_set)) {
      nu = nu_set[nu_idx]
      # cross validation
      split_set <- sample(rep(1:folds_num, length.out = num_t))
      cindex = rep(NA,folds_num)
      for(k in 1:folds_num){
        data_t1<-data_t[split_set==k,] # validation
        data_t2<-data_t[split_set!=k,] # train
        data_ts<-rbind(data_s,data_t2)
        
        # obtain weights
        theta_ini<-rep(0,p+1) # set initial value
        Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,1:p])) # add intercept term
        Xt_new<-as.matrix(cbind(rep(1,nrow(data_t2)),data_t2[,1:p])) # add intercept term
        theta = rep(NA, p+1)
        tryCatch({
          w_opt <- nlm(density_opt(Xs_new,Xt_new,lambda), theta_ini, iterlim = 500)
          theta<-w_opt$estimate
        }, error = function(e) {
          message("lambda = ", lambda, ": error in nlm")
        })
        weights<-exp(Xs_new%*%matrix(theta,ncol=1))
        
        # fit weighted Cox regression
        tryCatch({
          cox_train<- rms::cph(cox_formula,
                               data = data_ts,weights=c(nu*weights,rep(1,dim(data_t2)[1])),
                               method = "breslow", x=TRUE,y=TRUE,surv=TRUE)
          cindex[k]<-1-Hmisc::rcorr.cens(exp(as.matrix(data_t1[,1:p])%*%stats::coef(cox_train)),
                                         survival::Surv(data_t1$time,data_t1$status))['C Index']
        }, error = function(e) {
          message("lambda = ", lambda, ": results in extreme weights and failed Cox regression")
        })
      }
      cv_record[lam_idx, nu_idx] = mean(cindex)
    }
  }
  
  ind = arrayInd(which.max(cv_record), .dim=dim(cv_record))
  lam_best<-lam_set[ind[1]]
  nu_best = nu_set[ind[2]]
  message("optimal lambda: ", lam_best)
  message("optimal nu: ", nu_best)
  
  # obtain weights
  if(!is.na(lam_best)) {
    theta_ini<-rep(0,p+1) # set initial value
    Xs_new<-as.matrix(cbind(rep(1,num_s),data_s[,1:p])) # add intercept term
    Xt_new<-as.matrix(cbind(rep(1,num_t),data_t[,1:p])) # add intercept term
    w_opt <- nlm(density_opt(Xs_new,Xt_new,lam_best), theta_ini, iterlim = 500)
    theta<-w_opt$estimate
    weights<-exp(Xs_new%*%matrix(theta,ncol=1))
  } else {
    weights = rep(1.0, num_s)
  }
  
  if(is.na(nu_best)) {
    # ratio in sample size
    nu_best = num_t/num_s
  }

  return(list(weights = weights,
              nu = nu_best))
}

run_CoxTL_ms<-function(data_t,data_s_list,p,folds_num=5,
                    lam_set = 10^(seq(-4,1,1)),
                    nu_set=10^(seq(-4,1,1))){
  nsc = length(data_s_list)
  source_sizes = vapply(data_s_list, nrow, integer(1))
  all_sc_res = lapply(seq_len(nsc), run_CoxTL_ss,
                      data_t = data_t, data_s_list = data_s_list, p = p,
                      folds_num = folds_num, lam_set = lam_set, nu_set = nu_set)
  all_weights = lapply(all_sc_res, function(x) x$weights)
  all_nu = lapply(all_sc_res, function(x) x$nu)
  
  data_t = as.data.frame(data_t)
  strata_levels = 0:nsc
  data_t$R = factor(0, levels = strata_levels)
  stratified_sources = lapply(seq_len(nsc), function(sc) {
    source_data = as.data.frame(data_s_list[[sc]])
    source_data$R = factor(sc, levels = strata_levels)
    source_data
  })
  data_s = as.data.frame(do.call(rbind, stratified_sources))
  num_t<-dim(data_t)[1]
  num_s<-dim(data_s)[1]
  strat <- rms::strat
  cox_formula = .CoxTL_formula(colnames(data_t)[seq_len(p)])
  
  for(sc in 1:nsc) {
    # tune nu for source sc while fixing nu at initial values for all other sources 
    cv_record = rep(NA, length(nu_set))
    for(nu_idx in 1:length(nu_set)) {
      nu_sc = nu_set[nu_idx]
      all_nu[[sc]] = nu_sc
      
      # cross validation
      split_set <- sample(rep(1:folds_num, length.out = num_t))
      cindex = rep(NA,folds_num)
      for(k in 1:folds_num){
        data_t1<-data_t[split_set==k,] # validation
        data_t2<-data_t[split_set!=k,] # train
        data_ts<-rbind(data_s,data_t2)
        cox_weights = unlist(lapply(1:nsc, function(x) all_nu[[x]]*all_weights[[x]]))
        
        # remove subjects with 0 weights
        rm_ind = which(cox_weights==0.0)
        if(length(rm_ind)) {
          data_ts = data_ts[-rm_ind,]
          cox_weights = cox_weights[-rm_ind]
        }
        
        # fit weighted Cox regression
        tryCatch({
          cox_train<- rms::cph(cox_formula,
                               data = data_ts,weights=c(cox_weights,rep(1,dim(data_t2)[1])),
                               method = "breslow", x=TRUE,y=TRUE,surv=TRUE)
          cindex[k]<-1-Hmisc::rcorr.cens(exp(as.matrix(data_t1[,1:p])%*%stats::coef(cox_train)),
                                         survival::Surv(data_t1$time,data_t1$status))['C Index']
        }, error = function(e) {
          message("nu = ", nu_sc, ": weighted Cox regression failed")
        })
      }
      cv_record[nu_idx] = mean(cindex, na.rm = T)
    }
    ind = which.max(cv_record)
    if(length(ind)) {
      nu_sc = nu_set[ind]
      all_nu[[sc]] = nu_sc
    } else {
      # ratio in sample size
      all_nu[[sc]] = num_t/source_sizes[sc]
    }
  }
  
  message("best nu: ")
  print(unlist(all_nu))
  cox_weights = unlist(lapply(1:nsc, function(x) all_nu[[x]]*all_weights[[x]]))
 
  # remove subjects with 0 weights
  rm_ind = which(cox_weights==0.0)
  if(length(rm_ind)) {
    data_s = data_s[-rm_ind,]
    cox_weights = cox_weights[-rm_ind]
  }
  cox_final<- rms::cph(cox_formula,
                       data=rbind(data_s,data_t),
                       weights=c(cox_weights,rep(1,num_t)),
                       method = "breslow",x=TRUE,y=TRUE,surv=TRUE)
  return(cox_final)
}
