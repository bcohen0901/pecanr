#' Partial Eta-Squared for Linear Mixed Models
#'
#' Calculates partial eta-squared effect sizes for fixed effects in linear mixed
#' models with crossed or nested random effects.
#'
#' @param model A fitted model object from \code{lme4::lmer()}.
#' @param effect Character string specifying the fixed effect to analyze.
#'   Must match a name in \code{fixef(model)}.
#' @param data Data frame used to fit the model.
#' @param design Character string: either \code{"crossed"} or \code{"nested"}.
#' @param subj_var Character string specifying the subject/participant grouping
#'   variable. Retained for backward compatibility; prefer \code{cross_vars}.
#' @param item_var Character string specifying the item/stimulus grouping
#'   variable. Retained for backward compatibility; prefer \code{cross_vars}.
#' @param cross_vars Character vector of ALL crossed grouping variable names
#'   (e.g., \code{c("subject", "item", "rater")}). Supersedes \code{subj_var}
#'   and \code{item_var} when provided. Supports any number of crossed factors.
#' @param nest_vars Character vector specifying nesting variables from lowest to
#'   highest level (required for nested designs).
#' @param effect_level Character string specifying the level at which the effect
#'   varies (e.g., "L1", "L2"). If NULL, will be detected automatically.
#' @param var_x Optional numeric. Pre-computed variance of the predictor (or
#'   interaction product). If supplied, overrides the internal \code{var()}
#'   calculation from \code{data}. Useful when raw data are unavailable but the
#'   predictor variance is known (e.g., from design: a +/-1 binary predictor has
#'   \code{var_x = 1}).
#' @param operative Logical. If TRUE, calculates operative effect size (excludes
#'   variance components that don't contribute to SE). Default is FALSE.
#' @param verbose Logical. If TRUE, prints detailed results. Default is TRUE.
#'
#' @return An object of class \code{"eta2p_lmm"} containing:
#' \item{eta2p}{Partial eta-squared value}
#' \item{variance_effect}{Variance explained by the effect}
#' \item{variance_error}{Error variance}
#' \item{effect}{Name of the effect}
#' \item{design}{Design type ("crossed" or "nested")}
#' \item{operative}{Whether operative effect size was calculated}
#' \item{variance_components}{List of individual variance components}
#' \item{...}{Additional design-specific information}
#'
#' @details
#' The function implements a variance decomposition approach for computing partial
#' eta-squared in mixed models. Random slope variances are translated to the
#' outcome scale using the formula:
#' \deqn{\sigma^2_{slope}(Y) = \sigma^2_b \times \sigma^2_X}
#'
#' For interaction effects, the variance of the predictor is calculated as
#' the variance of the actual product term (e.g., var(X1 * X2)), not the
#' product of individual variances. This correctly accounts for centering,
#' scaling, and correlation between predictors. The \code{var_x} argument
#' allows bypassing this computation when the variance is known a priori.
#'
#' For \strong{general effect sizes} (default), all variance components are
#' included in the denominator. For \strong{operative effect sizes}
#' (\code{operative = TRUE}), only variance components that contribute to the
#' standard error of the effect are included.
#'
#' Crossed designs support any number of grouping factors via \code{cross_vars}.
#' The two-argument form (\code{subj_var} + \code{item_var}) is retained for
#' backward compatibility and is equivalent to
#' \code{cross_vars = c(subj_var, item_var)}.
#'
#' @references
#' Correll, J., Mellinger, C., McClelland, G. H., & Judd, C. M. (2020).
#' Avoid Cohen's 'Small', 'Medium', and 'Large' for Power Analysis.
#' \emph{Trends in Cognitive Sciences}, 24(3), 200-207.
#'
#' Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
#' multilevel models: An integrative framework for defining R-squared measures.
#' \emph{Psychological Methods}, 24(3), 309-338.
#'
#' @examples
#' \dontrun{
#' library(lme4)
#'
#' # Two crossed factors (backward-compatible call)
#' model <- lmer(y ~ condition + (1 | subject) + (1 | item),
#'               data = crossed_data)
#' eta2p(model, "condition", crossed_data,
#'       design = "crossed",
#'       subj_var = "subject",
#'       item_var = "item")
#'
#' # Three crossed factors using cross_vars
#' model3 <- lmer(y ~ condition + (1 | subject) + (1 | item) + (1 | rater),
#'               data = three_way_data)
#' eta2p(model3, "condition", three_way_data,
#'       design = "crossed",
#'       cross_vars = c("subject", "item", "rater"))
#'
#' # Supply predictor variance directly (no raw data needed for this step)
#' eta2p(model, "condition", crossed_data,
#'       design = "crossed",
#'       cross_vars = c("subject", "item"),
#'       var_x = 1)   # +/-1 binary predictor
#' }
#'
#' @export
eta2p <- function(model, effect, data,
                  design     = c("crossed", "nested"),
                  subj_var   = NULL,
                  item_var   = NULL,
                  cross_vars = NULL,
                  nest_vars  = NULL,
                  effect_level = NULL,
                  var_x      = NULL,
                  operative  = FALSE,
                  verbose    = TRUE) {

  design <- match.arg(design)


  # INPUT VALIDATION


  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required but not installed.")
  }

  if (!inherits(model, "lmerMod")) {
    stop("'model' must be a fitted lmer model (class 'lmerMod')")
  }

  if (model@optinfo$conv$opt != 0) {
    warning("Model did not converge. Results may be unreliable.")
  }

  if (lme4::isSingular(model)) {
    warning("Model has singular fit (some random effects are zero or perfectly ",
            "correlated). Results may be unreliable.")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.null(var_x) && (!is.numeric(var_x) || length(var_x) != 1 || var_x < 0)) {
    stop("'var_x' must be a single non-negative numeric value.")
  }


  # RESOLVE CROSSED GROUPING VARIABLES
  # cross_vars supersedes subj_var / item_var; fall back for backward compat.


  if (design == "crossed") {

    if (!is.null(cross_vars)) {
      if (length(cross_vars) < 1) {
        stop("'cross_vars' must name at least one grouping variable.")
      }
      missing_cv <- cross_vars[!cross_vars %in% colnames(data)]
      if (length(missing_cv) > 0) {
        stop("cross_vars variable(s) not found in data: ",
             paste(missing_cv, collapse = ", "))
      }
    } else {
      if (is.null(subj_var)) {
        stop("For crossed designs, supply either 'cross_vars' or 'subj_var'.")
      }
      if (!subj_var %in% colnames(data)) {
        stop("Subject variable '", subj_var, "' not found in data")
      }
      if (!is.null(item_var) && !item_var %in% colnames(data)) {
        stop("Item variable '", item_var, "' not found in data")
      }
      cross_vars <- c(subj_var, if (!is.null(item_var)) item_var)
    }

    subj_var <- cross_vars[1]
    item_var  <- if (length(cross_vars) >= 2) cross_vars[2] else NULL

  }

  if (design == "nested" && is.null(nest_vars)) {
    stop("For nested designs, 'nest_vars' is required.")
  }

  if (design == "nested") {
    missing_vars <- nest_vars[!nest_vars %in% colnames(data)]
    if (length(missing_vars) > 0) {
      stop("Nesting variable(s) not found in data: ",
           paste(missing_vars, collapse = ", "))
    }
  }


  # GET MODEL COMPONENTS


  vc             <- as.data.frame(lme4::VarCorr(model))
  fixed_effects  <- lme4::fixef(model)

  if (!effect %in% names(fixed_effects)) {
    stop("Effect '", effect, "' not found in model fixed effects.\n",
         "Available effects: ", paste(names(fixed_effects), collapse = ", "))
  }


  # CALCULATE VARIANCE EXPLAINED BY THE EFFECT


  B           <- fixed_effects[effect]
  effect_vars <- strsplit(effect, ":")[[1]]

  if (!is.null(var_x)) {
    sigma2_X <- var_x

  } else if (length(effect_vars) == 1) {
    if (effect_vars == "(Intercept)") {
      sigma2_X <- 1
    } else {
      if (!effect_vars %in% colnames(data)) {
        stop("Variable '", effect_vars, "' not found in data. ",
             "Supply 'var_x' directly if raw data are unavailable.")
      }
      sigma2_X <- var(data[[effect_vars]], na.rm = TRUE)
    }

  } else {
    missing_vars <- effect_vars[!effect_vars %in% colnames(data)]
    if (length(missing_vars) > 0) {
      stop("Variable(s) not found in data: ",
           paste(missing_vars, collapse = ", "),
           ". Supply 'var_x' directly if raw data are unavailable.")
    }
    interaction_product <- Reduce("*", lapply(effect_vars, function(v) data[[v]]))
    sigma2_X <- var(interaction_product, na.rm = TRUE)
  }

  variance_effect <- B^2 * sigma2_X


  # DETECT WITHIN/BETWEEN STRUCTURE (operative only)


  within_between <- NULL
  if (operative && design == "crossed") {
    within_between <- detect_within_between(data, effect_vars,
                                            subj_var, item_var)
  }


  # CALCULATE ERROR VARIANCE


  if (design == "crossed") {
    result <- calc_error_crossed(vc, data,
                                 cross_vars     = cross_vars,
                                 operative      = operative,
                                 within_between = within_between,
                                 effect_vars    = effect_vars)
  } else {
    if (is.null(effect_level)) {
      effect_level <- detect_effect_level(data, effect, effect_vars, nest_vars)
    }
    result <- calc_error_nested(vc, data, nest_vars, effect_level,
                                operative   = operative,
                                effect_vars = effect_vars)
  }

  variance_error <- result$variance_error


  # CALCULATE PARTIAL ETA-SQUARED


  eta2p <- variance_effect / (variance_effect + variance_error)


  # BUILD OUTPUT


  if (design == "crossed") {
    n_per_factor <- sapply(cross_vars, function(v) length(unique(data[[v]])))
    names(n_per_factor) <- cross_vars

    output <- list(
      eta2p              = eta2p,
      variance_effect    = variance_effect,
      variance_error     = variance_error,
      effect             = effect,
      design             = design,
      operative          = operative,
      within_between     = within_between,
      cross_vars         = cross_vars,
      n_per_factor       = n_per_factor,
      n_subj             = n_per_factor[subj_var],
      n_item             = if (!is.null(item_var)) n_per_factor[item_var] else NA,
      variance_components = result$components
    )

  } else {
    n_levels <- sapply(nest_vars, function(v) length(unique(data[[v]])))
    names(n_levels) <- nest_vars

    output <- list(
      eta2p              = eta2p,
      variance_effect    = variance_effect,
      variance_error     = variance_error,
      effect             = effect,
      design             = design,
      operative          = operative,
      effect_level       = effect_level,
      n_levels           = n_levels,
      variance_components = result$components
    )
  }

  class(output) <- c("eta2p_lmm", "list")

  if (verbose) print(output)

  return(output)
}


#' Batch Calculate Partial Eta-Squared for Multiple Effects
#'
#' Calculates partial eta-squared for all fixed effects in a model.
#'
#' @inheritParams eta2p
#'
#' @return A data frame with one row per effect containing eta-squared values
#'   and variance components.
#'
#' @export
batch_eta2p <- function(model, data,
                        design     = c("crossed", "nested"),
                        subj_var   = NULL,
                        item_var   = NULL,
                        cross_vars = NULL,
                        nest_vars  = NULL,
                        operative  = FALSE,
                        verbose    = FALSE) {

  design <- match.arg(design)

  effects <- names(lme4::fixef(model))
  effects <- effects[effects != "(Intercept)"]

  if (length(effects) == 0) stop("No fixed effects found (excluding intercept)")

  results <- lapply(effects, function(eff) {
    tryCatch({
      res <- eta2p(model, eff, data,
                   design     = design,
                   subj_var   = subj_var,
                   item_var   = item_var,
                   cross_vars = cross_vars,
                   nest_vars  = nest_vars,
                   operative  = operative,
                   verbose    = FALSE)

      if (design == "nested") {
        data.frame(effect = eff, eta2p = res$eta2p,
                   effect_level = res$effect_level,
                   variance_effect = res$variance_effect,
                   variance_error  = res$variance_error,
                   type = if (operative) "operative" else "general",
                   stringsAsFactors = FALSE)
      } else {
        df <- data.frame(effect = eff, eta2p = res$eta2p,
                         variance_effect = res$variance_effect,
                         variance_error  = res$variance_error,
                         type = if (operative) "operative" else "general",
                         stringsAsFactors = FALSE)
        if (operative && !is.null(res$within_between)) {
          df$within_subj <- res$within_between$subj
          df$within_item <- res$within_between$item
        }
        df
      }
    }, error = function(e) {
      warning("Error calculating eta2p for '", eff, "': ", e$message, call. = FALSE)
      data.frame(effect = eff, eta2p = NA,
                 effect_level    = if (design == "nested") NA else NULL,
                 variance_effect = NA,
                 variance_error  = NA,
                 type = if (operative) "operative" else "general",
                 stringsAsFactors = FALSE)
    })
  })

  results_df <- do.call(rbind, results)
  results_df <- results_df[order(results_df$eta2p, decreasing = TRUE, na.last = TRUE), ]
  rownames(results_df) <- NULL

  if (verbose) print(results_df)

  return(results_df)
}



# INTERNAL HELPER FUNCTIONS (not exported)


detect_within_between <- function(data, effect_vars, subj_var, item_var) {

  check_var <- effect_vars[1]

  if (!check_var %in% colnames(data)) {
    warning("Cannot detect within/between structure: '", check_var,
            "' not found in data.")
    return(list(subj = "unknown", item = "unknown"))
  }

  subj_within <- check_within_factor(data, check_var, subj_var)
  item_within <- if (!is.null(item_var)) {
    check_within_factor(data, check_var, item_var)
  } else {
    FALSE
  }

  list(
    subj = if (subj_within) "within" else "between",
    item = if (item_within) "within" else if (!is.null(item_var)) "between" else "none"
  )
}

check_within_factor <- function(data, effect_var, group_var) {
  var_within <- tapply(data[[effect_var]], data[[group_var]], var, na.rm = TRUE)
  any(var_within > 1e-10, na.rm = TRUE)
}

slope_contribution <- function(slope_var_name, vc_row_vcov, data) {

  if (grepl(":", slope_var_name)) {
    parts <- strsplit(slope_var_name, ":")[[1]]
    if (!all(parts %in% colnames(data))) {
      warning("Slope variable '", slope_var_name,
              "' components not found in data - skipping.")
      return(0)
    }
    product <- Reduce("*", lapply(parts, function(v) data[[v]]))
    sigma2_X <- var(product, na.rm = TRUE)
  } else {
    if (!slope_var_name %in% colnames(data)) {
      warning("Slope variable '", slope_var_name, "' not found in data - skipping.")
      return(0)
    }
    sigma2_X <- var(data[[slope_var_name]], na.rm = TRUE)
  }

  vc_row_vcov * sigma2_X
}

slope_matches_effect <- function(slope_var_name, effect_vars) {
  slope_parts <- if (grepl(":", slope_var_name)) {
    strsplit(slope_var_name, ":")[[1]]
  } else {
    slope_var_name
  }
  all(slope_parts %in% effect_vars) && all(effect_vars %in% slope_parts)
}

calc_error_crossed <- function(vc, data, cross_vars,
                               operative      = FALSE,
                               within_between = NULL,
                               effect_vars    = NULL) {

  residual_var <- vc$vcov[vc$grp == "Residual"]
  residual_var <- if (length(residual_var) == 0) 0 else residual_var[1]

  factor_intercept_vars <- setNames(numeric(length(cross_vars)), cross_vars)
  factor_slope_vars     <- setNames(vector("list", length(cross_vars)), cross_vars)

  for (cv in cross_vars) {

    int_row <- vc[vc$grp == cv & vc$var1 == "(Intercept)" & is.na(vc$var2), ]
    factor_intercept_vars[cv] <- if (nrow(int_row) == 0) 0 else int_row$vcov[1]

    slp_rows <- vc[vc$grp == cv &
                     !is.na(vc$var1) &
                     vc$var1 != "(Intercept)" &
                     is.na(vc$var2), ]

    contributions <- setNames(numeric(nrow(slp_rows)),
                              if (nrow(slp_rows) > 0) slp_rows$var1 else character(0))

    for (i in seq_len(nrow(slp_rows))) {
      nm   <- slp_rows$var1[i]
      val  <- slope_contribution(nm, slp_rows$vcov[i], data)
      contributions[nm] <- val
    }

    factor_slope_vars[[cv]] <- contributions
  }

  if (operative && !is.null(within_between) && !is.null(effect_vars)) {

    operative_slope_total <- 0
    for (cv in cross_vars) {
      for (nm in names(factor_slope_vars[[cv]])) {
        if (slope_matches_effect(nm, effect_vars)) {
          operative_slope_total <- operative_slope_total + factor_slope_vars[[cv]][nm]
        }
      }
    }

    variance_error <- residual_var + operative_slope_total

    if (within_between$subj == "between") {
      variance_error <- variance_error + factor_intercept_vars[cross_vars[1]]
    }
    if (length(cross_vars) >= 2 && within_between$item == "between") {
      variance_error <- variance_error + factor_intercept_vars[cross_vars[2]]
    }
    if (length(cross_vars) > 2) {
      for (cv in cross_vars[-(1:2)]) {
        variance_error <- variance_error + factor_intercept_vars[cv]
      }
    }

  } else {
    variance_error <- residual_var +
      sum(factor_intercept_vars) +
      sum(unlist(factor_slope_vars))
  }

  components <- list(
    residual          = residual_var,
    factor_intercepts = factor_intercept_vars,
    factor_slopes     = factor_slope_vars,
    slopes_total      = sum(unlist(factor_slope_vars)),
    subj_intercept    = factor_intercept_vars[cross_vars[1]],
    item_intercept    = if (length(cross_vars) >= 2) factor_intercept_vars[cross_vars[2]] else 0
  )

  list(variance_error = variance_error, components = components)
}

calc_error_nested <- function(vc, data, nest_vars, effect_level,
                              operative = FALSE, effect_vars = NULL) {

  n_levels  <- length(nest_vars)
  level_num <- as.integer(gsub("L", "", effect_level))

  residual_var <- vc$vcov[vc$grp == "Residual"]
  residual_var <- if (length(residual_var) == 0) 0 else residual_var[1]

  level_intercept_vars <- sapply(nest_vars, function(v) {
    row <- vc[vc$grp == v & vc$var1 == "(Intercept)" & is.na(vc$var2), ]
    if (nrow(row) == 0) {
      nested_grp <- grep(paste0("^", v, ":"), vc$grp, value = TRUE)
      if (length(nested_grp) > 0) {
        row <- vc[vc$grp == nested_grp[1] & vc$var1 == "(Intercept)" & is.na(vc$var2), ]
      }
    }
    if (nrow(row) == 0) 0 else row$vcov[1]
  })
  names(level_intercept_vars) <- nest_vars

  level_slope_vars <- sapply(nest_vars, function(v) {
    slp_rows <- vc[vc$grp == v & !is.na(vc$var1) &
                     vc$var1 != "(Intercept)" & is.na(vc$var2), ]
    if (nrow(slp_rows) == 0) {
      nested_grp <- grep(paste0("^", v, ":"), vc$grp, value = TRUE)
      if (length(nested_grp) > 0) {
        slp_rows <- vc[vc$grp == nested_grp[1] & !is.na(vc$var1) &
                         vc$var1 != "(Intercept)" & is.na(vc$var2), ]
      }
    }
    if (nrow(slp_rows) == 0) return(0)

    total <- 0
    for (i in seq_len(nrow(slp_rows))) {
      nm <- slp_rows$var1[i]
      if (operative && !is.null(effect_vars) && !slope_matches_effect(nm, effect_vars)) {
        next
      }
      total <- total + slope_contribution(nm, slp_rows$vcov[i], data)
    }
    total
  })
  names(level_slope_vars) <- nest_vars

  if (level_num == 1) {
    variance_error <- residual_var +
      sum(level_intercept_vars) +
      sum(level_slope_vars)
  } else {
    idx <- level_num:n_levels
    variance_error <- sum(level_intercept_vars[idx]) +
      sum(level_slope_vars[idx])
  }

  components <- list(
    residual         = residual_var,
    level_intercepts = level_intercept_vars,
    level_slopes     = level_slope_vars,
    effect_level     = effect_level,
    levels_included  = if (level_num == 1) c("residual", nest_vars) else nest_vars[level_num:n_levels]
  )

  list(variance_error = variance_error, components = components)
}

detect_effect_level <- function(data, effect, effect_vars, nest_vars) {

  check_var <- effect_vars[1]
  n_levels  <- length(nest_vars)

  if (!check_var %in% colnames(data)) {
    warning("Cannot detect effect level: '", check_var,
            "' not found in data. Defaulting to L1.")
    return("L1")
  }

  for (i in seq_along(nest_vars)) {
    grp_var    <- nest_vars[i]
    counts     <- table(data[[grp_var]])
    if (all(counts == 1)) next
    var_within <- tapply(data[[check_var]], data[[grp_var]], var, na.rm = TRUE)
    if (any(var_within > 1e-10, na.rm = TRUE)) {
      return(paste0("L", max(1, i - 1)))
    }
  }

  paste0("L", n_levels)
}
