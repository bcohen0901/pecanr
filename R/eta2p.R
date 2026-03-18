#' Partial Eta-Squared for Linear Mixed Models
#'
#' Calculates partial eta-squared effect sizes for fixed effects in linear mixed
#' models with crossed, nested, or mixed (crossed-and-nested) random effects.
#'
#' @param model A fitted model object from \code{lme4::lmer()}.
#' @param effect Character string specifying the fixed effect to analyze.
#'   Must match a name in \code{fixef(model)}.
#' @param data Data frame used to fit the model.
#' @param design Character string: \code{"crossed"}, \code{"nested"}, or
#'   \code{"mixed"}. Use \code{"mixed"} when the model contains both crossed and
#'   nested random effects (e.g., photos nested within models, but both crossed
#'   with participants).
#' @param subj_var Character string specifying the subject/participant grouping
#'   variable. Retained for backward compatibility; prefer \code{cross_vars}.
#' @param item_var Character string specifying the item/stimulus grouping
#'   variable. Retained for backward compatibility; prefer \code{cross_vars}.
#' @param cross_vars Character vector of grouping variable names that are
#'   \emph{crossed} with each other and with the nested hierarchy
#'   (e.g., \code{c("participant")} or \code{c("participant", "rater")}).
#'   Supersedes \code{subj_var} and \code{item_var} when provided.
#' @param nest_vars Character vector of nesting variables ordered from lowest
#'   to highest level (e.g., \code{c("photo", "model")} when photos are nested
#'   within models). Required for \code{design = "nested"} or
#'   \code{design = "mixed"}.
#' @param effect_level Character string specifying the level at which the effect
#'   varies (e.g., \code{"L1"}, \code{"L2"}). If \code{NULL}, detected
#'   automatically from \code{nest_vars}.
#' @param var_x Optional numeric. Pre-computed variance of the predictor (or
#'   interaction product). If supplied, overrides the internal \code{var()}
#'   calculation from \code{data}. Useful when raw data are unavailable but the
#'   predictor variance is known (e.g., from design: a +/-1 binary predictor has
#'   \code{var_x = 1}).
#' @param operative Logical. If \code{TRUE}, calculates operative effect size
#'   (excludes variance components that do not contribute to the SE of the
#'   effect). Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, prints detailed results. Default is
#'   \code{TRUE}.
#'
#' @return An object of class \code{"eta2p_lmm"} containing:
#' \item{eta2p}{Partial eta-squared value.}
#' \item{variance_effect}{Variance explained by the effect.}
#' \item{variance_error}{Error variance (denominator).}
#' \item{effect}{Name of the effect.}
#' \item{design}{Design type: \code{"crossed"}, \code{"nested"}, or
#'   \code{"mixed"}.}
#' \item{operative}{Whether operative effect size was calculated.}
#' \item{variance_components}{List of individual variance components. For
#'   \code{"mixed"} designs this is a list with \code{$crossed} and
#'   \code{$nested} sub-lists.}
#' \item{within_between}{(Mixed/crossed designs) Named list or vector
#'   classifying each grouping factor as \code{"within"} or \code{"between"}
#'   with respect to the focal effect. Only populated when
#'   \code{operative = TRUE}.}
#' \item{cross_vars}{(Crossed/mixed designs) Crossed grouping variable names.}
#' \item{nest_vars}{(Nested/mixed designs) Nested grouping variable names.}
#' \item{n_per_factor}{(Crossed designs) Number of units per crossed factor.}
#' \item{n_cross}{(Mixed designs) Number of units per crossed factor.}
#' \item{n_nested}{(Mixed designs) Number of units per nested factor.}
#' \item{n_levels}{(Nested designs) Number of units per nested level.}
#' \item{effect_level}{(Nested/mixed designs) Detected or supplied effect
#'   level.}
#'
#' @details
#' \subsection{Variance decomposition}{
#' The function implements a variance decomposition approach for computing
#' partial eta-squared in mixed models. Random slope variances are translated
#' to the outcome scale using:
#' \deqn{\sigma^2_{\text{slope}}(Y) = \sigma^2_b \times \sigma^2_X}
#' For interaction effects, \eqn{\sigma^2_X} is computed as the variance of
#' the actual product term, correctly accounting for centering, scaling, and
#' predictor correlation. The \code{var_x} argument allows bypassing this
#' computation when the variance is known a priori.
#' }
#'
#' \subsection{General vs. operative effect sizes}{
#' For \strong{general effect sizes} (default), all variance components appear
#' in the denominator. For \strong{operative effect sizes}
#' (\code{operative = TRUE}), only components that contribute to the standard
#' error of the effect are included. In mixed designs, each grouping factor
#' (both crossed and nested) is independently classified as
#' \code{"within"} or \code{"between"} the focal effect; between-subjects
#' intercept variances are excluded from the operative denominator.
#' }
#'
#' \subsection{Mixed designs}{
#' Use \code{design = "mixed"} when the model contains both crossed and nested
#' random effects simultaneously. A canonical example is participants viewing
#' multiple photos of each model (stimulus): photos are nested within models,
#' but both levels are crossed with participants.
#'
#' The corresponding \pkg{lme4} model and \code{eta2p} call would be:
#' \preformatted{
#' fit <- lmer(y ~ x + (1 | participant) + (1 | model) + (1 | photo:model),
#'             data = dat)
#' eta2p(fit, "x", dat,
#'       design     = "mixed",
#'       cross_vars = "participant",
#'       nest_vars  = c("photo", "model"))
#' }
#' Residual variance is attributed to the crossed side and counted exactly
#' once; the nested calculator's residual is subtracted to prevent
#' double-counting.
#' }
#'
#' \subsection{Crossed designs}{
#' Supports any number of grouping factors via \code{cross_vars}. The
#' two-argument form (\code{subj_var} + \code{item_var}) is retained for
#' backward compatibility and is equivalent to
#' \code{cross_vars = c(subj_var, item_var)}.
#' }
#'
#' @references
#' Correll, J., Mellinger, C., McClelland, G. H., & Judd, C. M. (2020).
#' Avoid Cohen's 'Small', 'Medium', and 'Large' for Power Analysis.
#' \emph{Trends in Cognitive Sciences}, 24(3), 200--207.
#' \doi{10.1016/j.tics.2019.12.009}
#'
#' Correll, J., Mellinger, C., & Pedersen, E. J. (2022). Flexible approaches
#' for estimating partial eta squared in mixed-effects models with crossed
#' random factors. \emph{Behavior Research Methods}, 54, 1626--1642.
#' \doi{10.3758/s13428-021-01687-2}
#'
#' Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
#' multilevel models: An integrative framework for defining R-squared measures.
#' \emph{Psychological Methods}, 24(3), 309--338.
#' \doi{10.1037/met0000184}
#'
#' @examples
#' \donttest{
#' library(lme4)
#'
#' # --- Two crossed factors (backward-compatible call) ---
#' set.seed(42)
#' crossed_data <- data.frame(
#'   y         = rnorm(120),
#'   condition = rep(c(-0.5, 0.5), 60),
#'   subject   = factor(rep(1:20, each = 6)),
#'   item      = factor(rep(1:6, 20))
#' )
#' model_c <- lmer(y ~ condition + (1 | subject) + (1 | item),
#'                 data = crossed_data)
#' eta2p(model_c, "condition", crossed_data,
#'       design   = "crossed",
#'       subj_var = "subject",
#'       item_var = "item")
#'
#' # --- Three crossed factors using cross_vars ---
#' set.seed(42)
#' three_way_data <- data.frame(
#'   y         = rnorm(180),
#'   condition = rep(c(-0.5, 0.5), 90),
#'   subject   = factor(rep(1:20, each = 9)),
#'   item      = factor(rep(rep(1:6, each = 3), 10)),
#'   rater     = factor(rep(1:3, 60))
#' )
#' model_3 <- lmer(y ~ condition +
#'                   (1 | subject) + (1 | item) + (1 | rater),
#'                 data = three_way_data)
#' eta2p(model_3, "condition", three_way_data,
#'       design     = "crossed",
#'       cross_vars = c("subject", "item", "rater"))
#'
#' # --- Mixed design: photos nested within models, crossed with participants ---
#' # Participants each view multiple photos of multiple models.
#' # Photos are nested within models (each photo belongs to one model),
#' # but both levels are crossed with participants.
#' set.seed(42)
#' n_subj  <- 30; n_model <- 10; n_photo <- 3
#' mixed_data <- expand.grid(
#'   participant = factor(seq_len(n_subj)),
#'   photo_id    = factor(seq_len(n_model * n_photo))
#' )
#' mixed_data$model_id <- factor(
#'   rep(seq_len(n_model), each = n_photo)[as.integer(mixed_data$photo_id)]
#' )
#' mixed_data$x <- rnorm(nrow(mixed_data))
#' mixed_data$y <- rnorm(nrow(mixed_data))
#'
#' model_m <- lmer(y ~ x +
#'                   (1 | participant) +
#'                   (1 | model_id) +
#'                   (1 | photo_id:model_id),
#'                 data = mixed_data)
#' eta2p(model_m, "x", mixed_data,
#'       design     = "mixed",
#'       cross_vars = "participant",
#'       nest_vars  = c("photo_id", "model_id"))
#'
#' # --- Supply predictor variance directly (var_x) ---
#' eta2p(model_c, "condition", crossed_data,
#'       design     = "crossed",
#'       cross_vars = c("subject", "item"),
#'       var_x      = 1)   # +/-1 binary predictor: var = 1 by design
#' }
#'
#' @export
eta2p <- function(model, effect, data,
                  design       = c("crossed", "nested", "mixed"),
                  subj_var     = NULL,
                  item_var     = NULL,
                  cross_vars   = NULL,
                  nest_vars    = NULL,
                  effect_level = NULL,
                  var_x        = NULL,
                  operative    = FALSE,
                  verbose      = TRUE) {

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

  if (!is.null(var_x) &&
      (!is.numeric(var_x) || length(var_x) != 1 || var_x < 0)) {
    stop("'var_x' must be a single non-negative numeric value.")
  }


  # RESOLVE GROUPING VARIABLES PER DESIGN

  # --- crossed ---
  if (design == "crossed") {

    if (!is.null(cross_vars)) {
      if (length(cross_vars) < 1)
        stop("'cross_vars' must name at least one grouping variable.")
      missing_cv <- cross_vars[!cross_vars %in% colnames(data)]
      if (length(missing_cv) > 0)
        stop("cross_vars variable(s) not found in data: ",
             paste(missing_cv, collapse = ", "))
    } else {
      if (is.null(subj_var))
        stop("For crossed designs, supply either 'cross_vars' or 'subj_var'.")
      if (!subj_var %in% colnames(data))
        stop("Subject variable '", subj_var, "' not found in data")
      if (!is.null(item_var) && !item_var %in% colnames(data))
        stop("Item variable '", item_var, "' not found in data")
      cross_vars <- c(subj_var, if (!is.null(item_var)) item_var)
    }

    subj_var <- cross_vars[1]
    item_var <- if (length(cross_vars) >= 2) cross_vars[2] else NULL
  }

  # --- nested ---
  if (design == "nested") {
    if (is.null(nest_vars))
      stop("For nested designs, 'nest_vars' is required.")
    missing_vars <- nest_vars[!nest_vars %in% colnames(data)]
    if (length(missing_vars) > 0)
      stop("Nesting variable(s) not found in data: ",
           paste(missing_vars, collapse = ", "))
  }

  # --- mixed ---
  if (design == "mixed") {

    # Must supply both sets of variables
    if (is.null(cross_vars) && is.null(subj_var))
      stop("For mixed designs, supply 'cross_vars' (or 'subj_var') AND 'nest_vars'.")
    if (is.null(nest_vars))
      stop("For mixed designs, 'nest_vars' is required.")

    # Resolve cross_vars (same logic as crossed branch)
    if (is.null(cross_vars)) {
      cross_vars <- c(subj_var, if (!is.null(item_var)) item_var)
    }
    missing_cv <- cross_vars[!cross_vars %in% colnames(data)]
    if (length(missing_cv) > 0)
      stop("cross_vars variable(s) not found in data: ",
           paste(missing_cv, collapse = ", "))

    # Validate nest_vars
    missing_nv <- nest_vars[!nest_vars %in% colnames(data)]
    if (length(missing_nv) > 0)
      stop("nest_vars variable(s) not found in data: ",
           paste(missing_nv, collapse = ", "))

    # Warn if any variable appears in both sets
    overlap <- intersect(cross_vars, nest_vars)
    if (length(overlap) > 0)
      warning("Variable(s) appear in both 'cross_vars' and 'nest_vars': ",
              paste(overlap, collapse = ", "),
              ". This is unusual -- please verify your design specification.")

    subj_var <- cross_vars[1]
    item_var <- if (length(cross_vars) >= 2) cross_vars[2] else NULL
  }


  # GET MODEL COMPONENTS

  vc            <- as.data.frame(lme4::VarCorr(model))
  fixed_effects <- lme4::fixef(model)

  if (!effect %in% names(fixed_effects)) {
    stop("Effect '", effect, "' not found in model fixed effects.\n",
         "Available effects: ", paste(names(fixed_effects), collapse = ", "))
  }


  # VARIANCE EXPLAINED BY THE EFFECT

  B           <- fixed_effects[effect]
  effect_vars <- strsplit(effect, ":")[[1]]

  if (!is.null(var_x)) {
    sigma2_X <- var_x

  } else if (length(effect_vars) == 1) {
    if (effect_vars == "(Intercept)") {
      sigma2_X <- 1
    } else {
      if (!effect_vars %in% colnames(data))
        stop("Variable '", effect_vars, "' not found in data. ",
             "Supply 'var_x' directly if raw data are unavailable.")
      sigma2_X <- var(data[[effect_vars]], na.rm = TRUE)
    }

  } else {
    missing_vars <- effect_vars[!effect_vars %in% colnames(data)]
    if (length(missing_vars) > 0)
      stop("Variable(s) not found in data: ",
           paste(missing_vars, collapse = ", "),
           ". Supply 'var_x' directly if raw data are unavailable.")
    interaction_product <- Reduce("*", lapply(effect_vars, function(v) data[[v]]))
    sigma2_X <- var(interaction_product, na.rm = TRUE)
  }

  variance_effect <- B^2 * sigma2_X


  # WITHIN/BETWEEN DETECTION (operative only)
  # For mixed designs, every grouping factor (crossed AND nested) is
  # classified; this drives the operative denominator correctly.

  within_between <- NULL

  if (operative) {
    if (design == "crossed") {
      within_between <- detect_within_between(data, effect_vars,
                                              subj_var, item_var)
    } else if (design == "mixed") {
      within_between <- detect_within_between_mixed(data, effect_vars,
                                                    cross_vars, nest_vars)
    }
    # nested: operative handled inside calc_error_nested via effect_level
  }


  # ERROR VARIANCE

  if (design == "crossed") {
    result <- calc_error_crossed(vc, data,
                                 cross_vars     = cross_vars,
                                 operative      = operative,
                                 within_between = within_between,
                                 effect_vars    = effect_vars)

  } else if (design == "nested") {
    if (is.null(effect_level))
      effect_level <- detect_effect_level(data, effect, effect_vars, nest_vars)
    result <- calc_error_nested(vc, data, nest_vars, effect_level,
                                operative   = operative,
                                effect_vars = effect_vars)

  } else {  # "mixed"
    if (is.null(effect_level))
      effect_level <- detect_effect_level(data, effect, effect_vars, nest_vars)
    result <- calc_error_mixed(vc, data,
                               cross_vars     = cross_vars,
                               nest_vars      = nest_vars,
                               effect_level   = effect_level,
                               operative      = operative,
                               within_between = within_between,
                               effect_vars    = effect_vars)
  }

  variance_error <- result$variance_error


  # PARTIAL ETA-SQUARED

  eta2p_val <- variance_effect / (variance_effect + variance_error)


  # BUILD OUTPUT

  if (design == "crossed") {

    n_per_factor <- sapply(cross_vars,
                           function(v) length(unique(data[[v]])))
    names(n_per_factor) <- cross_vars

    output <- list(
      eta2p               = eta2p_val,
      variance_effect     = variance_effect,
      variance_error      = variance_error,
      effect              = effect,
      design              = design,
      operative           = operative,
      within_between      = within_between,
      cross_vars          = cross_vars,
      n_per_factor        = n_per_factor,
      n_subj              = n_per_factor[subj_var],
      n_item              = if (!is.null(item_var)) n_per_factor[item_var] else NA,
      variance_components = result$components
    )

  } else if (design == "nested") {

    n_levels <- sapply(nest_vars, function(v) length(unique(data[[v]])))
    names(n_levels) <- nest_vars

    output <- list(
      eta2p               = eta2p_val,
      variance_effect     = variance_effect,
      variance_error      = variance_error,
      effect              = effect,
      design              = design,
      operative           = operative,
      effect_level        = effect_level,
      n_levels            = n_levels,
      variance_components = result$components
    )

  } else {  # "mixed"

    n_cross  <- sapply(cross_vars, function(v) length(unique(data[[v]])))
    n_nested <- sapply(nest_vars,  function(v) length(unique(data[[v]])))
    names(n_cross)  <- cross_vars
    names(n_nested) <- nest_vars

    output <- list(
      eta2p               = eta2p_val,
      variance_effect     = variance_effect,
      variance_error      = variance_error,
      effect              = effect,
      design              = design,
      operative           = operative,
      within_between      = within_between,
      cross_vars          = cross_vars,
      nest_vars           = nest_vars,
      effect_level        = effect_level,
      n_cross             = n_cross,
      n_nested            = n_nested,
      variance_components = result$components
    )
  }

  class(output) <- c("eta2p_lmm", "list")

  if (verbose) print(output)

  return(output)
}


#' Batch Calculate Partial Eta-Squared for Multiple Effects
#'
#' Calculates partial eta-squared for all fixed effects (excluding the
#' intercept) in a model.
#'
#' @inheritParams eta2p
#'
#' @return A data frame with one row per effect containing:
#' \item{effect}{Name of the fixed effect.}
#' \item{eta2p}{Partial eta-squared value.}
#' \item{effect_level}{Effect level (nested/mixed designs only).}
#' \item{variance_effect}{Variance attributed to the effect.}
#' \item{variance_error}{Error variance used in the denominator.}
#' \item{type}{\code{"general"} or \code{"operative"}.}
#' \item{within_subj, within_item}{(Crossed/mixed, operative only) Whether the
#'   effect is within subjects / within items.}
#'
#' Rows are sorted by \code{eta2p} in decreasing order.
#'
#' @export
batch_eta2p <- function(model, data,
                        design       = c("crossed", "nested", "mixed"),
                        subj_var     = NULL,
                        item_var     = NULL,
                        cross_vars   = NULL,
                        nest_vars    = NULL,
                        operative    = FALSE,
                        verbose      = FALSE) {

  design <- match.arg(design)

  effects <- names(lme4::fixef(model))
  effects <- effects[effects != "(Intercept)"]

  if (length(effects) == 0)
    stop("No fixed effects found (excluding intercept)")

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

      df <- data.frame(
        effect          = eff,
        eta2p           = res$eta2p,
        variance_effect = res$variance_effect,
        variance_error  = res$variance_error,
        type            = if (operative) "operative" else "general",
        stringsAsFactors = FALSE
      )

      if (design %in% c("nested", "mixed")) {
        df$effect_level <- res$effect_level
      }

      if (operative && !is.null(res$within_between)) {
        wb <- res$within_between
        # within_between is always a named list/vector keyed by variable name
        # (covers crossed and mixed designs; nested does not populate it)
        for (nm in names(wb)) {
          df[[paste0("within_", nm)]] <- wb[[nm]]
        }
      }

      df

    }, error = function(e) {
      warning("Error calculating eta2p for '", eff, "': ",
              e$message, call. = FALSE)
      data.frame(
        effect          = eff,
        eta2p           = NA_real_,
        variance_effect = NA_real_,
        variance_error  = NA_real_,
        type            = if (operative) "operative" else "general",
        stringsAsFactors = FALSE
      )
    })
  })

  results_df <- do.call(rbind, results)
  results_df <- results_df[
    order(results_df$eta2p, decreasing = TRUE, na.last = TRUE), ]
  rownames(results_df) <- NULL

  if (verbose) print(results_df)

  return(results_df)
}


# INTERNAL HELPERS (not exported)

# Within/between detection

#' Classify crossed grouping factors as within or between (crossed design)
#' Keys are the actual variable names so the operative branch in
#' calc_error_crossed can look them up consistently with the mixed design path.
#' @noRd
detect_within_between <- function(data, effect_vars, subj_var, item_var) {

  check_var <- effect_vars[1]

  if (!check_var %in% colnames(data)) {
    warning("Cannot detect within/between structure: '", check_var,
            "' not found in data.")
    out <- list("unknown", if (!is.null(item_var)) "unknown")
    names(out) <- c(subj_var, item_var)
    return(out)
  }

  out <- list()
  out[[subj_var]] <- if (check_within_factor(data, check_var, subj_var))
    "within" else "between"
  if (!is.null(item_var))
    out[[item_var]] <- if (check_within_factor(data, check_var, item_var))
      "within" else "between"
  out
}

#' Classify ALL grouping factors (crossed + nested) as within or between
#' (mixed design, used for operative effect sizes)
#'
#' Returns a named character vector with "within" or "between" for every
#' factor in cross_vars and nest_vars. For nested factors the classification
#' is based on whether the predictor varies within groups defined by that
#' factor.
#' @noRd
detect_within_between_mixed <- function(data, effect_vars,
                                        cross_vars, nest_vars) {

  check_var <- effect_vars[1]

  all_vars <- c(cross_vars, nest_vars)
  wb       <- character(length(all_vars))
  names(wb) <- all_vars

  if (!check_var %in% colnames(data)) {
    warning("Cannot detect within/between structure: '", check_var,
            "' not found in data. All factors set to 'unknown'.")
    wb[] <- "unknown"
    return(wb)
  }

  for (gv in all_vars) {
    if (!gv %in% colnames(data)) {
      wb[gv] <- "unknown"
      next
    }
    wb[gv] <- if (check_within_factor(data, check_var, gv)) "within" else "between"
  }

  wb
}

#' Returns TRUE when the effect variable varies within groups of group_var
#' @noRd
check_within_factor <- function(data, effect_var, group_var) {
  var_within <- tapply(data[[effect_var]], data[[group_var]], var, na.rm = TRUE)
  any(var_within > 1e-10, na.rm = TRUE)
}


# Slope helpers

#' Translate a random slope variance component to the outcome scale
#' (sigma2_b * sigma2_X).  Returns 0 and warns if variables are missing.
#' @noRd
slope_contribution <- function(slope_var_name, vc_row_vcov, data) {

  if (grepl(":", slope_var_name)) {
    parts <- strsplit(slope_var_name, ":")[[1]]
    if (!all(parts %in% colnames(data))) {
      warning("Slope variable '", slope_var_name,
              "' components not found in data - skipping.")
      return(0)
    }
    product  <- Reduce("*", lapply(parts, function(v) data[[v]]))
    sigma2_X <- var(product, na.rm = TRUE)
  } else {
    if (!slope_var_name %in% colnames(data)) {
      warning("Slope variable '", slope_var_name,
              "' not found in data - skipping.")
      return(0)
    }
    sigma2_X <- var(data[[slope_var_name]], na.rm = TRUE)
  }

  vc_row_vcov * sigma2_X
}

#' Returns TRUE when the slope variable matches the focal effect exactly
#' (same terms, any order)
#' @noRd
slope_matches_effect <- function(slope_var_name, effect_vars) {
  slope_parts <- if (grepl(":", slope_var_name)) {
    strsplit(slope_var_name, ":")[[1]]
  } else {
    slope_var_name
  }
  all(slope_parts %in% effect_vars) && all(effect_vars %in% slope_parts)
}


# Error-variance calculators

#' Compute error variance for fully crossed designs
#' @noRd
calc_error_crossed <- function(vc, data, cross_vars,
                               operative      = FALSE,
                               within_between = NULL,
                               effect_vars    = NULL) {

  residual_var <- vc$vcov[vc$grp == "Residual"]
  residual_var <- if (length(residual_var) == 0) 0 else residual_var[1]

  factor_intercept_vars <- setNames(numeric(length(cross_vars)), cross_vars)
  factor_slope_vars     <- setNames(vector("list", length(cross_vars)), cross_vars)

  for (cv in cross_vars) {
    int_row <- vc[vc$grp == cv &
                    vc$var1 == "(Intercept)" &
                    is.na(vc$var2), ]
    factor_intercept_vars[cv] <- if (nrow(int_row) == 0) 0 else int_row$vcov[1]

    slp_rows <- vc[vc$grp == cv &
                     !is.na(vc$var1) &
                     vc$var1 != "(Intercept)" &
                     is.na(vc$var2), ]

    contributions <- setNames(
      numeric(nrow(slp_rows)),
      if (nrow(slp_rows) > 0) slp_rows$var1 else character(0)
    )

    for (i in seq_len(nrow(slp_rows))) {
      nm                <- slp_rows$var1[i]
      contributions[nm] <- slope_contribution(nm, slp_rows$vcov[i], data)
    }

    factor_slope_vars[[cv]] <- contributions
  }

  if (operative && !is.null(within_between) && !is.null(effect_vars)) {

    # For operative: only matching slopes, and between-factor intercepts
    operative_slope_total <- 0
    for (cv in cross_vars) {
      for (nm in names(factor_slope_vars[[cv]])) {
        if (slope_matches_effect(nm, effect_vars)) {
          operative_slope_total <- operative_slope_total +
            factor_slope_vars[[cv]][nm]
        }
      }
    }

    variance_error <- residual_var + operative_slope_total

    # Add intercept variance for any factor where the effect is "between"
    for (cv in cross_vars) {
      status <- within_between[[cv]]   # works for both list and named vector
      if (!is.null(status) && identical(status, "between")) {
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
    item_intercept    = if (length(cross_vars) >= 2)
      factor_intercept_vars[cross_vars[2]] else 0
  )

  list(variance_error = variance_error, components = components)
}


#' Compute error variance for fully nested designs
#' @noRd
calc_error_nested <- function(vc, data, nest_vars, effect_level,
                              operative   = FALSE,
                              effect_vars = NULL) {

  n_levels  <- length(nest_vars)
  level_num <- as.integer(gsub("L", "", effect_level))

  residual_var <- vc$vcov[vc$grp == "Residual"]
  residual_var <- if (length(residual_var) == 0) 0 else residual_var[1]

  level_intercept_vars <- sapply(nest_vars, function(v) {
    row <- vc[vc$grp == v & vc$var1 == "(Intercept)" & is.na(vc$var2), ]
    if (nrow(row) == 0) {
      nested_grp <- grep(paste0("^", v, ":"), vc$grp, value = TRUE)
      if (length(nested_grp) > 0) {
        row <- vc[vc$grp == nested_grp[1] &
                    vc$var1 == "(Intercept)" &
                    is.na(vc$var2), ]
      }
    }
    if (nrow(row) == 0) 0 else row$vcov[1]
  })
  names(level_intercept_vars) <- nest_vars

  level_slope_vars <- sapply(nest_vars, function(v) {
    slp_rows <- vc[vc$grp == v &
                     !is.na(vc$var1) &
                     vc$var1 != "(Intercept)" &
                     is.na(vc$var2), ]
    if (nrow(slp_rows) == 0) {
      nested_grp <- grep(paste0("^", v, ":"), vc$grp, value = TRUE)
      if (length(nested_grp) > 0) {
        slp_rows <- vc[vc$grp == nested_grp[1] &
                         !is.na(vc$var1) &
                         vc$var1 != "(Intercept)" &
                         is.na(vc$var2), ]
      }
    }
    if (nrow(slp_rows) == 0) return(0)

    total <- 0
    for (i in seq_len(nrow(slp_rows))) {
      nm <- slp_rows$var1[i]
      if (operative && !is.null(effect_vars) &&
          !slope_matches_effect(nm, effect_vars)) next
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
    idx            <- level_num:n_levels
    variance_error <- sum(level_intercept_vars[idx]) +
      sum(level_slope_vars[idx])
  }

  components <- list(
    residual         = residual_var,
    level_intercepts = level_intercept_vars,
    level_slopes     = level_slope_vars,
    effect_level     = effect_level,
    levels_included  = if (level_num == 1)
      c("residual", nest_vars)
    else
      nest_vars[level_num:n_levels]
  )

  list(variance_error = variance_error, components = components)
}


#' Compute error variance for mixed (crossed-and-nested) designs
#'
#' Residual variance is attributed to the crossed side and subtracted from the
#' nested calculator's output to avoid double-counting.
#'
#' For general effect sizes all components are summed.
#' For operative effect sizes, \code{within_between} (a named vector covering
#' BOTH cross_vars and nest_vars, produced by
#' \code{detect_within_between_mixed}) drives which intercept variances are
#' included in the denominator.
#'
#' @noRd
calc_error_mixed <- function(vc, data, cross_vars, nest_vars, effect_level,
                             operative      = FALSE,
                             within_between = NULL,
                             effect_vars    = NULL) {

  # --- residual ---
  residual_var <- vc$vcov[vc$grp == "Residual"]
  residual_var <- if (length(residual_var) == 0) 0 else residual_var[1]

  # --- crossed intercepts and slopes ---
  cross_int   <- setNames(numeric(length(cross_vars)), cross_vars)
  cross_slope <- setNames(vector("list", length(cross_vars)), cross_vars)

  for (cv in cross_vars) {
    int_row <- vc[vc$grp == cv &
                    vc$var1 == "(Intercept)" &
                    is.na(vc$var2), ]
    cross_int[cv] <- if (nrow(int_row) == 0) 0 else int_row$vcov[1]

    slp_rows <- vc[vc$grp == cv &
                     !is.na(vc$var1) &
                     vc$var1 != "(Intercept)" &
                     is.na(vc$var2), ]

    contrib <- setNames(
      numeric(nrow(slp_rows)),
      if (nrow(slp_rows) > 0) slp_rows$var1 else character(0)
    )
    for (i in seq_len(nrow(slp_rows))) {
      nm         <- slp_rows$var1[i]
      contrib[nm] <- slope_contribution(nm, slp_rows$vcov[i], data)
    }
    cross_slope[[cv]] <- contrib
  }

  # --- nested intercepts and slopes (reuse nested logic) ---
  nested_result <- calc_error_nested(vc, data, nest_vars, effect_level,
                                     operative   = FALSE,   # handle manually below
                                     effect_vars = effect_vars)
  nest_int   <- nested_result$components$level_intercepts
  nest_slope <- nested_result$components$level_slopes
  # Nested calculator already included residual; we'll add it once ourselves.
  nested_variance_no_resid <- nested_result$variance_error - residual_var

  if (operative && !is.null(within_between)) {

    # Operative slopes: only slopes matching the focal effect
    op_cross_slopes <- 0
    for (cv in cross_vars) {
      for (nm in names(cross_slope[[cv]])) {
        if (slope_matches_effect(nm, effect_vars))
          op_cross_slopes <- op_cross_slopes + cross_slope[[cv]][nm]
      }
    }

    # Operative nested slopes from nested calculator
    nested_op <- calc_error_nested(vc, data, nest_vars, effect_level,
                                   operative   = TRUE,
                                   effect_vars = effect_vars)
    op_nest_slopes_and_int <- nested_op$variance_error - residual_var

    # Intercepts: include only for "between" factors
    op_cross_int <- 0
    for (cv in cross_vars) {
      if (identical(within_between[[cv]], "between"))
        op_cross_int <- op_cross_int + cross_int[cv]
    }

    # For nested factors, override the nested operative calculator's intercept
    # handling with our within_between classification.
    # Re-compute nested intercept contribution based on within_between.
    op_nest_int <- 0
    level_num   <- as.integer(gsub("L", "", effect_level))
    n_nest      <- length(nest_vars)
    relevant_nv <- if (level_num == 1) nest_vars else nest_vars[level_num:n_nest]
    for (nv in relevant_nv) {
      if (identical(within_between[[nv]], "between"))
        op_nest_int <- op_nest_int + nest_int[nv]
    }

    # Operative nested slopes (only matching slopes at relevant levels)
    op_nest_slopes <- 0
    for (nv in relevant_nv) {
      slp_rows <- vc[vc$grp == nv &
                       !is.na(vc$var1) &
                       vc$var1 != "(Intercept)" &
                       is.na(vc$var2), ]
      if (nrow(slp_rows) == 0) {
        nested_grp <- grep(paste0("^", nv, ":"), vc$grp, value = TRUE)
        if (length(nested_grp) > 0)
          slp_rows <- vc[vc$grp == nested_grp[1] &
                           !is.na(vc$var1) &
                           vc$var1 != "(Intercept)" &
                           is.na(vc$var2), ]
      }
      for (i in seq_len(nrow(slp_rows))) {
        nm <- slp_rows$var1[i]
        if (slope_matches_effect(nm, effect_vars))
          op_nest_slopes <- op_nest_slopes +
            slope_contribution(nm, slp_rows$vcov[i], data)
      }
    }

    variance_error <- residual_var +
      op_cross_int + op_cross_slopes +
      op_nest_int  + op_nest_slopes

  } else {
    # General: residual + all crossed components + all nested components
    # (residual counted once)
    variance_error <- residual_var +
      sum(cross_int) +
      sum(unlist(cross_slope)) +
      nested_variance_no_resid
  }

  components <- list(
    crossed = list(
      residual          = residual_var,
      factor_intercepts = cross_int,
      factor_slopes     = cross_slope
    ),
    nested = list(
      level_intercepts = nest_int,
      level_slopes     = nest_slope,
      effect_level     = effect_level
    )
  )

  list(variance_error = variance_error, components = components)
}


# Effect-level detection (nested / mixed)

#' Auto-detect the level at which a fixed effect varies in the nesting hierarchy
#' @noRd
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
