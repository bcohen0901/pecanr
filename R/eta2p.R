#' Partial Eta-Squared for Linear Mixed Models
#'
#' Calculates partial eta-squared effect sizes for fixed effects in linear mixed
#' models with crossed, nested, or mixed (crossed-and-nested) random effects.
#'
#' @param model A fitted model object from \code{lme4::lmer()}.
#' @param effect Character string specifying the fixed effect to analyze.
#'   May be given either as the model coefficient name (e.g.
#'   \code{"Group75yr"}) or as the underlying variable name (e.g.
#'   \code{"Group"}); a variable name that maps to a single coefficient is
#'   resolved automatically (with a message). A variable that maps to several
#'   coefficients (a multi-level factor) raises an informative error pointing to
#'   \code{eta2p_omnibus()} (for one factor-level value) or
#'   \code{batch_eta2p()} (for per-coefficient values). Factor predictors and
#'   factor random slopes are read from the model design matrix / model frame,
#'   so no manual recoding to numeric is required.
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
#' @param partial_predictors Logical. If \code{TRUE} (the default), the variance
#'   attributed to a single (non-interaction) predictor is its \emph{unique}
#'   (semipartial) variance -- the variance remaining after removing the part
#'   linearly predictable from the other fixed-effect predictors, equivalently
#'   \eqn{\mathrm{Var}(X) \times \mathrm{tol}(X)} where \eqn{\mathrm{tol}} is the
#'   predictor's tolerance. This yields a measure that reflects only the unique
#'   contribution of the predictor and declines as its redundancy with the
#'   others increases. If \code{FALSE}, the predictor's \emph{total} variance is
#'   used instead (the raw decomposition). For centered, orthogonal designs the
#'   two are identical; they differ only when predictors are correlated. See
#'   Details. Has no effect on interaction terms (whose components are always
#'   centered) or when \code{var_x} is supplied.
#' @param ci Logical. If \code{TRUE}, computes a confidence interval for
#'   \code{eta2p} by parametric bootstrap (\code{lme4::bootMer}): new responses
#'   are simulated from the fitted model, the model is refit, and \code{eta2p}
#'   is recomputed on each replicate, with percentile limits taken from the
#'   resulting distribution. This is the appropriate interval for a ratio of
#'   REML variance components, which has no closed-form sampling distribution.
#'   Each replicate refits the model, so this is opt-in and can be slow on large
#'   models. Default is \code{FALSE}.
#' @param ci_level Numeric confidence level for the bootstrap interval (default
#'   \code{0.95}).
#' @param n_boot Integer number of bootstrap replicates when \code{ci = TRUE}
#'   (default \code{1000}). Use a smaller value (e.g. 100-200) to gauge runtime
#'   on large models before committing to the full count.
#' @param seed Optional integer seed for reproducible bootstrap intervals.
#' @param verbose Logical. If \code{TRUE}, prints detailed results. Default is
#'   \code{FALSE}.
#'
#' @return An object of class \code{"eta2p_lmm"} containing:
#' \item{eta2p}{Partial eta-squared value.}
#' \item{variance_effect}{Variance explained by the effect.}
#' \item{variance_error}{Error variance (denominator).}
#' \item{effect}{Name of the effect.}
#' \item{design}{Design type: \code{"crossed"}, \code{"nested"}, or
#'   \code{"mixed"}.}
#' \item{operative}{Whether operative effect size was calculated.}
#' \item{ci_level, ci_lower, ci_upper, boot_n, boot_method}{(Only when
#'   \code{ci = TRUE}) The confidence level, the lower and upper percentile
#'   bootstrap limits, the number of usable bootstrap replicates, and the
#'   method label (\code{"parametric (bootMer)"}).}
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
#' partial eta-squared in mixed models. The variance a fixed effect explains is
#' \deqn{\mathrm{Var}_{\text{explained}} = b^2 \times \sigma^2_X,}
#' where \eqn{b} is the (partial) fixed-effect coefficient and \eqn{\sigma^2_X}
#' is the variance of the predictor column. Random slope variances are
#' translated to the outcome scale by the same principle,
#' \deqn{\sigma^2_{\text{slope}}(Y) = \sigma^2_b \times \sigma^2_X.}
#' For an interaction, \eqn{\sigma^2_X} is the variance of the product of the
#' \emph{mean-centered} constituent predictors (see "Centering" below). The
#' \code{var_x} argument bypasses this computation when the variance is known a
#' priori.
#' }
#'
#' \subsection{Centering and intended scope}{
#' These quantities are exact when predictors are \strong{mean-centered}, which
#' is standard and recommended practice for the experimental and contrast-coded
#' designs the method primarily targets (sum/contrast coding centers factors
#' automatically). Two points:
#' \itemize{
#'   \item \strong{Interactions: components are centered automatically.} For an
#'     interaction the function mean-centers each \emph{constituent} predictor
#'     and then takes the variance of their product (it does \emph{not} center
#'     the product itself). This makes the interaction effect size invariant to
#'     the location of the constituents: re-centering a predictor by a constant
#'     leaves the model fit and the interaction coefficient unchanged, and now
#'     also leaves the effect size unchanged. (The variance of the \emph{raw}
#'     product depends on the constituents' means and so is not well defined as
#'     an effect size; this is why the centering is applied. The constituent
#'     predictors should still be centered before fitting for the coefficients
#'     themselves to be interpretable, but the effect-size value no longer
#'     depends on it.)
#'   \item \strong{Main effects / random slopes use} \eqn{\sigma^2_X}\strong{,
#'     which equals} \eqn{E[X^2]} \strong{only when} \eqn{X} \strong{is
#'     centered.} The function does not center main-effect predictors (doing so
#'     would not change \eqn{\sigma^2_X = \mathrm{Var}(X)}, which is already
#'     location-invariant), but the random-slope translation
#'     \eqn{\sigma^2_b \sigma^2_X} is exact only for centered predictors;
#'     center continuous predictors before fitting.
#' }
#' Under correlated predictors, each \eqn{b} is a partial coefficient (it
#' controls for the other predictors). By default (\code{partial_predictors =
#' TRUE}) it is multiplied by the predictor's \emph{unique} variance -- the
#' variance remaining after the part linearly predictable from the other
#' predictors is removed -- so the variance attributed to the predictor is
#' \deqn{b^2 \, \mathrm{Var}(X \mid X_{others}) = b^2 \, \mathrm{Var}(X) \,
#'       \mathrm{tol}(X),}
#' where \eqn{\mathrm{tol}(X) = 1 - R^2_{X \sim X_{others}}} is the predictor's
#' tolerance. (The two forms are algebraically identical.) This yields a
#' semipartial variance-explained measure: it reflects only the unique
#' contribution of the predictor and declines as its redundancy with the others
#' increases, matching the unique-variance target estimated by an explained-
#' error (model-comparison) approach.
#'
#' Setting \code{partial_predictors = FALSE} instead uses the predictor's
#' \emph{total} variance \eqn{\sigma^2_X}. The resulting per-predictor values
#' are partial (conditional) effect sizes that do not isolate unique variance
#' when predictors are correlated. For centered, orthogonal designs the two
#' options coincide exactly, because the covariance among predictors is zero, so
#' the choice matters only under collinearity.
#'
#' The tolerance correction is also what makes the interaction handling
#' location-invariant: centering an interaction's components is equivalent to
#' the \eqn{\mathrm{Var} \times \mathrm{tol}} of the product, which is constant
#' across re-centering. The choice is independent of the operative option, which
#' concerns the error denominator rather than the numerator. See Rights & Sterba
#' (2019) on the partial-versus-total distinction.
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
#'
#' The operative option concerns only the denominator; it is independent of the
#' \code{partial_predictors} setting, which controls the numerator. The two
#' compose: an operative effect size uses the same numerator a general one would
#' (unique variance by default), paired with the reduced operative denominator.
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
#' \subsection{Confidence intervals}{
#' With \code{ci = TRUE}, a confidence interval is obtained by parametric
#' bootstrap via \code{lme4::bootMer}: responses are simulated from the fitted
#' model, the model is refit, and \code{eta2p} is recomputed on each replicate;
#' \code{ci_lower} and \code{ci_upper} are the percentile limits. Because
#' \code{eta2p} is a ratio of REML variance components with no closed-form
#' sampling distribution, this bootstrap is the appropriate interval rather than
#' an analytic (e.g. noncentral) formula. Set \code{seed} for reproducibility,
#' and start with a small \code{n_boot} to gauge runtime on large models.
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
#' #  Two crossed factors (backward-compatible call)
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
#' #  Three crossed factors using cross_vars
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
#' #  Mixed design: photos nested within models, crossed with participants
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
#' #  Supply predictor variance directly (var_x)
#' eta2p(model_c, "condition", crossed_data,
#'       design     = "crossed",
#'       cross_vars = c("subject", "item"),
#'       var_x      = 1)   # +/-1 binary predictor: var = 1 by design
#'
#' #  Factor predictor (no manual recoding needed)
#' crossed_data$grp <- factor(rep(c("A", "B"), 60))
#' model_f <- lmer(y ~ grp + (1 | subject) + (1 | item), data = crossed_data)
#' # accepts the variable name "grp"; resolves to the coefficient automatically
#' eta2p(model_f, "grp", crossed_data,
#'       design = "crossed", cross_vars = c("subject", "item"))
#'
#' #  Operative effect size
#' eta2p(model_c, "condition", crossed_data,
#'       design = "crossed", cross_vars = c("subject", "item"),
#'       operative = TRUE)
#'
#' #  Total-variance (raw) numerator instead of the default unique variance
#' eta2p(model_c, "condition", crossed_data,
#'       design = "crossed", cross_vars = c("subject", "item"),
#'       partial_predictors = FALSE)
#'
#' #  Confidence interval by parametric bootstrap
#' eta2p(model_c, "condition", crossed_data,
#'       design = "crossed", cross_vars = c("subject", "item"),
#'       ci = TRUE, n_boot = 200, seed = 1)
#'
#' #  Omnibus effect size for a multi-level factor
#' crossed_data$emo <- factor(rep(c("a", "b", "c", "d"), 30))
#' model_e <- lmer(y ~ emo + (1 | subject) + (1 | item), data = crossed_data)
#' eta2p_omnibus(model_e, "emo", crossed_data,
#'               design = "crossed", cross_vars = c("subject", "item"))
#'
#' #  All fixed effects at once (per-coefficient)
#' batch_eta2p(model_c, crossed_data,
#'             design = "crossed", cross_vars = c("subject", "item"))
#' }
#'
#' @seealso \code{\link{eta2p_omnibus}} for a single factor-level value for a
#'   multi-df factor or interaction (matching an omnibus \emph{F}/chi-square
#'   test); \code{\link{batch_eta2p}} to compute per-coefficient values for all
#'   fixed effects at once.
#'
#' @importFrom stats model.frame var
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
                  partial_predictors = TRUE,
                  ci           = FALSE,
                  ci_level     = 0.95,
                  n_boot       = 1000,
                  seed         = NULL,
                  verbose      = FALSE) {

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

  #  crossed
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

  #  nested
  if (design == "nested") {
    if (is.null(nest_vars))
      stop("For nested designs, 'nest_vars' is required.")
    missing_vars <- nest_vars[!nest_vars %in% colnames(data)]
    if (length(missing_vars) > 0)
      stop("Nesting variable(s) not found in data: ",
           paste(missing_vars, collapse = ", "))
  }

  #  mixed
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
    # The user may have passed a variable name (e.g. "Age") whose model
    # coefficient is contrast-coded (e.g. "Age1" or "Ageyr75"). Resolve it: a
    # factor with k levels produces k-1 coefficients whose names start with the
    # variable name. If exactly one matches, use it; if several, the factor is
    # multi-level and maps to several effects, so guide the user.
    cand <- names(fixed_effects)[startsWith(names(fixed_effects), effect)]
    cand <- setdiff(cand, "(Intercept)")
    # Exclude interaction coefficients (those containing ':'): a bare variable
    # name like "Day" should resolve to its main-effect contrast column, not to
    # "Day1:StimulusType1". Interaction effects must be requested explicitly.
    cand <- cand[!grepl(":", cand)]
    if (length(cand) == 1L) {
      message("Note: interpreting effect '", effect,
              "' as model coefficient '", cand, "'.")
      effect <- cand
    } else if (length(cand) > 1L) {
      stop("Effect '", effect, "' maps to several model coefficients: ",
           paste(cand, collapse = ", "),
           ".\nThis is a multi-level factor; request one coefficient at a time ",
           "(e.g. eta2p(model, \"", cand[1], "\", ...)), or use batch_eta2p().")
    } else {
      stop("Effect '", effect, "' not found in model fixed effects.\n",
           "Available effects: ", paste(names(fixed_effects), collapse = ", "))
    }
  }


  # VARIANCE EXPLAINED BY THE EFFECT

  B           <- fixed_effects[effect]
  effect_vars <- strsplit(effect, ":")[[1]]

  # The model design matrix holds exactly the (possibly contrast-coded,
  # possibly interaction) columns the coefficients refer to. We read predictor
  # variances from it when the effect is not a plain numeric column in `data`,
  # which is what lets factor predictors work without manual recoding. This does
  # not change the estimand: for a numeric predictor present in `data`, the
  # design-matrix column is identical, so var() is unchanged.
  Xmat <- tryCatch(lme4::getME(model, "X"), error = function(e) NULL)

  if (!is.null(var_x)) {
    sigma2_X <- var_x

  } else if (length(effect_vars) == 1) {
    if (effect_vars == "(Intercept)") {
      sigma2_X <- 1
    } else if (effect_vars %in% colnames(data) &&
               is.numeric(data[[effect_vars]])) {
      sigma2_X <- var(data[[effect_vars]], na.rm = TRUE)
    } else if (!is.null(Xmat) && effect %in% colnames(Xmat)) {
      sigma2_X <- var(Xmat[, effect], na.rm = TRUE)
    } else {
      stop("Could not determine the variance of '", effect, "'. ",
           "It is not a numeric column in `data` and no matching column was ",
           "found in the model design matrix. Supply 'var_x' directly, or pass ",
           "the model coefficient name.")
    }

  } else {
    # Interaction. Per the centering rule (Josh, 6-16-26): mean-center each
    # CONSTITUENT predictor, then take the variance of their PRODUCT -- do NOT
    # center the product itself. This makes the interaction effect size
    # invariant to the location (means) of the constituents; the raw product
    # variance is not, because re-centering a constituent by a constant changes
    # Var(X1*X2) while leaving the model fit and the interaction coefficient
    # unchanged. With centered components the interaction variance is well
    # defined and matches the unique-variance (DEE) target.
    #
    # We obtain each constituent's column (centered), multiply, and take var().
    # Constituents may be numeric data columns or factor contrast columns held
    # in the design matrix.
    get_constituent <- function(v) {
      if (v %in% colnames(data) && is.numeric(data[[v]])) return(data[[v]])
      # factor / contrast-coded constituent: find its design-matrix column(s).
      # A constituent of an interaction maps to a main-effect column whose name
      # starts with the variable name (e.g. "Group" -> "Group75yr").
      if (!is.null(Xmat)) {
        if (v %in% colnames(Xmat)) return(Xmat[, v])
        hits <- colnames(Xmat)[startsWith(colnames(Xmat), v) &
                                 !grepl(":", colnames(Xmat))]
        if (length(hits) == 1) return(Xmat[, hits])
      }
      NULL
    }

    cols_list <- lapply(effect_vars, get_constituent)

    if (all(!vapply(cols_list, is.null, logical(1)))) {
      # mean-center each constituent, then form the product
      centered <- lapply(cols_list, function(z) z - mean(z, na.rm = TRUE))
      interaction_product <- Reduce("*", centered)
      sigma2_X <- var(interaction_product, na.rm = TRUE)
    } else if (!is.null(Xmat) && effect %in% colnames(Xmat)) {
      # fallback: constituents not individually resolvable; use the model's
      # product column as-is (assumes the user centered before fitting).
      warning("Interaction '", effect, "': could not resolve all constituents ",
              "to center them; using the model's product column directly. ",
              "Ensure the constituent predictors were mean-centered before ",
              "fitting, or supply 'var_x'.")
      sigma2_X <- var(Xmat[, effect], na.rm = TRUE)
    } else {
      stop("Could not determine the variance of interaction '", effect, "'. ",
           "Supply 'var_x' directly, or ensure the interaction column is in the ",
           "model design matrix.")
    }
  }

  # OPTIONAL: semipartial (residualized) variance for correlated predictors.
  # When partial_predictors = TRUE, replace sigma2_X (the focal predictor's
  # TOTAL variance) with its UNIQUE variance: the variance remaining after the
  # part linearly predictable from the other fixed-effect predictors is removed.
  # The coefficient B is already partial; multiplying it by the residualized
  # variance yields a semipartial variance-explained measure that declines as a
  # predictor's redundancy with the others increases (matching the unique-
  # variance / DEE target). Off by default; for centered, orthogonal designs it
  # is identical to the default. Applied to single (non-interaction) predictors
  # only; for interactions, supply 'var_x' if a residualized product variance is
  # required.
  if (partial_predictors && is.null(var_x) && length(effect_vars) == 1 &&
      effect_vars != "(Intercept)" && !is.null(Xmat)) {
    focal <- if (effect %in% colnames(Xmat)) effect else {
      hits <- colnames(Xmat)[startsWith(colnames(Xmat), effect_vars) &
                               !grepl(":", colnames(Xmat))]
      if (length(hits) == 1) hits else NA_character_
    }
    others <- setdiff(colnames(Xmat), c("(Intercept)", focal))
    if (!is.na(focal) && length(others) > 0) {
      x_focal <- Xmat[, focal]
      Xo      <- Xmat[, others, drop = FALSE]
      # residual of focal regressed on the other predictors
      x_resid <- tryCatch(
        stats::residuals(stats::lm.fit(cbind(1, Xo), x_focal)),
        error = function(e) NULL)
      if (!is.null(x_resid)) {
        sigma2_X <- stats::var(x_resid)
      } else {
        warning("partial_predictors = TRUE: could not residualize '", effect,
                "'; using total variance.")
      }
    }
    # If there are no other predictors, unique == total: sigma2_X unchanged.
  }

  variance_effect <- B^2 * sigma2_X


  # WITHIN/BETWEEN DETECTION (operative only)
  # For mixed designs, every grouping factor (crossed AND nested) is
  # classified; this drives the operative denominator correctly.

  within_between <- NULL

  if (operative) {
    # Detection needs the underlying data variable name(s), but `effect_vars`
    # may hold contrast-coded coefficient names (e.g. "Age1" for factor "Age").
    # Map each back to the data column whose name prefixes it.
    resolve_to_data_var <- function(ev) {
      if (ev %in% colnames(data)) return(ev)
      hit <- colnames(data)[vapply(colnames(data),
                                   function(cn) startsWith(ev, cn), logical(1))]
      # choose the longest matching column name (most specific)
      if (length(hit)) hit[which.max(nchar(hit))] else ev
    }
    detect_vars <- vapply(effect_vars, resolve_to_data_var, character(1))

    if (design == "crossed") {
      within_between <- detect_within_between(data, detect_vars,
                                              subj_var, item_var)
    } else if (design == "mixed") {
      within_between <- detect_within_between_mixed(data, detect_vars,
                                                    cross_vars, nest_vars)
    }
    # nested: operative handled inside calc_error_nested via effect_level
  }


  # ERROR VARIANCE

  if (design == "crossed") {
    result <- calc_error_crossed(vc, data,
                                 cross_vars     = cross_vars,
                                 model          = model,
                                 operative      = operative,
                                 within_between = within_between,
                                 effect_vars    = effect_vars)

  } else if (design == "nested") {
    if (is.null(effect_level))
      effect_level <- detect_effect_level(data, effect, effect_vars, nest_vars)
    result <- calc_error_nested(vc, data, nest_vars, effect_level,
                                model       = model,
                                operative   = operative,
                                effect_vars = effect_vars)

  } else {  # "mixed"
    if (is.null(effect_level))
      effect_level <- detect_effect_level(data, effect, effect_vars, nest_vars)
    result <- calc_error_mixed(vc, data,
                               cross_vars     = cross_vars,
                               nest_vars      = nest_vars,
                               effect_level   = effect_level,
                               model          = model,
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

  # PARAMETRIC BOOTSTRAP CONFIDENCE INTERVAL (opt-in)
  # eta2p is a ratio of REML variance components, which has no closed-form
  # sampling distribution. We therefore obtain a CI by parametric bootstrap:
  # simulate new responses from the fitted model (lme4::bootMer), refit, and
  # recompute eta2p on each replicate, then take percentile limits. This makes
  # no distributional assumption beyond the model itself and respects the actual
  # estimand. It is opt-in because each replicate refits the model.
  if (ci) {
    if (!requireNamespace("lme4", quietly = TRUE))
      stop("Package 'lme4' is required for bootstrap confidence intervals.")
    if (!is.null(seed)) set.seed(seed)

    boot_stat <- function(fitted_model) {
      r <- tryCatch(
        eta2p(fitted_model, effect, data,
              design       = design,
              subj_var     = subj_var,
              item_var     = item_var,
              cross_vars   = cross_vars,
              nest_vars    = nest_vars,
              effect_level = effect_level,
              var_x        = var_x,
              operative    = operative,
              partial_predictors = partial_predictors,
              ci           = FALSE,          # prevent recursion
              verbose      = FALSE),
        error = function(e) NULL)
      if (is.null(r)) NA_real_ else as.numeric(r$eta2p)
    }

    bb <- tryCatch(
      suppressWarnings(lme4::bootMer(model, boot_stat, nsim = n_boot,
                                     type = "parametric",
                                     use.u = FALSE)),
      error = function(e) NULL)

    if (!is.null(bb)) {
      vals  <- bb$t[, 1]
      vals  <- vals[is.finite(vals)]
      alpha <- 1 - ci_level
      qs    <- stats::quantile(vals, c(alpha / 2, 1 - alpha / 2),
                               na.rm = TRUE, names = FALSE)
      output$ci_level   <- ci_level
      output$ci_lower   <- qs[1]
      output$ci_upper   <- qs[2]
      output$boot_n     <- length(vals)
      output$boot_method <- "parametric (bootMer)"
    } else {
      output$ci_lower <- NA_real_
      output$ci_upper <- NA_real_
      warning("Bootstrap failed; confidence interval not available.")
    }
  }

  if (verbose) print(output)

  return(output)
}


#' Omnibus (factor-level) Partial Eta-Squared for a Multi-df Effect
#'
#' For a multi-level factor or a multi-df interaction, computes a single
#' factor-level partial eta-squared corresponding to the omnibus test of that
#' effect (e.g. the chi-square or F test on all of its degrees of freedom),
#' rather than one value per contrast. The variance attributed to the effect is
#' the variance of the summed fitted contribution of all the design-matrix
#' columns belonging to the effect (which correctly includes the covariances
#' among those columns); the error denominator is the same one eta2p() uses.
#'
#' @param model A fitted lme4 model.
#' @param effect The variable/effect name. For a factor main effect, the bare
#'   variable name (e.g. "rating_type"). For an interaction, the colon form of
#'   the term as it appears in the design matrix (e.g. "group:rating_type").
#' @param data The data frame used to fit the model.
#' @param ... Passed to eta2p() to obtain the error denominator (design,
#'   cross_vars, subj_var, item_var, nest_vars, operative, etc.).
#'
#' @return An object of class \code{"eta2p_omnibus"}: a list containing
#' \item{eta2p}{The omnibus (factor-level) partial eta-squared.}
#' \item{variance_effect}{Variance of the summed fitted contribution of all the
#'   effect's design-matrix columns (includes covariances among them).}
#' \item{variance_error}{Error variance denominator (from \code{eta2p()}).}
#' \item{effect}{The effect name requested.}
#' \item{coefficients}{The design-matrix column names aggregated.}
#' \item{n_df}{Number of columns aggregated (the effect's degrees of freedom;
#'   matches the df of the corresponding omnibus test).}
#' \item{type}{\code{"omnibus"}.}
#'
#' @details
#' \code{batch_eta2p()} and a per-coefficient call to \code{eta2p()} return one
#' value per contrast (per design-matrix column). For a multi-level factor or a
#' multi-df interaction, none of those single values corresponds to the omnibus
#' test of the whole effect (the \emph{F} or chi-square test on all of its
#' degrees of freedom). This function fills that gap.
#'
#' The omnibus variance explained is computed as the variance of the
#' \emph{summed} fitted contribution of all the effect's columns,
#' \eqn{\mathrm{Var}(\sum_j b_j x_j)}, which correctly includes the covariances
#' among the dummy/contrast columns. Simply adding the per-coefficient
#' numerators would drop those covariances and is incorrect. The error
#' denominator is the same one \code{eta2p()} uses, so the result remains a
#' partial eta-squared (error-variance denominator), \emph{not} a total-variance
#' \eqn{R^2}; expect it to differ from \code{r2glmm::r2beta()}, which targets the
#' latter.
#'
#' @examples
#' \donttest{
#' library(lme4)
#' set.seed(1)
#' d <- expand.grid(subject = factor(1:30), item = factor(1:12),
#'                  emo = factor(c("a", "b", "c", "d")))
#' d$y <- as.integer(d$emo) + rnorm(30)[d$subject] + rnorm(nrow(d))
#' m <- lmer(y ~ emo + (1 | subject) + (1 | item), data = d)
#' # one factor-level value to pair with the omnibus test of `emo`:
#' eta2p_omnibus(m, "emo", d, design = "crossed",
#'               cross_vars = c("subject", "item"))
#' }
#'
#' @seealso \code{\link{eta2p}} for single-coefficient values and
#'   \code{\link{batch_eta2p}} for all coefficients at once.
#'
#' @export
eta2p_omnibus <- function(model, effect, data, ...) {
  X <- lme4::getME(model, "X")
  b <- lme4::fixef(model)

  # Identify the design-matrix columns belonging to this effect.
  if (grepl(":", effect)) {
    # interaction: match columns whose set of ':'-separated variable stems equals
    # the requested term's stems (order-independent), so "group:rating_type"
    # collects group?:rating_type? dummy columns.
    want <- sort(strsplit(effect, ":")[[1]])
    col_is_term <- vapply(colnames(X), function(cn) {
      if (!grepl(":", cn)) return(FALSE)
      parts <- strsplit(cn, ":")[[1]]
      # strip trailing factor-level text by matching each requested stem as a prefix
      all(vapply(want, function(w) any(startsWith(parts, w)), logical(1))) &&
        length(parts) == length(want)
    }, logical(1))
    cols <- colnames(X)[col_is_term]
  } else {
    cols <- colnames(X)[startsWith(colnames(X), effect) & !grepl(":", colnames(X))]
    cols <- setdiff(cols, "(Intercept)")
  }

  if (length(cols) == 0)
    stop("No design-matrix columns found for effect '", effect, "'. ",
         "Available columns: ", paste(colnames(X), collapse = ", "))

  # Omnibus variance explained = variance of the summed fitted contribution of
  # all columns for this effect (includes their covariances).
  contrib  <- as.numeric(X[, cols, drop = FALSE] %*% b[cols])
  var_eff  <- stats::var(contrib)

  # Reuse eta2p()'s error denominator (identical across the effect's columns).
  ref <- eta2p(model, cols[1], data, ..., verbose = FALSE)
  var_err <- ref$variance_error

  out <- list(
    eta2p           = var_eff / (var_eff + var_err),
    variance_effect = var_eff,
    variance_error  = var_err,
    effect          = effect,
    coefficients    = cols,
    n_df            = length(cols),
    type            = "omnibus"
  )
  class(out) <- c("eta2p_omnibus", "list")
  out
}


#' Batch Calculate Partial Eta-Squared for Multiple Effects
#'
#' Calculates partial eta-squared for all fixed effects (excluding the
#' intercept) in a model, returning one row per model coefficient. For a
#' multi-level factor or multi-df interaction this yields one value
#' \emph{per contrast} (per design-matrix column), not a single factor-level
#' value; for the omnibus effect size of such a term, use
#' \code{\link{eta2p_omnibus}}.
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
                        partial_predictors = TRUE,
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
                   partial_predictors = partial_predictors,
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
  x <- data[[effect_var]]
  # A predictor "varies within" a grouping factor if, within at least one level
  # of that factor, the predictor takes more than one distinct value. This works
  # for factors, characters, and numerics alike (var() is undefined on factors).
  xnum <- if (is.factor(x) || is.character(x)) as.integer(as.factor(x)) else x
  n_distinct_within <- tapply(xnum, data[[group_var]],
                              function(z) length(unique(z[!is.na(z)])))
  any(n_distinct_within > 1, na.rm = TRUE)
}


# Slope helpers

#' Translate a random slope variance component to the outcome scale
#' (sigma2_b * sigma2_X).  Returns 0 and warns if variables are missing.
#' @noRd
slope_contribution <- function(slope_var_name, vc_row_vcov, data, model = NULL) {

  # Resolve the predictor column for this random slope. lme4 names a factor's
  # random slope with its contrast-coded name (e.g. "Ageyr75"), which is NOT a
  # column in `data`; it lives in the model design matrix. We look there first,
  # falling back to the data frame for plain numeric predictors. This mirrors the
  # fixed-effect variance lookup and is what makes factor random slopes work.
  Xmat <- if (!is.null(model)) tryCatch(lme4::getME(model, "X"),
                                        error = function(e) NULL) else NULL

  lookup_var <- function(nm) {
    if (!is.null(Xmat) && nm %in% colnames(Xmat))
      return(var(Xmat[, nm], na.rm = TRUE))
    if (nm %in% colnames(data) && is.numeric(data[[nm]]))
      return(var(data[[nm]], na.rm = TRUE))
    # Factor random slope: lme4 may name it with a level (e.g. "Ageyr75") and
    # code it as a dummy, independent of the fixed-effect contrast. Reconstruct
    # that dummy column from the model frame: find the factor whose name prefixes
    # `nm` and whose level completes it.
    if (!is.null(model)) {
      mf <- tryCatch(model.frame(model), error = function(e) NULL)
      if (!is.null(mf)) {
        for (v in names(mf)) {
          col <- mf[[v]]
          if (is.factor(col) && startsWith(nm, v)) {
            lev <- sub(paste0("^", v), "", nm)
            if (lev %in% levels(col))
              return(var(as.numeric(col == lev), na.rm = TRUE))
          }
        }
      }
    }
    NA_real_
  }

  if (grepl(":", slope_var_name)) {
    # Interaction random slope. The RE term name (e.g. "Dayd2:StimulusTypes2")
    # is typically dummy-coded and will not match the (possibly sum-coded)
    # fixed-effect product column in X. Resolve each component with the same
    # logic used for single slopes (design matrix, data, or reconstructed factor
    # dummy), then take the variance of their product. Components are mean-
    # centered first, consistent with the interaction handling for fixed effects.
    parts <- strsplit(slope_var_name, ":")[[1]]
    cols  <- lapply(parts, function(p) {
      if (!is.null(Xmat) && p %in% colnames(Xmat)) return(Xmat[, p])
      if (p %in% colnames(data) && is.numeric(data[[p]])) return(data[[p]])
      # reconstruct factor dummy from the model frame (e.g. "Dayd2" -> Day=="d2")
      if (!is.null(model)) {
        mf <- tryCatch(model.frame(model), error = function(e) NULL)
        if (!is.null(mf)) for (v in names(mf)) {
          col <- mf[[v]]
          if (is.factor(col) && startsWith(p, v)) {
            lev <- sub(paste0("^", v), "", p)
            if (lev %in% levels(col)) return(as.numeric(col == lev))
          }
        }
      }
      NULL
    })
    if (any(vapply(cols, is.null, logical(1)))) {
      warning("Interaction slope '", slope_var_name,
              "' components could not be resolved - skipping.")
      return(0)
    }
    centered <- lapply(cols, function(z) z - mean(z, na.rm = TRUE))
    product  <- Reduce("*", centered)
    sigma2_X <- var(product, na.rm = TRUE)
  } else {
    sigma2_X <- lookup_var(slope_var_name)
    if (is.na(sigma2_X)) {
      warning("Slope variable '", slope_var_name,
              "' not found in model or data - skipping.")
      return(0)
    }
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
calc_error_crossed <- function(vc, data, cross_vars, model = NULL,
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
      contributions[nm] <- slope_contribution(nm, slp_rows$vcov[i], data, model = model)
    }

    factor_slope_vars[[cv]] <- contributions
  }

  if (operative && !is.null(within_between) && !is.null(effect_vars)) {

    # Operative denominator (per the general/operative definition): keep the
    # residual, ALL random-slope contributions, and the intercept variance only
    # of factors the predictor varies BETWEEN. Intercepts of factors the
    # predictor varies WITHIN are dropped, because those between-group
    # differences are differenced out by the within-factor structure of the
    # effect's test and so do not enter its standard error.
    all_slope_total <- sum(unlist(factor_slope_vars))

    variance_error <- residual_var + all_slope_total

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
calc_error_nested <- function(vc, data, nest_vars, effect_level, model = NULL,
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
      total <- total + slope_contribution(nm, slp_rows$vcov[i], data, model = model)
    }
    total
  })
  names(level_slope_vars) <- nest_vars

  if (operative) {
    # Operative: keep residual + ALL slopes, plus intercepts only of nesting
    # levels the predictor varies BETWEEN. A level the predictor varies WITHIN
    # has its intercept dropped (its between-group differences are differenced
    # out by the within-level structure of the effect's test). This mirrors the
    # crossed-design operative rule.
    variance_error <- residual_var + sum(level_slope_vars)
    for (v in nest_vars) {
      varies_within <- check_within_factor(data, effect_vars[1], v)
      if (!varies_within)
        variance_error <- variance_error + level_intercept_vars[v]
    }
  } else if (level_num == 1) {
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
calc_error_mixed <- function(vc, data, cross_vars, nest_vars, effect_level, model = NULL,
                             operative      = FALSE,
                             within_between = NULL,
                             effect_vars    = NULL) {

  #  residual
  residual_var <- vc$vcov[vc$grp == "Residual"]
  residual_var <- if (length(residual_var) == 0) 0 else residual_var[1]

  #  crossed intercepts and slopes
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
      contrib[nm] <- slope_contribution(nm, slp_rows$vcov[i], data, model = model)
    }
    cross_slope[[cv]] <- contrib
  }

  #  nested intercepts and slopes (reuse nested logic)
  nested_result <- calc_error_nested(vc, data, nest_vars, effect_level,
                                     model       = model,
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
            slope_contribution(nm, slp_rows$vcov[i], data, model = model)
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
