# pecanr 0.3.0

## New features

* New function `eta2p_omnibus()` computes a single factor-level partial
  eta-squared for a multi-level factor or multi-df interaction, corresponding
  to the omnibus (multi-df) test of that effect rather than to individual
  contrasts. The variance attributed to the effect is the variance of the
  summed fitted contribution of all of the effect's design-matrix columns,
  which correctly includes the covariances among them; the error denominator
  matches `eta2p()`. Use this when pairing an effect size with an omnibus F or
  chi-square test.
* `eta2p()` and `batch_eta2p()` now accept factor predictors directly. The
  focal effect may be given either as a model coefficient name (e.g.
  `"Group75yr"`) or as the underlying variable name (e.g. `"Group"`); a
  variable that maps to a single coefficient is resolved automatically, and a
  multi-level factor raises an informative error pointing to `eta2p_omnibus()`
  or `batch_eta2p()`. Predictor variances and factor random slopes are read
  from the model design matrix and model frame, so manual recoding of factors
  to numeric is no longer required.
* `eta2p()` gains an optional confidence interval via parametric bootstrap
  (`ci = TRUE`), with `ci_level`, `n_boot`, and `seed` arguments. Because the
  effect size is a ratio of REML variance components with no closed-form
  sampling distribution, the interval is obtained by simulating from the fitted
  model (`lme4::bootMer`), refitting, and recomputing the effect size on each
  replicate. This is opt-in, as each replicate refits the model.
* `eta2p()` and `batch_eta2p()` gain a `partial_predictors` argument. When
  `TRUE`, the variance attributed to a single (non-interaction) predictor is its
  unique (semipartial) variance -- residualized on the other fixed-effect
  predictors -- rather than its total variance, yielding a measure that declines
  with a predictor's redundancy under correlation. Defaults to `FALSE`; for
  centered, orthogonal designs the two are identical.

## Changes

* Interaction effect sizes are now computed from the mean-centered constituent
  predictors: each component of an interaction is centered before forming the
  product whose variance enters the numerator (the product itself is not
  centered). This makes the interaction effect size invariant to the location
  (means) of its constituent predictors, which the variance of the raw product
  is not. Effect sizes for main effects are unaffected, as is the value for
  interactions whose components were already centered (e.g. contrast-coded
  factors). Note that for models with uncentered continuous interaction
  components, interaction effect sizes may differ from those produced by
  earlier versions.
* The default for `verbose` in `eta2p()` is now `FALSE` (previously `TRUE`),
  matching `batch_eta2p()`. Pass `verbose = TRUE` to restore printed output.

## Bug fixes

* Operative effect sizes (`operative = TRUE`) now work correctly for factor
  predictors. Within/between classification previously matched the focal
  effect's name directly against the data columns, which failed for a factor
  whose coefficient name (e.g. `"Age75"`) differs from its variable name
  (e.g. `"Age"`); the grouping factor could then be misclassified and the
  wrong intercept variances included in or excluded from the operative
  denominator. Coefficient names are now resolved back to their underlying
  variable before classification.
* Fixed a bug in which the random slope of an interaction term (e.g. a random
  slope on `x1:x2`) was silently omitted from the error denominator when its
  random-effect name did not match the fixed-effect product column. The
  contribution of interaction random slopes is now resolved from the model
  frame and included.

# pecanr 0.2.0

## New features

* `eta2p()` and `batch_eta2p()` now support `design = "mixed"` for models with
  both crossed and nested random effects simultaneously. A canonical example is
  participants viewing multiple photos of each model: photos are nested within
  models, but both levels are crossed with participants. Supply both
  `cross_vars` and `nest_vars` to use this design.

## Bug fixes

* Fixed a bug in operative effect size calculation for crossed designs.
  `detect_within_between()` previously used hardcoded `$subj` and `$item` keys
  to classify grouping factors as within or between, which caused intercept
  variances to be silently omitted from the operative denominator. Keys are now
  indexed by actual variable name, so the correct components are always
  included.

## Breaking changes

* `batch_eta2p()` output columns for within/between status are now named
  `within_<varname>` (e.g., `within_participant`, `within_item`) rather than
  the hardcoded `within_subj` and `within_item`. Code that references these
  columns by name will need to be updated.
* Operative effect sizes with 3 or more crossed factors now correctly gate each
  factor's intercept variance on its within/between status. Previously, third
  and higher factors were always included in the operative denominator
  regardless of whether the effect was within or between those factors.

# pecanr 0.1.2

* Initial CRAN release.
* `eta2p()` computes partial eta-squared for a single fixed effect in a fitted
  `lmer` model, supporting crossed and nested random effects structures.
* `batch_eta2p()` computes partial eta-squared for all fixed effects in a model.
* Crossed designs support any number of grouping factors via `cross_vars`.
* Nested designs support automatic effect-level detection.
* Operative effect sizes available via `operative = TRUE`.
* Random slope variances are translated to the outcome scale using
  sigma^2_slope x sigma^2_X, correctly accounting for predictor scaling and
  interaction terms.
