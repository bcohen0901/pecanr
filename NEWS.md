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
  σ²_slope × σ²_X, correctly accounting for predictor scaling and
  interaction terms.
