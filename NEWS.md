#  pecanr 0.1.0

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
