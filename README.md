# pecanr

<!-- badges: start -->
[![CRAN
status](https://www.r-pkg.org/badges/version/pecanr)](https://CRAN.R-project.org/package=pecanr)
[![R-CMD-check](https://github.com/bcohen0901/pecanr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bcohen0901/pecanr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**pecanr** computes partial eta-squared (eta2p) effect sizes for fixed
effects in linear mixed models fitted with `lme4`. It correctly handles
crossed, nested, and mixed (crossed-and-nested) random effects
structures – including random slopes – using a variance decomposition
approach that translates slope variances to the outcome scale.

## Why pecanr?

**pecanr** accounts for:

- **Crossed designs** – any number of grouping factors (subjects, items,
  raters, etc.)
- **Nested designs** – hierarchical structures with automatic level
  detection
- **Mixed designs** – grouping factors that are simultaneously nested
  within some variables and crossed with others (e.g., photos nested
  within models, but both crossed with participants)
- **Random slopes** – translated to the outcome scale via sigma^2_slope
  x sigma^2_X
- **Factor predictors** – supplied directly, by variable name or
  coefficient name, with no manual recoding to numeric
- **Omnibus effect sizes** – a single factor-level value for a
  multi-level factor or multi-df interaction, matching the omnibus
  F or chi-square test
- **Operative effect sizes** – excluding variance components that don’t
  contribute to the standard error of the tested effect
- **Unique (semipartial) variance** – optional residualization for
  correlated predictors via `partial_predictors = TRUE`
- **Confidence intervals** – optional, via parametric bootstrap

## Installation

You can install the development version of pecanr from
[GitHub](https://github.com/bcohen0901/pecanr):

``` r
# install.packages("pak")
pak::pak("bcohen0901/pecanr")
```

Once on CRAN:

``` r
install.packages("pecanr")
```

## Usage

### Crossed design (subjects x items)

``` r
library(lme4)
library(pecanr)

model <- lmer(y ~ condition + (1 | subject) + (1 | item), data = my_data)

eta2p(model, effect = "condition", data = my_data,
      design     = "crossed",
      cross_vars = c("subject", "item"))
```

### Three crossed factors

``` r
model3 <- lmer(y ~ condition + (1 | subject) + (1 | item) + (1 | rater),
               data = my_data)

eta2p(model3, effect = "condition", data = my_data,
      design     = "crossed",
      cross_vars = c("subject", "item", "rater"))
```

### Nested design

``` r
model_nested <- lmer(y ~ treatment + (1 | school/class), data = my_data)

eta2p(model_nested, effect = "treatment", data = my_data,
      design    = "nested",
      nest_vars = c("class", "school"))
```

### Mixed design (nested-and-crossed)

Use `design = "mixed"` when some grouping factors are nested within
others but all levels are crossed with additional factors. A common
example is participants viewing multiple photos of each model: photos
are nested within models, but both levels are crossed with participants.

``` r
model_mixed <- lmer(y ~ x + (1 | participant) + (1 | model) + (1 | photo:model),
                    data = my_data)

eta2p(model_mixed, effect = "x", data = my_data,
      design     = "mixed",
      cross_vars = "participant",
      nest_vars  = c("photo", "model"))
```

### Factor predictors

Factor predictors work directly – pass either the variable name or a
specific coefficient name. A factor with more than two levels maps to
several coefficients, so request one contrast at a time, use
`eta2p_omnibus()` for the whole factor, or use `batch_eta2p()`.

``` r
model_f <- lmer(rating ~ group + (1 | subject) + (1 | item), data = my_data)

# variable name is resolved to its coefficient automatically
eta2p(model_f, effect = "group", data = my_data,
      design     = "crossed",
      cross_vars = c("subject", "item"))
```

### Omnibus effect size for a multi-level factor or interaction

`eta2p_omnibus()` returns a single factor-level partial eta-squared for
a multi-df effect, corresponding to its omnibus F or chi-square test
(rather than to individual contrasts).

``` r
model_e <- lmer(rating ~ group * emotion + (1 | subject) + (1 | item),
                data = my_data)

eta2p_omnibus(model_e, effect = "emotion", data = my_data,
              design     = "crossed",
              cross_vars = c("subject", "item"))

# interactions: give the term in colon form
eta2p_omnibus(model_e, effect = "group:emotion", data = my_data,
              design     = "crossed",
              cross_vars = c("subject", "item"))
```

### Batch over all effects

``` r
batch_eta2p(model, data = my_data,
            design     = "crossed",
            cross_vars = c("subject", "item"))
```

### Operative effect sizes

``` r
eta2p(model, effect = "condition", data = my_data,
      design     = "crossed",
      cross_vars = c("subject", "item"),
      operative  = TRUE)
```

### Unique (semipartial) variance for correlated predictors

By default the variance attributed to a predictor uses its total variance,
so the per-predictor effect sizes are partial (conditional) and are exact
for orthogonal designs but do not partition variance when predictors are
correlated. Set `partial_predictors = TRUE` to instead use each predictor's
unique variance (residualized on the others), which declines as redundancy
increases. For centered, orthogonal designs the two are identical.

``` r
eta2p(model, effect = "x1", data = my_data,
      design     = "crossed",
      cross_vars = c("subject", "item"),
      partial_predictors = TRUE)
```

### Confidence intervals

A parametric bootstrap confidence interval is available via `ci = TRUE`.
Because the effect size is a ratio of REML variance components, the
interval is obtained by simulating from the fitted model, refitting, and
recomputing the effect size on each replicate. This refits the model
`n_boot` times, so it is opt-in and can be slow on large models.

``` r
eta2p(model, effect = "condition", data = my_data,
      design     = "crossed",
      cross_vars = c("subject", "item"),
      ci         = TRUE,
      n_boot     = 1000,
      seed       = 1)
```

### A note on interactions and centering

For an interaction, the variance entering the numerator is computed from
the **mean-centered** constituent predictors automatically (the product
itself is not centered). As a result the interaction effect size is
invariant to the location of its constituent predictors – you do not
need to center them yourself for the effect size to be well defined.
Centering continuous predictors before fitting is still good practice
for a separate reason: it makes the lower-order (main-effect)
coefficients interpretable in models that contain interactions.

## References

Correll, J., Mellinger, C., McClelland, G. H., & Judd, C. M. (2020).
Avoid Cohen’s “Small”, “Medium”, and “Large” for Power Analysis. *Trends
in Cognitive Sciences*, 24(3), 200-207.

Correll, J., Mellinger, C., & Pedersen, E. J. (2022). Flexible
approaches for estimating partial eta squared in mixed-effects models
with crossed random factors. *Behavior Research Methods*, 54, 1626-1642.

Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
multilevel models: An integrative framework for defining R-squared
measures. *Psychological Methods*, 24(3), 309-338.
