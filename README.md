
# pecanr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/pecanr)](https://CRAN.R-project.org/package=pecanr)
[![R-CMD-check](https://github.com/bcohen0901/pecanr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bcohen0901/pecanr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**pecanr** computes partial eta-squared (eta2p) effect sizes for fixed
effects in linear mixed models fitted with `lme4`. It correctly handles
crossed and nested random effects structures – including random slopes –
using a variance decomposition approach that translates slope variances
to the outcome scale.

## Why pecanr?

**pecanr** accounts for:

- **Crossed designs** – any number of grouping factors (subjects, items,
  raters, etc.)
- **Nested designs** – hierarchical structures with automatic level
  detection
- **Random slopes** – translated to the outcome scale via sigma^2_slope
  x sigma^2_X
- **Operative effect sizes** – excluding variance components that don’t
  contribute to the standard error of the tested effect

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

## References

Correll, J., Mellinger, C., McClelland, G. H., & Judd, C. M. (2020).
Avoid Cohen’s “Small”, “Medium”, and “Large” for Power Analysis. *Trends
in Cognitive Sciences*, 24(3), 200-207.

Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
multilevel models: An integrative framework for defining R-squared
measures. *Psychological Methods*, 24(3), 309-338.
