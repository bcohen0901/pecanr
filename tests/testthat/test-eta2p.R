library(testthat)
library(pecanr)
library(lme4)


# SHARED TEST DATA


# Crossed design: subjects x items
set.seed(42)
n_subj <- 30
n_item <- 20

crossed_data <- expand.grid(
  subject   = factor(1:n_subj),
  item      = factor(1:n_item)
)
crossed_data$condition <- ifelse(as.integer(crossed_data$item) <= 10, -1, 1)
crossed_data$y <- 2 * crossed_data$condition +
  rnorm(n_subj * n_item, 0, 1) +
  rep(rnorm(n_subj, 0, 0.5), each = n_item) +
  rep(rnorm(n_item, 0, 0.5), times = n_subj)

# Nested design: students within classrooms within schools
set.seed(42)
nested_data <- data.frame(
  school    = factor(rep(1:5, each = 40)),
  class     = factor(rep(1:20, each = 10)),
  student   = factor(1:200),
  treatment = rep(c(-1, 1), 100)
)
nested_data$y <- 1.5 * nested_data$treatment +
  rep(rnorm(5,  0, 0.8), each = 40) +
  rep(rnorm(20, 0, 0.4), each = 10) +
  rnorm(200, 0, 1)

# Fitted models
model_crossed <- lmer(y ~ condition + (1 | subject) + (1 | item),
                      data = crossed_data, REML = TRUE)

model_nested <- lmer(y ~ treatment + (1 | school) + (1 | class),
                     data = nested_data, REML = TRUE)

model_slopes <- suppressWarnings(
  lmer(y ~ condition + (1 + condition | subject) + (1 | item),
       data = crossed_data, REML = TRUE)
)


# INPUT VALIDATION


test_that("errors on non-lmer model", {
  lm_model <- lm(y ~ condition, data = crossed_data)
  expect_error(
    eta2p(lm_model, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "item")),
    "lmerMod"
  )
})

test_that("errors on missing effect name", {
  expect_error(
    eta2p(model_crossed, "nonexistent", crossed_data,
          design = "crossed", cross_vars = c("subject", "item")),
    "not found in model fixed effects"
  )
})

test_that("errors on missing cross_vars for crossed design", {
  expect_error(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed"),
    "cross_vars"
  )
})

test_that("errors on missing nest_vars for nested design", {
  expect_error(
    eta2p(model_nested, "treatment", nested_data,
          design = "nested"),
    "nest_vars"
  )
})

test_that("errors on invalid var_x", {
  expect_error(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "item"),
          var_x = -1),
    "non-negative"
  )
  expect_error(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "item"),
          var_x = c(1, 2)),
    "single"
  )
})

test_that("errors on cross_vars not found in data", {
  expect_error(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "nonexistent")),
    "not found in data"
  )
})

test_that("errors on nest_vars not found in data", {
  expect_error(
    eta2p(model_nested, "treatment", nested_data,
          design = "nested", nest_vars = c("class", "nonexistent")),
    "not found in data"
  )
})


# RETURN VALUE


test_that("eta2p returns correct class", {
  result <- eta2p(model_crossed, "condition", crossed_data,
                  design = "crossed", cross_vars = c("subject", "item"),
                  verbose = FALSE)
  expect_s3_class(result, "eta2p_lmm")
})

test_that("eta2p value is between 0 and 1", {
  result <- eta2p(model_crossed, "condition", crossed_data,
                  design = "crossed", cross_vars = c("subject", "item"),
                  verbose = FALSE)
  expect_gte(result$eta2p, 0)
  expect_lte(result$eta2p, 1)
})

test_that("eta2p output contains required fields", {
  result <- eta2p(model_crossed, "condition", crossed_data,
                  design = "crossed", cross_vars = c("subject", "item"),
                  verbose = FALSE)
  expect_true(all(c("eta2p", "variance_effect", "variance_error",
                    "effect", "design", "operative",
                    "variance_components") %in% names(result)))
})

test_that("variance_effect + variance_error equals denominator", {
  result <- eta2p(model_crossed, "condition", crossed_data,
                  design = "crossed", cross_vars = c("subject", "item"),
                  verbose = FALSE)
  expected_eta2p <- result$variance_effect /
    (result$variance_effect + result$variance_error)
  expect_equal(result$eta2p, expected_eta2p, tolerance = 1e-10)
})


# CROSSED DESIGN


test_that("crossed design runs with two factors", {
  expect_no_error(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "item"),
          verbose = FALSE)
  )
})

test_that("backward-compatible subj_var/item_var interface works", {
  result_new <- eta2p(model_crossed, "condition", crossed_data,
                      design = "crossed",
                      cross_vars = c("subject", "item"),
                      verbose = FALSE)
  result_old <- eta2p(model_crossed, "condition", crossed_data,
                      design = "crossed",
                      subj_var = "subject", item_var = "item",
                      verbose = FALSE)
  expect_equal(unname(result_new$eta2p), unname(result_old$eta2p))
})

test_that("crossed design with three factors runs", {
  set.seed(42)
  three_data <- expand.grid(
    subject = factor(1:20),
    item    = factor(1:10),
    rater   = factor(1:3)
  )
  three_data$condition <- ifelse(as.integer(three_data$item) <= 5, -1, 1)
  three_data$y <- 1.5 * three_data$condition +
    rnorm(nrow(three_data), 0, 1)

  model3 <- suppressWarnings(
    lmer(y ~ condition + (1 | subject) + (1 | item) + (1 | rater),
         data = three_data, REML = TRUE)
  )

  expect_no_error(
    suppressWarnings(
      eta2p(model3, "condition", three_data,
            design = "crossed",
            cross_vars = c("subject", "item", "rater"),
            verbose = FALSE)
    )
  )
})

test_that("var_x override gives same result as computed variance", {
  computed_var <- var(crossed_data$condition)
  result_auto <- eta2p(model_crossed, "condition", crossed_data,
                       design = "crossed", cross_vars = c("subject", "item"),
                       verbose = FALSE)
  result_manual <- eta2p(model_crossed, "condition", crossed_data,
                         design = "crossed", cross_vars = c("subject", "item"),
                         var_x = computed_var, verbose = FALSE)
  expect_equal(unname(result_auto$eta2p), unname(result_manual$eta2p),
               tolerance = 1e-10)
})

test_that("model with random slopes runs without error", {
  expect_no_error(
    suppressWarnings(
      eta2p(model_slopes, "condition", crossed_data,
            design = "crossed", cross_vars = c("subject", "item"),
            verbose = FALSE)
    )
  )
})

test_that("random slopes increase error variance", {
  result_no_slopes <- eta2p(model_crossed, "condition", crossed_data,
                            design = "crossed",
                            cross_vars = c("subject", "item"),
                            verbose = FALSE)
  result_slopes <- suppressWarnings(
    eta2p(model_slopes, "condition", crossed_data,
          design = "crossed",
          cross_vars = c("subject", "item"),
          verbose = FALSE)
  )
  expect_gte(result_slopes$variance_error, result_no_slopes$variance_error)
})


# NESTED DESIGN


test_that("nested design runs without error", {
  expect_no_error(
    eta2p(model_nested, "treatment", nested_data,
          design = "nested", nest_vars = c("class", "school"),
          verbose = FALSE)
  )
})

test_that("nested design eta2p is between 0 and 1", {
  result <- eta2p(model_nested, "treatment", nested_data,
                  design = "nested", nest_vars = c("class", "school"),
                  verbose = FALSE)
  expect_gte(result$eta2p, 0)
  expect_lte(result$eta2p, 1)
})

test_that("nested design output contains effect_level", {
  result <- eta2p(model_nested, "treatment", nested_data,
                  design = "nested", nest_vars = c("class", "school"),
                  verbose = FALSE)
  expect_true("effect_level" %in% names(result))
  expect_match(result$effect_level, "^L[0-9]+$")
})


# OPERATIVE EFFECT SIZES


test_that("operative = TRUE runs without error", {
  expect_no_error(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "item"),
          operative = TRUE, verbose = FALSE)
  )
})

test_that("operative eta2p >= general eta2p for within-subjects effect", {
  result_general   <- eta2p(model_crossed, "condition", crossed_data,
                            design = "crossed",
                            cross_vars = c("subject", "item"),
                            operative = FALSE, verbose = FALSE)
  result_operative <- eta2p(model_crossed, "condition", crossed_data,
                            design = "crossed",
                            cross_vars = c("subject", "item"),
                            operative = TRUE, verbose = FALSE)
  expect_gte(result_operative$eta2p, result_general$eta2p)
})


# BATCH FUNCTION


test_that("batch_eta2p returns a data frame", {
  result <- batch_eta2p(model_crossed, crossed_data,
                        design = "crossed",
                        cross_vars = c("subject", "item"),
                        verbose = FALSE)
  expect_s3_class(result, "data.frame")
})

test_that("batch_eta2p returns one row per fixed effect", {
  result <- batch_eta2p(model_crossed, crossed_data,
                        design = "crossed",
                        cross_vars = c("subject", "item"),
                        verbose = FALSE)
  n_effects <- length(names(lme4::fixef(model_crossed))) - 1
  expect_equal(nrow(result), n_effects)
})

test_that("batch_eta2p contains required columns", {
  result <- batch_eta2p(model_crossed, crossed_data,
                        design = "crossed",
                        cross_vars = c("subject", "item"),
                        verbose = FALSE)
  expect_true(all(c("effect", "eta2p", "variance_effect",
                    "variance_error", "type") %in% colnames(result)))
})

test_that("batch_eta2p matches single eta2p call", {
  batch  <- batch_eta2p(model_crossed, crossed_data,
                        design = "crossed",
                        cross_vars = c("subject", "item"),
                        verbose = FALSE)
  single <- eta2p(model_crossed, "condition", crossed_data,
                  design = "crossed",
                  cross_vars = c("subject", "item"),
                  verbose = FALSE)
  expect_equal(batch$eta2p[batch$effect == "condition"],
               unname(single$eta2p), tolerance = 1e-10)
})


# VERBOSE OUTPUT


test_that("verbose = FALSE suppresses output", {
  expect_silent(
    eta2p(model_crossed, "condition", crossed_data,
          design = "crossed", cross_vars = c("subject", "item"),
          verbose = FALSE)
  )
})
