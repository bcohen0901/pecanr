# tests/testthat/test-double-bar.R
#
# Regression tests for the <= 0.3.0 bug in which random-effect components
# stored under suffixed group names ("subj.1", "subj.2", ... -- produced by
# the double-bar syntax (x || g) or by explicit (1|g) + (0+x|g) terms) were
# silently dropped from the error variance, inflating eta2p.
#
# Strategy: for a fitted model, reconstruct the expected error variance BY
# HAND from VarCorr() -- summing across ALL group rows including suffixed
# names -- and require eta2p()'s variance_error to match it exactly. This is
# a sharper test than comparing | vs || fits (those estimate different
# G structures, so their variance components legitimately differ slightly);
# an approximate |-vs-|| sanity check is included last.

library(testthat)
library(lme4)
library(pecanr)

# hand computation: residual + all intercept VCs + all slope VCs * Var(col),
# reading slope predictor columns from the model design matrix
hand_error_general <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  X  <- lme4::getME(model, "X")
  own <- vc[is.na(vc$var2), ]
  tot <- 0
  for (i in seq_len(nrow(own))) {
    if (own$grp[i] == "Residual") { tot <- tot + own$vcov[i]; next }
    v <- own$var1[i]
    if (is.na(v) || v == "(Intercept)") tot <- tot + own$vcov[i]
    else tot <- tot + own$vcov[i] * stats::var(X[, v])
  }
  tot
}

hand_error_within <- function(model) {   # residual + slopes only
  vc <- as.data.frame(VarCorr(model))
  X  <- lme4::getME(model, "X")
  own <- vc[is.na(vc$var2), ]
  tot <- 0
  for (i in seq_len(nrow(own))) {
    if (own$grp[i] == "Residual") { tot <- tot + own$vcov[i]; next }
    v <- own$var1[i]
    if (!is.na(v) && v != "(Intercept)")
      tot <- tot + own$vcov[i] * stats::var(X[, v])
  }
  tot
}

make_crossed_data <- function(seed = 1, ns = 40, nt = 12) {
  set.seed(seed)
  d <- data.frame(subj = factor(rep(seq_len(ns), each = nt)))
  d$x <- rnorm(ns * nt); d$z <- rnorm(ns * nt)
  u0 <- rnorm(ns); u1 <- rnorm(ns, 0, 0.4); u2 <- rnorm(ns, 0, 0.3)
  d$y <- 0.3 * d$x + 0.2 * d$z +
    u0[d$subj] + u1[d$subj] * d$x + u2[d$subj] * d$z + rnorm(ns * nt)
  d
}

test_that("double-bar single slope: slope variance enters the error", {
  d <- make_crossed_data()
  m <- suppressWarnings(lmer(y ~ x + (x || subj), data = d))
  vc <- as.data.frame(VarCorr(m))
  # precondition: this model really does store the slope under a suffixed name
  expect_true(any(grepl("^subj\\.[0-9]+$", vc$grp)))

  r_gen <- suppressWarnings(suppressMessages(
    eta2p(m, "x", d, design = "crossed", cross_vars = "subj")))
  expect_equal(r_gen$variance_error, hand_error_general(m), tolerance = 1e-8)
  expect_gt(sum(unlist(r_gen$variance_components$factor_slopes)), 0)

  r_op <- suppressWarnings(suppressMessages(
    eta2p(m, "x", d, design = "crossed", cross_vars = "subj",
          operative = TRUE)))
  expect_equal(r_op$variance_error, hand_error_within(m), tolerance = 1e-8)
})

test_that("explicit (1|g) + (0+x|g) form: identical handling", {
  d <- make_crossed_data(seed = 2)
  m <- suppressWarnings(lmer(y ~ x + (1 | subj) + (0 + x | subj), data = d))
  r <- suppressWarnings(suppressMessages(
    eta2p(m, "x", d, design = "crossed", cross_vars = "subj")))
  expect_equal(r$variance_error, hand_error_general(m), tolerance = 1e-8)
})

test_that("three uncorrelated slopes (x + z + x:z || subj): all included", {
  d <- make_crossed_data(seed = 3)
  m <- suppressWarnings(lmer(y ~ x * z + (x + z + x:z || subj), data = d))
  vc <- as.data.frame(VarCorr(m))
  expect_true(sum(grepl("^subj", vc$grp)) >= 4)   # subj, subj.1, subj.2, subj.3
  r <- suppressWarnings(suppressMessages(
    eta2p(m, "x", d, design = "crossed", cross_vars = "subj")))
  # hand computation includes ALL slope rows (x, z, and the x:z product);
  # pecanr centers the x:z product's constituents, and with x, z generated
  # centered the two agree to sampling-noise-free numerical precision on the
  # single-df slope columns; interaction slope enters via its own path, so
  # compare against the components list instead of the scalar for x:z.
  # inspect the per-factor slope vector directly (unlist() would prefix
  # names with the grouping-factor name, e.g. "subj.x")
  slopes <- r$variance_components$factor_slopes[["subj"]]
  expect_true(all(c("x", "z") %in% names(slopes)))
  expect_true(any(grepl(":", names(slopes))))     # interaction slope present
  expect_gt(sum(slopes), 0)
})

test_that("nested design with (x1 || stud): student slope enters the error", {
  set.seed(4)
  nc <- 20; np <- 8; no <- 6
  d <- data.frame(
    class = factor(rep(seq_len(nc), each = np * no)),
    stud  = factor(rep(seq_len(nc * np), each = no)),
    x1    = rnorm(nc * np * no))
  u_c <- rnorm(nc, 0, 0.7); u_s <- rnorm(nc * np, 0, 0.8)
  u_s1 <- rnorm(nc * np, 0, 0.3)
  d$y <- 0.35 * d$x1 + u_c[d$class] + u_s[d$stud] + u_s1[d$stud] * d$x1 +
    rnorm(nrow(d))
  m <- suppressWarnings(lmer(y ~ x1 + (x1 || stud) + (1 | class), data = d))
  vc <- as.data.frame(VarCorr(m))
  expect_true(any(grepl("^stud\\.[0-9]+$", vc$grp)))
  r <- suppressWarnings(suppressMessages(
    eta2p(m, "x1", d, design = "nested", nest_vars = c("stud", "class"))))
  expect_equal(r$variance_error, hand_error_general(m), tolerance = 1e-8)
  expect_gt(r$variance_components$level_slopes[["stud"]], 0)
})

test_that("mixed design with (x || participant): slope enters the error", {
  set.seed(5)
  ns <- 30; nm <- 8; nph <- 3
  d <- expand.grid(participant = factor(seq_len(ns)),
                   photo = factor(seq_len(nm * nph)))
  d$stimmodel <- factor(rep(seq_len(nm), each = nph)[as.integer(d$photo)])
  d$x <- rnorm(nrow(d))
  u_p <- rnorm(ns, 0, 0.9); u_p1 <- rnorm(ns, 0, 0.3)
  u_m <- rnorm(nm, 0, 0.6); u_ph <- rnorm(nm * nph, 0, 0.5)
  d$y <- 0.35 * d$x + u_p[d$participant] + u_p1[d$participant] * d$x +
    u_m[d$stimmodel] + u_ph[d$photo] + rnorm(nrow(d))
  m <- suppressWarnings(lmer(
    y ~ x + (x || participant) + (1 | stimmodel) + (1 | photo:stimmodel),
    data = d))
  vc <- as.data.frame(VarCorr(m))
  expect_true(any(grepl("^participant\\.[0-9]+$", vc$grp)))
  r <- suppressWarnings(suppressMessages(
    eta2p(m, "x", d, design = "mixed", cross_vars = "participant",
          nest_vars = c("photo", "stimmodel"))))
  expect_equal(r$variance_error, hand_error_general(m), tolerance = 1e-8)
  expect_gt(sum(unlist(
    r$variance_components$crossed$factor_slopes[["participant"]])), 0)
})

test_that("sanity: | and || fits give close (not identical) eta2p", {
  # different G structures -> components differ slightly; this catches only
  # gross errors like the dropped-slope bug (which produced ~2x differences)
  d <- make_crossed_data(seed = 6)
  m1 <- suppressWarnings(lmer(y ~ x + (x |  subj), data = d))
  m2 <- suppressWarnings(lmer(y ~ x + (x || subj), data = d))
  e1 <- suppressWarnings(suppressMessages(
    eta2p(m1, "x", d, design = "crossed", cross_vars = "subj")))$eta2p
  e2 <- suppressWarnings(suppressMessages(
    eta2p(m2, "x", d, design = "crossed", cross_vars = "subj")))$eta2p
  expect_lt(abs(e1 - e2) / e1, 0.15)
})

test_that("grouping names containing dots do not false-match", {
  # "subj.extra" must NOT be absorbed into "subj": the suffix rule requires
  # the post-dot text to be purely numeric
  d <- make_crossed_data(seed = 7)
  d$subj.extra <- factor(rep(seq_len(20), length.out = nrow(d)))
  u <- rnorm(20, 0, 0.5)
  d$y <- d$y + u[d$subj.extra]
  m <- suppressWarnings(lmer(y ~ x + (x || subj) + (1 | subj.extra), data = d))
  r <- suppressWarnings(suppressMessages(
    eta2p(m, "x", d, design = "crossed", cross_vars = c("subj", "subj.extra"))))
  # subj's intercept must not double-count subj.extra's
  vc <- as.data.frame(VarCorr(m))
  int_subj_only <- sum(vc$vcov[pecanr:::match_re_groups(vc$grp, "subj") &
                                 vc$var1 == "(Intercept)" & is.na(vc$var2)])
  expect_equal(unname(r$variance_components$factor_intercepts["subj"]),
               int_subj_only, tolerance = 1e-10)
  expect_equal(r$variance_error, hand_error_general(m), tolerance = 1e-8)
})
