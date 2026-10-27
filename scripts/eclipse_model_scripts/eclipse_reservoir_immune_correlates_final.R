# Reservoir and Immune Correlates: Final Analysis ####
#
# Goal: show that model-inferred individual parameters (t_a for reservoir size, N_p for immune
# avidity) are associated with independent experimental measurements of the quantities they
# represent -- data that were not used to fit the model.
#
# Every analysis below uses the same "draw and pool" strategy:
#   - Individual parameters are SAMPLED FROM THE CONDITIONAL DISTRIBUTION (10 draws per animal),
#     not the EBEs, which shrink toward the population mean and bias associations toward the
#     null (Lavielle & Ribba 2016, Pharm Res; Ayral et al. 2021, CPT:PSP).
#   - The analysis is run separately on each draw (exactly one row per animal, so no
#     pseudoreplication), and the per-draw estimates are pooled with Rubin's rules
#     (mice::pool.scalar, rule = "rubin1987"). The between-draw variance carries the uncertainty
#     in each animal's true parameter into the final SE and p-value. No latent-variable
#     (brms me()) model is used: it applies a second population prior on top of draws that
#     already reflect the NLME's prior, which shrinks the covariate further (roughly 6-14% for a
#     typical animal, up to ~45% for the least certain).
#
# Analyses:
#   1. Kendall's tau (with an Akritas-Theil-Sen / Theil-Sen line) for reservoir size (censored)
#      and immune breadth, for each study and for the combined data.
#   2. Reservoir detection model (the detection part of the hurdle model), separately by study:
#      P(detectable) ~ t_a by logistic regression. Detection rates differ enormously by study
#      (17/18, 12/34 and 3/44 animals detectable in Studies 1-3, all with the same detection limit
#      of 3), so studies are never pooled and a study is only fit when it has enough animals in
#      each class (see min_per_class).
#   2b. Reservoir Tobit (left-censored) regression of log10 SIV DNA on t_a, separately by study.
#   2c. Reservoir censored negative binomial / zero-inflated negative binomial regression of SIV DNA on t_a (count-scale
#      alternative to the Tobit), separately by study.
#   3. Immune breadth beta-binomial regression, separately by study. Breadth is the number of
#      positive subpools out of a fixed, study-specific panel (see below).
#
# Reservoir: Studies 1-3 (PBMC). Study 4 uses a different assay (intact provirus) and is not part
# of this analysis. Immune: Study 2 at ART interruption and Study 4 post-MVA, as in the paper
# figure. Figures are saved to presentation figures; nothing here writes to paper_figures.

library(tidyverse)
library(here)
library(NADA2)
library(mblm)
library(mice)
library(glmmTMB)
library(censReg)
library(patchwork)

RERUN <- FALSE  # TRUE: rerun the bootstrapped Kendall's tau pooling (slow). FALSE: load cached results.
N_BOOT <- 200   # bootstrap replicates per draw for the Kendall's tau line/tau SEs

out_dir <- here("results_eclipse", "reservoir_immune_correlates_final")
# fig_dir <- here("presentation figures")
fig_dir <- out_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# Data setup ####

colors <- c("Control" = "black", "TLR7" = "#027dbb", "Vac" = "#fcba12", "Ab" = "#fa61c2",
            "TLR7 + Vac" = "#7ab148", "TLR7 + Ab" = "#7f48b1", "TLR7 + Vac + Ab" = "#a96e7c")

reservoir_covar <- read_csv(here("data", "derived", "reservoir_correlates.csv"), show_col_types = FALSE) %>%
  mutate(Treatment = factor(Treatment, levels = names(colors)))
immune_covar <- read_csv(here("data", "derived", "immune_correlates.csv"), show_col_types = FALSE) %>%
  mutate(Treatment = factor(Treatment, levels = names(colors)))

read_conditional <- function(study) {
  read_csv(here("results_eclipse", paste0("study", study), "best_model", "IndividualParameters",
                "simulatedIndividualParameters.txt"), show_col_types = FALSE) %>% 
    mutate(betas = betas*1e-7, NPs = NPs*1e5)
}
study1_conditional <- read_conditional(1)
study2_conditional <- read_conditional(2)
study3_conditional <- read_conditional(3)
study4_conditional <- read_conditional(4)

all_conditional_param <- rbind(
  study1_conditional %>% mutate(Study = 1) %>% select(-tlr7),
  study2_conditional %>% mutate(Study = 2) %>% select(-c(intTLR7VAC, tlr7, vac)),
  study3_conditional %>% mutate(Study = 3) %>% select(-c(inttlr7ab, tlr7, ab)),
  study4_conditional %>% mutate(Study = 4) %>% select(-c(ab, tlr7, vac))
)

# EBEs (conditional modes) from Monolix, used only to DISPLAY one point per animal in the Kendall
# figures; every statistic below still uses the 10 conditional-distribution draws.
ebe <- map_dfr(1:4, function(study) {
  read_csv(here("results_eclipse", paste0("study", study), "best_model", "IndividualParameters",
                "estimatedIndividualParameters.txt"), show_col_types = FALSE) %>%
    select(id, ta_mode, NPs_mode) %>% mutate(Study = study, NPs_mode = NPs_mode*1e5)
})

# Reservoir: one row per animal per draw. Animal ids are reused across studies, so animals are
# always identified by (Study, id).
reservoir_pbmc <- left_join(
  all_conditional_param %>% filter(Study != 4),
  reservoir_covar %>% mutate(log_SIV_DNA = log10(SIV_DNA)),
  relationship = "many-to-many"
) %>%
  filter(!is.na(SIV_DNA), cell_type == "PBMC") %>%
  mutate(log_ta = log10(ta), detected = 1 - Censored)

# Immune breadth: number of positive ELISPOT subpools out of a FIXED, STUDY-SPECIFIC panel.
# Panel sizes come from the source papers: Study 2 (Borducchi et al. 2016) 13 Gag + 26 Pol +
# 22 Env = 61 subpools; Study 4 (Walker-Sperling et al. 2022) 13 Gag + 26 Pol + 49 Env = 88.
# prop_positive puts the two studies on a comparable 0-1 scale (identical ranks to the raw count
# within a study, so per-study Kendall's tau is unaffected).
immune_breadth <- left_join(
  rbind(study2_conditional %>% select(-intTLR7VAC, -tlr7, -vac) %>% mutate(Study = 2),
        study4_conditional %>% select(-ab, -tlr7, -vac) %>% mutate(Study = 4)),
  immune_covar %>% select(Study, id, Timing, Treatment, observation = breadth_Total),
  relationship = "many-to-many"
) %>%
  filter(!is.na(observation),
         (Study == 2 & Timing == "ART interruption") | (Study == 4 & Timing == "Post MVA vaccination")) %>%
  mutate(log_NPs = log10(NPs),
         N_trials = ifelse(Study == 2, 61, 88),
         prop_positive = observation / N_trials)

n_animals <- function(data) nrow(distinct(data, Study, id))

reservoir_sets <- list(
  "Study 1" = reservoir_pbmc %>% filter(Study == 1),
  "Study 2" = reservoir_pbmc %>% filter(Study == 2),
  "Study 3" = reservoir_pbmc %>% filter(Study == 3),
  "Combined (Studies 1-3)" = reservoir_pbmc
)
immune_sets <- list(
  "Study 2" = immune_breadth %>% filter(Study == 2),
  "Study 4" = immune_breadth %>% filter(Study == 4),
  "Combined (Studies 2 + 4)" = immune_breadth,
  "Combined, excl. TLR7 + Ab" = immune_breadth %>% filter(Treatment != "TLR7 + Ab")
)


# Pooling helpers ####

# Rubin's-rules pooling of per-draw estimates and their variances; drops draws where the fit
# failed. n and k give dfcom = n - k (animals minus estimated parameters).
pool_scalar_ok <- function(ests, vars, n, k) {
  ok <- is.finite(ests) & is.finite(vars)
  mice::pool.scalar(ests[ok], vars[ok], n = n, k = k, rule = "rubin1987")
}

pool_summary <- function(pooled) {
  with(pooled, tibble(
    estimate = qbar,
    se       = sqrt(t),
    df       = df,
    p_value  = 2 * pt(-abs(qbar / sqrt(t)), df),
    lower    = qbar - qt(0.975, df) * sqrt(t),
    upper    = qbar + qt(0.975, df) * sqrt(t)
  ))
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", as.character(signif(p, 2)))
p_label <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2)))


# 1. Kendall's tau, draw and pool ####
# Per draw: Kendall's tau plus a robust line from the same rank-based fit. Tau and line SEs:
#   - reservoir (left-censored at the detection limit): NADA2::cenken() returns tau and the
#     Akritas-Theil-Sen line but no SEs, and no analytic alternative exists for censored data,
#     so SEs are bootstrapped (N_BOOT resamples of animals per draw).
#   - immune (uncensored): tau has an analytic asymptotic SE (cor.test); the Theil-Sen line
#     (mblm) needs a bootstrap because its summary() reports the MAD of the pairwise slopes, not
#     the SE of the median estimator.

pool_kendall_censored <- function(data, x_var, y_var, cens_var, label, n_boot = N_BOOT, seed = 1) {
  set.seed(seed)
  n_draws <- max(data$rep)
  res <- matrix(NA_real_, n_draws, 7,
                dimnames = list(NULL, c("tau", "tau_var", "int", "int_var", "slp", "slp_var", "cov")))
  fit_fun <- function(dd) {
    NADA2::cenken(y = dd[[y_var]], ycen = as.logical(dd[[cens_var]]), x = dd[[x_var]])
  }

  for (i in seq_len(n_draws)) {
    cat(sprintf("[%s] draw %d/%d\n", label, i, n_draws))
    d <- data[data$rep == i, ]
    f0 <- tryCatch(fit_fun(d), error = function(e) NULL)
    if (is.null(f0)) next
    boots <- t(replicate(n_boot, {
      fb <- tryCatch(fit_fun(d[sample(nrow(d), replace = TRUE), ]), error = function(e) NULL)
      if (is.null(fb)) rep(NA_real_, 3) else c(fb$tau, fb$intercept, fb$slope)
    }))
    res[i, ] <- c(f0$tau, var(boots[, 1], na.rm = TRUE), f0$intercept, var(boots[, 2], na.rm = TRUE),
                  f0$slope, var(boots[, 3], na.rm = TRUE), cov(boots[, 2], boots[, 3], use = "complete.obs"))
  }

  n <- n_animals(data)
  list(
    tau = pool_scalar_ok(res[, "tau"], res[, "tau_var"], n, k = 1),
    intercept = pool_scalar_ok(res[, "int"], res[, "int_var"], n, k = 1),
    slope = pool_scalar_ok(res[, "slp"], res[, "slp_var"], n, k = 1),
    n_animals = n, n_ok = sum(is.finite(res[, "tau"])), raw = as_tibble(res)
  )
}

pool_kendall_uncensored <- function(data, x_var, y_var, label, n_boot = N_BOOT, seed = 1) {
  set.seed(seed)
  n_draws <- max(data$rep)
  res <- matrix(NA_real_, n_draws, 7,
                dimnames = list(NULL, c("tau", "tau_var", "int", "int_var", "slp", "slp_var", "cov")))
  line_fun <- function(dd) {
    # classic Theil-Sen (median of all pairwise slopes), the line that goes with Kendall's tau;
    # repeated = TRUE would give Siegel's repeated medians, which behave differently with many ties/zeros
    coef(mblm(as.formula(paste(y_var, "~", x_var)), data = dd, repeated = FALSE))
  }

  for (i in seq_len(n_draws)) {
    cat(sprintf("[%s] draw %d/%d\n", label, i, n_draws))
    d <- data[data$rep == i, ]
    kt <- suppressWarnings(cor.test(d[[x_var]], d[[y_var]], method = "kendall", exact = FALSE))
    co <- line_fun(d)
    boots <- t(replicate(n_boot, {
      cb <- tryCatch(line_fun(d[sample(nrow(d), replace = TRUE), ]), error = function(e) NULL)
      if (is.null(cb)) rep(NA_real_, 2) else unname(cb)
    }))
    res[i, ] <- c(unname(kt$estimate), (unname(kt$estimate) / unname(kt$statistic))^2,
                  unname(co[1]), var(boots[, 1], na.rm = TRUE), unname(co[2]), var(boots[, 2], na.rm = TRUE),
                  cov(boots[, 1], boots[, 2], use = "complete.obs"))
  }

  n <- n_animals(data)
  list(
    tau = pool_scalar_ok(res[, "tau"], res[, "tau_var"], n, k = 1),
    intercept = pool_scalar_ok(res[, "int"], res[, "int_var"], n, k = 1),
    slope = pool_scalar_ok(res[, "slp"], res[, "slp_var"], n, k = 1),
    n_animals = n, n_ok = sum(is.finite(res[, "tau"])), raw = as_tibble(res), line_method = "theil-sen"
  )
}

kendall_file <- here(out_dir, "kendall_pooled.RData")

if (RERUN) {
  reservoir_kendall <- imap(reservoir_sets, ~ pool_kendall_censored(
    .x, "log_ta", "log_SIV_DNA", "Censored", label = paste("reservoir", .y)))
  immune_kendall <- imap(immune_sets, ~ pool_kendall_uncensored(
    .x, "log_NPs", "prop_positive", label = paste("immune", .y)))
  save(reservoir_kendall, immune_kendall, file = kendall_file)
} else {
  load(kendall_file)
  # Sets that are missing from the cache (or immune sets fitted with a different line estimator) are computed now
  # without redoing the cached ones; the reservoir bootstraps are the slow part, so each new set is saved once done.
  new_reservoir <- setdiff(names(reservoir_sets), names(reservoir_kendall))
  if (length(new_reservoir) > 0) {
    reservoir_kendall[new_reservoir] <- imap(reservoir_sets[new_reservoir], ~ pool_kendall_censored(
      .x, "log_ta", "log_SIV_DNA", "Censored", label = paste("reservoir", .y)))
    save(reservoir_kendall, immune_kendall, file = kendall_file)
  }
  current <- map_lgl(immune_sets, ~ FALSE)
  for (nm in names(immune_sets)) current[nm] <- identical(immune_kendall[[nm]]$line_method, "theil-sen")
  new_immune <- names(immune_sets)[!current]
  if (length(new_immune) > 0) {
    immune_kendall[new_immune] <- imap(immune_sets[new_immune], ~ pool_kendall_uncensored(
      .x, "log_NPs", "prop_positive", label = paste("immune", .y)))
    save(reservoir_kendall, immune_kendall, file = kendall_file)
  }
}

kendall_table <- bind_rows(
  imap_dfr(reservoir_kendall, ~ pool_summary(.x$tau) %>%
             mutate(outcome = "Reservoir (SIV DNA, censored)", dataset = .y, n_animals = .x$n_animals, draws_ok = .x$n_ok)),
  imap_dfr(immune_kendall, ~ pool_summary(.x$tau) %>%
             mutate(outcome = "Immune breadth (proportion positive)", dataset = .y, n_animals = .x$n_animals, draws_ok = .x$n_ok))
) %>% mutate(analysis = "Kendall's tau", term = "tau", note = NA_character_, .before = 1)
print(kendall_table %>% select(outcome, dataset, estimate, lower, upper, p_value, n_animals, draws_ok))


# 2. Reservoir detection model (the detection part of the hurdle model), draw and pool, by study ####
# P(detectable) ~ log10(t_a), logistic regression fit separately to each draw and each study.
# Detection rates differ enormously between studies, so studies are never pooled. A logistic slope
# needs enough animals in BOTH classes; if a study has fewer than min_per_class detected or
# undetected animals (e.g. Study 1: 1 undetectable animal; Study 3: 3 detectable animals) the fit
# is skipped and reported with the reason instead of forced. Lower min_per_class to see what the
# fit does anyway, but expect very wide intervals.

pool_detection_draws <- function(data, x_var = "log_ta", min_per_class = 5) {
  animals <- distinct(data, Study, id, detected)
  n_det <- sum(animals$detected == 1); n_not <- sum(animals$detected == 0)
  if (min(n_det, n_not) < min_per_class) {
    return(list(skipped = sprintf("too few animals in one class (%d detected, %d not detected; need >= %d each)",
                                  n_det, n_not, min_per_class)))
  }
  n_draws <- max(data$rep)
  res <- matrix(NA_real_, n_draws, 4, dimnames = list(NULL, c("int", "int_var", "slp", "slp_var")))
  for (i in seq_len(n_draws)) {
    d <- data[data$rep == i, ]
    d$x <- d[[x_var]]
    fit <- tryCatch(suppressWarnings(glm(detected ~ x, family = binomial, data = d)), error = function(e) NULL)
    if (is.null(fit) || !fit$converged) next
    co <- summary(fit)$coefficients
    res[i, ] <- c(co["(Intercept)", "Estimate"], co["(Intercept)", "Std. Error"]^2,
                  co["x", "Estimate"], co["x", "Std. Error"]^2)
  }
  n <- nrow(animals)
  list(
    intercept = pool_scalar_ok(res[, "int"], res[, "int_var"], n, k = 2),
    slope = pool_scalar_ok(res[, "slp"], res[, "slp_var"], n, k = 2),
    n_animals = n, n_ok = sum(is.finite(res[, "slp"]))
  )
}

detection_studies <- reservoir_sets[c("Study 1", "Study 2", "Study 3")]
detection_fits <- map(detection_studies, pool_detection_draws)

# One results row per fit; skipped fits are kept with the reason in `note`.
result_row <- function(fit, analysis, outcome, dataset, term) {
  if (!is.null(fit$skipped)) {
    return(tibble(analysis = analysis, outcome = outcome, dataset = dataset, term = term, note = fit$skipped))
  }
  pool_summary(fit$slope) %>%
    mutate(analysis = analysis, outcome = outcome, dataset = dataset, term = term,
           n_animals = fit$n_animals, draws_ok = fit$n_ok, note = NA_character_, .before = 1)
}

detection_table <- imap_dfr(detection_fits, ~ result_row(.x, "Reservoir detection model", "P(detectable) vs log10 t_a", .y,
                                                          "log-odds slope per log10 t_a")) %>%
  mutate(odds_ratio = exp(estimate))
print(detection_table %>% select(dataset, estimate, lower, upper, odds_ratio, p_value, n_animals, draws_ok, note))


# 2b. Reservoir Tobit regression, draw and pool, by study ####
# Left-censored (Tobit) regression of log10 SIV DNA on log10(t_a) via censReg, fit separately to
# each draw and each study, then pooled with Rubin's rules. Studies are never pooled: separate
# per-study fits fit far better than one shared fit (LRT on draw 1, p ~ 0.001) and the residual SD
# differs a lot between studies. The slope is on the latent (uncensored) log10 SIV DNA scale.
# sigma is the latent residual SD, pooled on the log scale (censReg's logSigma) and back-transformed.

pool_tobit_draws <- function(data, censor_value = 3, adjust_study = FALSE) {
  n_draws <- max(data$rep)
  form <- if (adjust_study) log_SIV_DNA ~ log_ta + factor(Study) else log_SIV_DNA ~ log_ta   # study intercepts (CSV only)
  n_par <- NA_integer_
  res <- matrix(NA_real_, n_draws, 7,
                dimnames = list(NULL, c("int", "int_var", "slp", "slp_var", "lsig", "lsig_var", "cov")))
  for (i in seq_len(n_draws)) {
    d <- data[data$rep == i, ]
    fit <- tryCatch(suppressWarnings(censReg::censReg(form, left = log10(censor_value), data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$estimate
    if (any(!is.finite(co[, "Std. error"]))) next
    res[i, ] <- c(co["(Intercept)", "Estimate"], co["(Intercept)", "Std. error"]^2,
                  co["log_ta", "Estimate"], co["log_ta", "Std. error"]^2,
                  co["logSigma", "Estimate"], co["logSigma", "Std. error"]^2,
                  vcov(fit)["(Intercept)", "log_ta"])
    n_par <- length(coef(fit))
  }
  n <- n_animals(data)
  list(
    intercept = pool_scalar_ok(res[, "int"], res[, "int_var"], n, k = n_par),
    slope = pool_scalar_ok(res[, "slp"], res[, "slp_var"], n, k = n_par),
    log_sigma = pool_scalar_ok(res[, "lsig"], res[, "lsig_var"], n, k = n_par),
    n_animals = n, n_ok = sum(is.finite(res[, "slp"])), raw = as_tibble(res)
  )
}

tobit_fits <- map(reservoir_sets[c("Study 1", "Study 2", "Study 3", "Combined (Studies 1-3)")], pool_tobit_draws)

tobit_table <- imap_dfr(tobit_fits, ~ result_row(.x, "Reservoir Tobit regression", "log10 SIV DNA vs log10 t_a", .y,
                                                  "slope per log10 t_a")) %>%
  mutate(sigma = map_dbl(tobit_fits, ~ exp(.x$log_sigma$qbar)))
print(tobit_table %>% select(dataset, estimate, se, lower, upper, p_value, sigma, n_animals, draws_ok, note))


# 2c. Reservoir censored (zero-inflated) negative binomial regression, draw and pool, by study ####
# A count-scale alternative to the Gaussian Tobit model. SIV DNA is an integer count per million cells (the smallest
# quantified value is 4), and animals at or below the detection limit (value <= 3) are left-censored. Per animal:
#   quantified:  (1 - pi) * NB(y; mu, theta)
#   censored:     pi + (1 - pi) * P_NB(Y <= 3; mu, theta)
# with log(mu) = b0 + b1 * log10(t_a) and NB size theta (variance mu + mu^2 / theta). pi is the zero-inflation
# probability (a "structural" zero: no measurable reservoir). Three versions are fit by maximum likelihood:
#   nb      no zero-inflation (pi = 0): a censored NB regression
#   zinb    pi constant
#   zinb_x  logit(pi) = g0 + g1 * log10(t_a)
# Each is fit to every draw and pooled with Rubin's rules (dfcom = n - number of parameters).
# Caution: below-limit values are only known to be <= 3, so a structural zero cannot be told apart from a small NB
# count; pi is identified by the SHAPE of the quantified values and the form of the model, not by observed zeros. The
# zero-inflated versions are therefore only fit where each side of the limit holds at least min_per_class animals, and
# draws whose fit lies on a parameter boundary or has a singular Hessian are dropped and counted.
# Slopes: b1 is the log rate ratio of mean SIV DNA (among animals that are not structural zeros) per log10 t_a; g1 is the
# log odds ratio of being a structural zero per log10 t_a; the average slope is d log E[SIV DNA] / d log10 t_a averaged
# over the animals (includes both parts, so it is one number comparable across the three versions).

zinb_models <- c(nb = "Reservoir censored NB regression",
                 zinb = "Reservoir censored ZINB regression (pi constant)",
                 zinb_x = "Reservoir censored ZINB regression (pi depends on t_a)")
log_sum_exp2 <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))

zinb_cens_fit <- function(d, model, center, limit = 3) {
  x <- d$log_ta - center; cens <- d$Censored == 1; y <- ifelse(cens, 0, round(d$SIV_DNA))
  np <- c(nb = 3, zinb = 4, zinb_x = 5)[[model]]
  nll <- function(par) {
    mu <- exp(par[1] + par[2] * x); theta <- exp(par[3])
    if (model == "nb") { lpi <- rep(-Inf, length(x)); l1mpi <- rep(0, length(x)) } else {
      eta <- par[4] + (if (model == "zinb_x") par[5] * x else 0)
      lpi <- plogis(eta, log.p = TRUE); l1mpi <- plogis(-eta, log.p = TRUE)
    }
    ll_det <- l1mpi + dnbinom(y, size = theta, mu = mu, log = TRUE)
    ll_cen <- log_sum_exp2(lpi, l1mpi + pnbinom(limit, size = theta, mu = mu, log.p = TRUE))
    v <- -sum(ifelse(cens, ll_cen, ll_det))
    if (is.finite(v)) v else 1e10
  }
  b0 <- log(mean(y[!cens])); fc <- mean(cens)
  starts <- lapply(c(0, -1, 2), function(lt) c(b0, 0, lt, if (np >= 4) qlogis(min(max(fc / 2, 0.05), 0.9)), if (np == 5) 0))
  if (np >= 4) starts <- c(starts, lapply(starts[1:2], function(s0) { s0[4] <- -2; s0 }), lapply(starts[1:2], function(s0) { s0[4] <- 1; s0 }))
  best <- NULL
  for (s0 in starts) {
    o <- tryCatch(suppressWarnings(optim(s0, nll, method = "BFGS", control = list(maxit = 3000, reltol = 1e-12))), error = function(e) NULL)
    if (!is.null(o) && o$value < 1e9 && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(list(ok = FALSE))
  H <- tryCatch(suppressWarnings(numDeriv::hessian(nll, best$par)), error = function(e) NULL)
  V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL)
  boundary <- any(abs(best$par[3:np]) > 10)   # theta -> 0 or Inf, or pi -> 0 or 1
  list(ok = !is.null(V) && all(is.finite(V)) && all(diag(V) > 0) && !boundary,
       par = best$par, V = V, nll = best$value, np = np, x = x)
}

# d log E[Y] / dx averaged over the animals; E[Y | x] = (1 - pi(x)) * mu(x)
zinb_avg_slope <- function(par, x, model) {
  if (model == "zinb_x") mean(par[2] - par[5] * plogis(par[4] + par[5] * x)) else par[2]
}
zinb_log10_mean <- function(par, xc, model) {
  l1mpi <- switch(model, nb = 0, zinb = plogis(-par[4], log.p = TRUE), zinb_x = plogis(-(par[4] + par[5] * xc), log.p = TRUE))
  (par[1] + par[2] * xc + l1mpi) / log(10)
}

pool_zinb_draws <- function(data, model, min_per_class = 5, limit = 3) {
  animals <- distinct(data, Study, id, Censored)
  n_cen <- sum(animals$Censored == 1); n_det <- sum(animals$Censored == 0)
  if (model != "nb" && min(n_cen, n_det) < min_per_class) {
    return(list(skipped = sprintf("zero-inflation not identifiable: %d below detection, %d quantified (need >= %d each)",
                                  n_cen, n_det, min_per_class)))
  }
  center <- mean(data$log_ta)
  fits <- map(seq_len(max(data$rep)), ~ zinb_cens_fit(data[data$rep == .x, ], model, center, limit))
  ok <- map_lgl(fits, "ok")
  if (!any(ok)) return(list(skipped = "no draw gave a well-defined fit (parameter on a boundary or singular Hessian)"))
  fits <- fits[ok]; np <- fits[[1]]$np; n <- n_animals(data)
  pool1 <- function(q, v) pool_scalar_ok(q, v, n, k = np)
  avg <- map(fits, function(f) {
    g <- numDeriv::grad(function(par) zinb_avg_slope(par, f$x, model), f$par)
    list(est = zinb_avg_slope(f$par, f$x, model), var = drop(t(g) %*% f$V %*% g))
  })
  list(count_slope = pool1(map_dbl(fits, ~ .x$par[2]), map_dbl(fits, ~ .x$V[2, 2])),
       zi_slope = if (model == "zinb_x") pool1(map_dbl(fits, ~ .x$par[5]), map_dbl(fits, ~ .x$V[5, 5])),
       avg_slope = pool1(map_dbl(avg, "est"), map_dbl(avg, "var")),
       model = model, center = center, fits = fits, n_animals = n, n_ok = sum(ok), n_draws = length(ok),
       aic = median(map_dbl(fits, ~ 2 * .x$nll + 2 * .x$np)), theta = median(map_dbl(fits, ~ exp(.x$par[3]))),
       pi = if (model != "nb") median(map_dbl(fits, ~ mean(plogis(.x$par[4] + (if (model == "zinb_x") .x$par[5] * .x$x else 0))))))
}

# Pooled 95% band for log10 of the expected SIV DNA, by the delta method per draw and Rubin's rules across draws
zinb_band <- function(fit, xs) {
  map_dfr(xs, function(x) {
    xc <- x - fit$center
    q <- map_dbl(fit$fits, ~ zinb_log10_mean(.x$par, xc, fit$model))
    v <- map_dbl(fit$fits, function(f) {
      g <- numDeriv::grad(function(par) zinb_log10_mean(par, xc, fit$model), f$par)
      drop(t(g) %*% f$V %*% g)
    })
    pool_summary(pool_scalar_ok(q, v, fit$n_animals, k = fit$fits[[1]]$np))
  })
}

zinb_fits <- map(set_names(c("Study 1", "Study 2", "Study 3", "Combined (Studies 1-3)")),
                 function(nm) map(set_names(names(zinb_models)), ~ pool_zinb_draws(reservoir_sets[[nm]], .x)))

zinb_table <- imap_dfr(zinb_fits, function(by_model, nm) {
  imap_dfr(by_model, function(f, model) {
    if (!is.null(f$skipped)) {
      return(tibble(analysis = zinb_models[[model]], outcome = "SIV DNA (per million cells) vs log10 t_a", dataset = nm,
                    term = "not fit", note = f$skipped))
    }
    one <- function(pooled, term) pool_summary(pooled) %>% mutate(term = term, ratio = exp(estimate))
    bind_rows(one(f$count_slope, "count slope (log rate ratio per log10 t_a)"),
              if (!is.null(f$zi_slope)) one(f$zi_slope, "zero-inflation slope (log odds per log10 t_a)"),
              one(f$avg_slope, "average slope of log expected SIV DNA per log10 t_a")) %>%
      mutate(analysis = zinb_models[[model]], outcome = "SIV DNA (per million cells) vs log10 t_a", dataset = nm,
             n_animals = f$n_animals, draws_ok = f$n_ok,
             note = sprintf("%d of %d draws well-defined; median AIC %.1f; median theta %.2f%s", f$n_ok, f$n_draws, f$aic, f$theta,
                            if (!is.null(f$pi)) sprintf("; median mean pi %.2f", f$pi) else ""),
             .before = 1)
  })
})
print(zinb_table %>% select(dataset, analysis, term, estimate, lower, upper, ratio, p_value, n_animals, note) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))


# 3. Immune breadth beta-binomial, draw and pool, by study ####
# Frequentist glmmTMB beta-binomial, cbind(positive subpools, negative subpools) ~ log10(N_p),
# fit to each draw and pooled. The beta-binomial (rather than binomial) allows the between-animal
# overdispersion the data show. dfcom uses the number of estimated parameters (intercept, slope,
# dispersion). Draws whose fit didn't converge properly (missing SEs) are dropped and counted.

pool_betabinom_draws <- function(data, x_var = "log_NPs", y_var = "observation", n_var = "N_trials", adjust_study = FALSE) {
  n_draws <- max(data$rep)
  form <- if (adjust_study) cbind(succ, fail) ~ x + factor(Study) else cbind(succ, fail) ~ x   # study effect (CSV only)
  res <- matrix(NA_real_, n_draws, 5, dimnames = list(NULL, c("int", "int_var", "slp", "slp_var", "cov")))
  n_par <- NA_integer_
  for (i in seq_len(n_draws)) {
    d <- data[data$rep == i, ]
    d$succ <- d[[y_var]]; d$fail <- d[[n_var]] - d[[y_var]]; d$x <- d[[x_var]]
    fit <- tryCatch(suppressWarnings(glmmTMB(form, family = betabinomial(), data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$coefficients$cond
    if (any(is.na(co[, "Std. Error"]))) next
    res[i, ] <- c(co["(Intercept)", "Estimate"], co["(Intercept)", "Std. Error"]^2,
                  co["x", "Estimate"], co["x", "Std. Error"]^2, vcov(fit)$cond["(Intercept)", "x"])
    n_par <- length(fixef(fit)$cond) + 1   # + dispersion
  }
  n <- n_animals(data)
  list(
    intercept = pool_scalar_ok(res[, "int"], res[, "int_var"], n, k = n_par),
    slope = pool_scalar_ok(res[, "slp"], res[, "slp_var"], n, k = n_par),
    n_animals = n, n_ok = sum(is.finite(res[, "slp"])), raw = as_tibble(res)
  )
}

immune_betabinom <- map(immune_sets[c("Study 2", "Study 4", "Combined (Studies 2 + 4)", "Combined, excl. TLR7 + Ab")], pool_betabinom_draws)

betabinom_table <- imap_dfr(immune_betabinom, ~ result_row(.x, "Beta-binomial regression",
                                                            "Immune breadth (positive subpools / N_trials)", .y,
                                                            "log-odds slope per log10 N_p")) %>%
  mutate(odds_ratio = exp(estimate))
print(betabinom_table %>% select(dataset, estimate, lower, upper, odds_ratio, p_value, n_animals, draws_ok))


# 3b. Immune breadth zero-inflated negative binomial, draw and pool, pooled studies ####
# glmmTMB nbinom2 count model with a constant zero-inflation probability, fit to the RAW count of positive subpools:
#   y ~ log10(N_p), ziformula ~ 1, no panel-size offset and no study term.
# Counts from the two studies (panels of 61 and 88 subpools) are modelled as they are, so this is a pooled-studies
# comparison to the beta-binomial, which uses the panel size as its denominator. Per draw it gives the log rate ratio
# per log10 N_p of the mean count among non-structural zeros; the pooled band is the expected count (1 - pi) * mu,
# computed per draw with glmmTMB's delta-method SEs on the log scale and combined with Rubin's rules.
pool_zinb_immune_draws <- function(data, x_var = "log_NPs", y_var = "observation", n_grid = 60) {
  xs <- seq(min(data[[x_var]]), max(data[[x_var]]), length.out = n_grid)
  fits <- map(seq_len(max(data$rep)), function(i) {
    d <- data[data$rep == i, ]; d$y <- d[[y_var]]; d$x <- d[[x_var]]
    f <- tryCatch(suppressWarnings(glmmTMB(y ~ x, ziformula = ~1, family = nbinom2, data = d)), error = function(e) NULL)
    if (is.null(f) || !isTRUE(f$sdr$pdHess)) return(NULL)
    pr <- tryCatch(suppressWarnings(predict(f, newdata = data.frame(x = xs), type = "response", se.fit = TRUE)), error = function(e) NULL)
    if (is.null(pr)) return(NULL)
    co <- summary(f)$coefficients$cond
    list(b = co["x", "Estimate"], se = co["x", "Std. Error"], aic = AIC(f), npar = length(f$fit$par),
         pi = plogis(fixef(f)$zi[[1]]), theta = sigma(f), log_pred = log(pr$fit), var_log = (pr$se.fit / pr$fit)^2)
  })
  ok <- !map_lgl(fits, is.null)
  if (!any(ok)) return(list(skipped = "no draw gave a well-defined fit"))
  fits <- fits[ok]; n <- n_animals(data); k <- fits[[1]]$npar
  band <- map_dfr(seq_along(xs), ~ pool_summary(pool_scalar_ok(map_dbl(fits, function(f) f$log_pred[.x]),
                                                                 map_dbl(fits, function(f) f$var_log[.x]), n, k)))
  list(slope = pool_scalar_ok(map_dbl(fits, "b"), map_dbl(fits, ~ .x$se^2), n, k),
       grid = tibble(x = xs, estimate = exp(band$estimate), lower = exp(band$lower), upper = exp(band$upper)),
       n_animals = n, n_ok = sum(ok), n_draws = length(ok), aic = median(map_dbl(fits, "aic")),
       pi = median(map_dbl(fits, "pi")), theta = median(map_dbl(fits, "theta")))
}

immune_zinb <- map(immune_sets[c("Combined (Studies 2 + 4)", "Combined, excl. TLR7 + Ab")], pool_zinb_immune_draws)

immune_zinb_table <- imap_dfr(immune_zinb, function(f, nm) {
  if (!is.null(f$skipped)) {
    return(tibble(analysis = "Immune ZINB regression (no offset)", outcome = "Immune breadth (positive subpools, raw count)",
                  dataset = nm, term = "not fit", note = f$skipped))
  }
  pool_summary(f$slope) %>%
    mutate(ratio = exp(estimate), analysis = "Immune ZINB regression (no offset)",
           outcome = "Immune breadth (positive subpools, raw count)", dataset = nm,
           term = "count slope (log rate ratio per log10 N_p)", n_animals = f$n_animals, draws_ok = f$n_ok,
           note = sprintf("%d of %d draws well-defined; median AIC %.1f; median theta %.2f; median pi %.2f; no panel-size offset, no study term",
                          f$n_ok, f$n_draws, f$aic, f$theta, f$pi), .before = 1)
})
print(immune_zinb_table %>% select(dataset, term, estimate, lower, upper, ratio, p_value, n_animals, note) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))


# Summary table ####

final_results <- bind_rows(kendall_table, detection_table, tobit_table, zinb_table, betabinom_table, immune_zinb_table) %>%
  select(analysis, outcome, dataset, term, estimate, lower, upper, p_value, n_animals, draws_ok, note)
write_csv(final_results, here(out_dir, "final_results_summary.csv"))
print(final_results %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))


# Figures ####
# Kendall, Tobit and beta-binomial figures: one point per animal at its EBE (conditional mode) with
# a bar spanning the range of its 10 conditional-distribution draws; the line/curve and statistic
# come from the pooled draw-and-pool fit. Open circles = below detection limit (reservoir) or zero
# breadth (immune). Colored by Treatment. The detection figure still plots every draw.

# Jitter is drawn ONCE per animal so the point and its bar share the same y (a geom-level jitter
# would move them independently). Reservoir: animals below the detection limit are spread just
# under it. Immune: many animals sit at exactly 0, so those are spread just BELOW 0 (never above,
# so anything plotted above 0 is a true positive count); the fitted line and tau use the
# unjittered values.
set.seed(1)
reservoir_animals <- reservoir_pbmc %>%
  group_by(Study, id, Treatment, log_SIV_DNA, Censored) %>%
  summarise(x_min = min(log_ta), x_max = max(log_ta), .groups = "drop") %>%
  left_join(ebe %>% transmute(Study, id, x = log10(ta_mode)), by = c("Study", "id")) %>%
  mutate(status = ifelse(Censored == 1, "Below detection limit", "Detected"),
         y = ifelse(status == "Detected", log_SIV_DNA, log_SIV_DNA - 0.25 + runif(n(), -0.1, 0.1)))

immune_animals <- immune_breadth %>%
  group_by(Study, id, Treatment, prop_positive, observation, N_trials) %>%
  summarise(x_min = min(log_NPs), x_max = max(log_NPs), .groups = "drop") %>%
  left_join(ebe %>% transmute(Study, id, x = log10(NPs_mode)), by = c("Study", "id")) %>%
  mutate(status = ifelse(prop_positive == 0, "Zero breadth", "Positive breadth"),
         jit = runif(n()),
         y = ifelse(prop_positive == 0, -0.03 * jit, prop_positive),          # proportion scale (combined panels)
         y_count = ifelse(observation == 0, -0.8 * jit, observation))         # count scale (single-study panels)

# Single-study immune panels show the observed number of positive subpools out of that study's
# fixed panel (lines/curves are the proportion fit times N_trials, so tau and p are unchanged).
# Combined panels mix panels of 61 and 88 subpools, so they keep the proportion.
immune_view <- function(animals, nm) {
  if (nm %in% c("Study 2", "Study 4")) {
    n_trials <- unique(animals$N_trials)
    list(animals = mutate(animals, y = y_count), scale = n_trials,
         ylab = sprintf("Positive subpools (of %d)", n_trials))
  } else {
    list(animals = animals, scale = 1, ylab = "Proportion of subpools positive")
  }
}

animal_sets <- function(animals, sets) map(sets, ~ semi_join(animals, .x, by = c("Study", "id")))

# Invisible points (size and alpha 0) that give every treatment a legend entry with a real key,
# so panels with different treatments share one identical, mergeable legend.
legend_keys_layer <- function(treatments, x, y) {
  geom_point(data = tibble(x = x, y = y, Treatment = factor(treatments, levels = names(colors))),
             aes(x, y, color = Treatment), size = 0, alpha = 0, inherit.aes = FALSE, show.legend = TRUE)
}

point_size <- 3   # one point size for every per-animal panel and boxplot
reservoir_shapes <- c("Detected" = 16, "Below detection limit" = 1)
immune_shapes <- c("Positive breadth" = 16, "Zero breadth" = 1)

# Pooled 95% confidence band for a fitted line a + b*x: at each x the per-draw estimate a_i + b_i*x and
# its variance (var_a + 2x cov_ab + x^2 var_b) are combined with Rubin's rules, so the band carries both
# within-draw and between-draw uncertainty. k sets dfcom = n - k as in each analysis's pooling. Needs the
# per-draw intercept-slope covariance; caches written before it was stored return NULL (rerun to add it).
line_band <- function(fit, xs, k) {
  raw <- fit$raw
  if (is.null(raw) || !"cov" %in% names(raw)) return(NULL)
  r <- raw[is.finite(raw$int) & is.finite(raw$slp) & is.finite(raw$cov), ]
  map_dfr(xs, ~ pool_summary(pool_scalar_ok(r$int + r$slp * .x, r$int_var + 2 * .x * r$cov + .x^2 * r$slp_var,
                                             fit$n_animals, k)))
}

# Line/curve over the animals' x-range plus its confidence ribbon. fun maps the linear predictor to the
# plotted y-scale; nonneg clips count-scale lines (and the band) at 0.
line_layers <- function(fit, animals, k, fun = identity, nonneg = FALSE) {
  xs <- seq(min(animals$x_min), max(animals$x_max), length.out = 100)
  ln <- tibble(x = xs, y = fun(fit$intercept$qbar + fit$slope$qbar * xs))
  band <- line_band(fit, xs, k)
  rib <- if (!is.null(band)) tibble(x = xs, ymin = fun(band$lower), ymax = fun(band$upper))
  if (nonneg) {
    keep <- ln$y >= 0
    ln <- ln[keep, ]
    if (!is.null(rib)) rib <- mutate(rib[keep, ], ymin = pmax(ymin, 0))
  }
  list(ribbon = if (!is.null(rib)) geom_ribbon(data = rib, aes(x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
                                               fill = "grey40", alpha = 0.25),
       line = geom_line(data = ln, aes(x, y), inherit.aes = FALSE, linewidth = 1))
}

# Shared point style for every per-animal figure: EBE point, bar = range of the animal's draws,
# open circle for below-limit / zero animals. all_treatments (optional) adds invisible legend keys
# so panels holding different treatments can share one identical legend.
integer_breaks <- function(view = NULL) {
  function(lims) {
    r <- if (is.null(view)) lims else view
    by <- if (diff(r) > 6) 2 else 1
    seq(ceiling(r[1] / by) * by, floor(r[2] / by) * by, by = by)
  }
}

animal_base <- function(animals, shapes, limit_line = FALSE, all_treatments = NULL, underlay = NULL, ylim = NULL, ylog = FALSE) {
  p <- ggplot(animals, aes(x, y, color = Treatment))
  if (!is.null(underlay)) p <- p + underlay
  p <- p +
    geom_linerange(aes(xmin = x_min, xmax = x_max), orientation = "y", alpha = 0.4, show.legend = FALSE)
  if (limit_line) p <- p + geom_hline(yintercept = log10(3), lty = "dashed", color = "gray")
  p <- p +
    geom_point(aes(shape = status), size = point_size, alpha = 0.7) +
    scale_shape_manual(values = shapes, name = NULL)
  # keep the axis on the data (a confidence band can be much wider than the points)
  if (is.null(ylim)) { r <- range(animals$y); ylim <- r + c(-1, 1) * 0.08 * diff(r) }
  # limits + oob_keep put the requested range on the scale itself, so ticks are drawn across the whole
  # range even where no data or band reach (coord_cartesian alone would drop those ticks)
  p <- p + if (ylog) {
    scale_y_continuous(limits = ylim, oob = scales::oob_keep, breaks = integer_breaks(ylim),
                       labels = scales::label_math(expr = 10^.x))
  } else {
    scale_y_continuous(limits = ylim, oob = scales::oob_keep)
  }
  p <- p + coord_cartesian(ylim = ylim) +
    scale_x_continuous(breaks = integer_breaks(), labels = scales::label_math(expr = 10^.x))
  if (is.null(all_treatments)) return(p + scale_color_manual(values = colors))
  p <- p + 
    legend_keys_layer(all_treatments, mean(animals$x), mean(animals$y)) +
    scale_color_manual(values = colors, limits = all_treatments) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = point_size)))  
    return(p)
}

keep_legend_on_last_panel <- function(panels) {
  n <- length(panels)
  panels[-n] <- map(panels[-n], ~ .x + guides(color = "none", shape = "none"))
  panels
}

kendall_panel <- function(animals, fit, title, xlab, ylab, shapes, limit_line = FALSE, scale = 1,
                          all_treatments = NULL, nonneg_line = FALSE, ylim = NULL, ylog = FALSE, show_n = TRUE,
                          show_stats = TRUE) {
  tab <- pool_summary(fit$tau)
  # count outcomes can't be negative: nonneg_line clips the line and band at 0
  ll <- line_layers(fit, animals, k = 1, fun = function(eta) scale * eta, nonneg = nonneg_line)
  animal_base(animals, shapes, limit_line, all_treatments, underlay = ll$ribbon, ylim = ylim, ylog = ylog) +
    ll$line +
    (if (show_stats) geom_text(data = tibble(label = sprintf("tau = %.2f\n%s", tab$estimate, p_label(tab$p_value))),
              aes(Inf, Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 4)) +
    labs(x = xlab, y = ylab, color = "Treatment", title = paste0(title, if (show_n) paste0(" (n = ", fit$n_animals, ")"))) +
    theme_bw()
}

# With several panels holding different treatments, patchwork can't merge their Treatment legends.
# Only the panel showing the most treatments keeps the color legend (the combined panel).
last_panel_keeps_color_legend <- function(panels) {
  keep <- which.max(map_int(panels, ~ n_distinct(.x$data$Treatment)))
  panels[-keep] <- map(panels[-keep], ~ .x + guides(color = "none"))
  panels
}

reservoir_kendall_plot <- wrap_plots(last_panel_keeps_color_legend(imap(animal_sets(reservoir_animals, reservoir_sets), ~ kendall_panel(
  .x, reservoir_kendall[[.y]], .y,
  bquote("Individual"~t['a']~"estimates"), "SIV DNA (per million cells)",
  shapes = reservoir_shapes, limit_line = TRUE, ylog = T))), ncol = 2) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") &
  plot_annotation(caption = "Points represent EBEs and line ranges indicate range of 10 samples from the conditional distribution.")

immune_kendall_plot <- wrap_plots(last_panel_keeps_color_legend(imap(animal_sets(immune_animals, immune_sets), function(animals, nm) {
  v <- immune_view(animals, nm)
  kendall_panel(v$animals, immune_kendall[[nm]], nm,
                bquote("Individual"~N['p']~"estimates"), v$ylab,
                shapes = immune_shapes, scale = v$scale)
})), ncol = 2) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") &
  labs(y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT") &
  plot_annotation(caption = "Points represent EBEs and line ranges indicate range of 10 samples from the conditional distribution.")

ggsave(file.path(fig_dir, "SIV_final_reservoir_kendall.pdf"), reservoir_kendall_plot, width = 10, height = 10, bg = "white")
ggsave(file.path(fig_dir, "SIV_final_immune_kendall.pdf"), immune_kendall_plot, width = 10, height = 10, bg = "white")

# Reservoir detection: one panel per study whose model was estimable
detection_panel <- function(data, fit, title, show_stats = TRUE) {
  d <- data %>% mutate(x = log_ta)
  tab <- pool_summary(fit$slope)
  curve <- tibble(x = seq(min(d$x), max(d$x), length.out = 200)) %>%
    mutate(prob = plogis(fit$intercept$qbar + fit$slope$qbar * x))
  ggplot(d, aes(x, detected, color = Treatment)) +
    geom_point(size = 2, alpha = 0.5, position = position_jitter(height = 0.03, width = 0)) +
    legend_keys_layer(levels(droplevels(reservoir_pbmc$Treatment)), mean(d$x), 0) +
    geom_line(data = curve, aes(x, prob), inherit.aes = FALSE, linewidth = 1) +
    (if (show_stats) geom_text(data = tibble(label = sprintf("slope = %.2f\nOR = %.2f per log10\n%s", tab$estimate,
                                            exp(tab$estimate), p_label(tab$p_value))),
              aes(Inf, 0.5, label = label), inherit.aes = FALSE, hjust = 1.05, size = 4)) +
    scale_color_manual(values = colors, limits = levels(droplevels(reservoir_pbmc$Treatment))) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    labs(x = bquote("Individual"~t['a']~"estimate (log10, all draws)"), y = "P(detectable reservoir)",
         color = "Treatment", title = paste0(title, ": detection (n = ", fit$n_animals, ")")) +
    theme_bw()
}

detection_panels <- imap(compact(map(detection_fits, ~ if (is.null(.x$skipped)) .x)),
                         ~ detection_panel(detection_studies[[.y]], detection_fits[[.y]], .y))
if (length(detection_panels) > 0) {
  reservoir_detection_plot <- wrap_plots(detection_panels, nrow = 1) + 
    scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = integer_breaks()) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "SIV_final_reservoir_detection.pdf"), reservoir_detection_plot,
         width = 5.5 * length(detection_panels), height = 5.5, bg = "white")
}

# Immune beta-binomial: EBE points with the pooled probability curve
betabinom_panel <- function(nm, all_treatments = levels(droplevels(immune_breadth$Treatment)), label = nm, ylim = NULL,
                            model_in_title = TRUE, show_n = TRUE, show_stats = TRUE) {
  v <- immune_view(animal_sets(immune_animals, immune_sets)[[nm]], nm)
  animals <- v$animals
  fit <- immune_betabinom[[nm]]
  tab <- pool_summary(fit$slope)
  ll <- line_layers(fit, animals, k = 3, fun = function(eta) v$scale * plogis(eta))
  animal_base(animals, immune_shapes, all_treatments = all_treatments, underlay = ll$ribbon, ylim = ylim) +
    ll$line +
    # geom_text(data = tibble(label = sprintf("slope = %.2f\nOR = %.2f per log10\n%s", tab$estimate,
    #                                         exp(tab$estimate), p_label(tab$p_value))),
    #           aes(Inf, Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 4) +
    (if (show_stats) geom_text(data = tibble(label = sprintf("OR = %.2f\n%s", 
                                            exp(tab$estimate), p_label(tab$p_value))),
              aes(Inf, Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 4)) +
    labs(x = bquote("Individual"~N['p']~"estimates"), y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT",
         color = "Treatment", title = paste0(label, if (model_in_title) ": beta-binomial", if (show_n) paste0(" (n = ", fit$n_animals, ")"))) +
    theme_bw()
}

immune_betabinom_plot <- wrap_plots(keep_legend_on_last_panel(map(c("Study 2", "Study 4"), betabinom_panel)), nrow = 1) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") & coord_cartesian(ylim = c(0,30)) &
  plot_annotation(caption = "Points represent EBEs and line ranges indicate range of 10 samples from the conditional distribution.")
ggsave(file.path(fig_dir, "SIV_final_immune_betabinom.pdf"), immune_betabinom_plot, width = 10, height = 5.5, bg = "white")

# Reservoir Tobit: EBE points with the pooled latent-scale line per study
tobit_panel <- function(nm, all_treatments = levels(droplevels(reservoir_pbmc$Treatment)), label = nm, ylim = NULL,
                        model_in_title = TRUE, ylog = TRUE, show_n = TRUE, show_stats = TRUE) {
  animals <- animal_sets(reservoir_animals, reservoir_sets)[[nm]]
  fit <- tobit_fits[[nm]]
  tab <- pool_summary(fit$slope)
  ll <- line_layers(fit, animals, k = 3)
  animal_base(animals, reservoir_shapes, limit_line = TRUE,
              all_treatments = all_treatments, underlay = ll$ribbon, ylim = ylim, ylog = ylog) +
    ll$line +
    (if (show_stats) geom_text(data = tibble(label = sprintf("slope = %.2f\n%s", tab$estimate, p_label(tab$p_value))),
              aes(Inf, Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 4)) +
    labs(x = bquote("Individual"~t['a']~"estimate"), y = "SIV DNA (per million cells)",
         color = "Treatment", title = paste0(label, if (model_in_title) ": Tobit", if (show_n) paste0(" (n = ", fit$n_animals, ")"))) +
    theme_bw() 
}

reservoir_tobit_plot <- wrap_plots(keep_legend_on_last_panel(map(c("Study 1", "Study 2", "Study 3"), tobit_panel, ylog = TRUE)), nrow = 1) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") & 
  plot_annotation(caption = "Points represent EBEs and line ranges indicate range of 10 samples from the conditional distribution.")
ggsave(file.path(fig_dir, "SIV_final_reservoir_tobit.pdf"), reservoir_tobit_plot, width = 16, height = 5.5, bg = "white")

# Reservoir censored (ZI)NB: EBE points with the pooled expected SIV DNA and its pooled 95% band. Studies with enough
# animals on both sides of the limit show the ZINB with pi depending on t_a; the others show the censored NB.
zinb_panel <- function(nm, model, all_treatments = levels(droplevels(reservoir_pbmc$Treatment)), label = nm, ylim = NULL,
                       show_n = TRUE, show_stats = TRUE, model_in_title = TRUE) {
  animals <- animal_sets(reservoir_animals, reservoir_sets)[[nm]]
  fit <- zinb_fits[[nm]][[model]]
  xs <- seq(min(animals$x_min), max(animals$x_max), length.out = 60)
  band <- zinb_band(fit, xs)
  tab <- pool_summary(fit$avg_slope)
  animal_base(animals, reservoir_shapes, limit_line = TRUE, all_treatments = all_treatments, ylim = ylim, ylog = TRUE,
              underlay = geom_ribbon(data = tibble(x = xs, ymin = band$lower, ymax = band$upper), aes(x, ymin = ymin, ymax = ymax),
                                     inherit.aes = FALSE, fill = "grey40", alpha = 0.25)) +
    geom_line(data = tibble(x = xs, y = band$estimate), aes(x, y), inherit.aes = FALSE, linewidth = 1) +
    (if (show_stats) geom_text(data = tibble(label = sprintf("avg slope = %.2f\n%s", tab$estimate, p_label(tab$p_value))),
                               aes(Inf, Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 4)) +
    labs(x = bquote("Individual"~t['a']~"estimate"), y = "SIV DNA (per million cells)", color = "Treatment",
         title = paste0(label, if (model_in_title) paste0(": ", c(nb = "censored NB", zinb = "censored ZINB", zinb_x = "censored ZINB (pi ~ t_a)")[[model]]),
                        if (show_n) paste0(" (n = ", fit$n_animals, ")"))) +
    theme_bw()
}
zinb_figure_models <- map_chr(zinb_fits[c("Study 1", "Study 2", "Study 3")], ~ if (is.null(.x$zinb_x$skipped)) "zinb_x" else "nb")
reservoir_zinb_plot <- wrap_plots(keep_legend_on_last_panel(imap(zinb_figure_models, ~ zinb_panel(.y, .x))), nrow = 1) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") &
  plot_annotation(caption = "Points represent EBEs and line ranges indicate range of 10 samples from the conditional distribution.")
ggsave(file.path(fig_dir, "SIV_final_reservoir_zinb.pdf"), reservoir_zinb_plot, width = 16, height = 5.5, bg = "white")


# Paper figure ####
# 2 x 4 layout, boxplot first then the association panel for each study. Top row (reservoir SIV DNA): Study 1
# boxplot, Study 1 fit, Study 2 boxplot, Study 2 fit. Bottom row (immune breadth): Study 2 boxplot, Study 2 fit,
# Study 4 boxplot, Study 4 fit.
# Every panel in a row shares the same y-axis limits (the Kendall/fit panel and boxplot of a study, and
# across studies), so the scatter and the group comparison can be read against each other; the top of
# the range leaves room for the significance brackets.
#
# Three versions are saved. The first two use ONE kind of method across both rows; "mixed" is for comparison:
#   - nonparametric: pooled Kendall's tau with a robust line (reservoir: censored Kendall / Akritas-Theil-Sen;
#     immune: Kendall / Theil-Sen).
#   - parametric: regression fits (reservoir: Tobit; immune: beta-binomial).
#   - mixed: nonparametric reservoir (censored Kendall / Akritas-Theil-Sen), parametric immune (beta-binomial).
# Group comparisons (the boxplot brackets): every pair of treatment groups within a study, Holm-adjusted across
# the pairs in each panel.
#   - reservoir (left-censored at the detection limit): Peto-Peto test, i.e. survminer::pairwise_survdiff
#     (survdiff with rho = 1) on data reflected with NADA::flip(NADA::Cen(...)), the Helsel/NADA approach for
#     left-censored data.
#   - immune (uncensored counts): stats::pairwise.wilcox.test (Wilcoxon rank-sum, normal approximation with tie
#     correction).
#   - Permutation where the large-sample approximation is not trusted: when the smaller group of a pair has fewer
#     than min_group_n (20) animals, the p-value is the fraction of 10,000 random relabelings of the animals whose
#     test statistic is at least as large as the observed one (same statistic as the asymptotic test). Holm is applied
#     to the p-values actually used. Every group here has 5-15 animals, so all comparisons use permutation; the
#     asymptotic p-values are kept in the output for reference.
# Study 3 is not shown: 41 of 44 animals are below detection, so no group comparison has any information.

library(ggbeeswarm)
library(survival)

sig_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

pairwise_group_tests <- function(df, value_var, censored_var = NULL, min_group_n = 20, n_perm = 10000, seed = 1) {
  df <- df %>% mutate(Treatment = droplevels(Treatment), value = .data[[value_var]])
  if (!is.null(censored_var)) {
    df$flipd <- NADA::flip(NADA::Cen(df$value, as.logical(df[[censored_var]]), type = "left"))
  }
  # pairwise p-value matrix -> one row per pair (group1 = the earlier factor level)
  tidy_p <- function(m) as.data.frame(as.table(m)) %>% filter(!is.na(Freq)) %>%
    transmute(group1 = as.character(Var2), group2 = as.character(Var1), p_asymptotic = Freq)
  # The test itself: an (unadjusted) p-value for every pair. The same function supplies the asymptotic
  # p-values below and is reused as the test statistic inside the permutation test.
  pairwise_p <- function(x) {
    suppressWarnings(
      if (is.null(censored_var)) {
        pairwise.wilcox.test(x$value, x$Treatment, p.adjust.method = "none", exact = FALSE)$p.value
      } else {
        survminer::pairwise_survdiff(flipd ~ Treatment, data = x, p.adjust.method = "none", rho = 1)$p.value
      })
  }
  asymptotic <- tidy_p(pairwise_p(df))

  # Permutation p-value for one pair, identical for both outcomes: keep only the two groups, shuffle their
  # labels, rerun the test, and count how often the shuffled result is at least as extreme as the observed
  # one. The test's p-value is the statistic (smaller = more extreme), which orders results exactly as the
  # underlying chi-square / rank-sum statistic does.
  perm_p <- function(pr) {
    d <- df[df$Treatment %in% pr, ]
    d$Treatment <- droplevels(d$Treatment)
    stat <- function(x) { p <- pairwise_p(x)[1, 1]; if (is.nan(p)) 1 else p }
    obs <- stat(d)
    set.seed(seed)
    sims <- replicate(n_perm, { x <- d; x$Treatment <- sample(x$Treatment); stat(x) })
    (1 + sum(sims <= obs * (1 + 1e-8), na.rm = TRUE)) / (n_perm + 1)
  }

  n_by_group <- table(df$Treatment)
  out <- asymptotic %>%
    mutate(n1 = as.integer(n_by_group[group1]), n2 = as.integer(n_by_group[group2]),
           method = ifelse(pmin(n1, n2) < min_group_n, "permutation", "asymptotic"), p = p_asymptotic)
  for (i in which(out$method == "permutation")) out$p[i] <- perm_p(c(out$group1[i], out$group2[i]))
  out %>% mutate(p_adj = p.adjust(p, method = "holm"), stars = sig_stars(p_adj))
}

# Bracket heights: each significant pair gets a bracket `pad` above the highest point among the groups it spans;
# brackets whose spans overlap are stacked `step` apart (narrow spans first). The heights are in data units of a
# linear axis, and the row's y-limit is then set from the highest bracket, so axes stay only as tall as needed.
bracket_table <- function(sig, groups, group_max, pad, step) {
  sig <- sig %>%
    filter(p_adj < 0.05) %>%
    mutate(x1 = match(group1, groups), x2 = match(group2, groups), lo = pmin(x1, x2), hi = pmax(x1, x2),
           span = hi - lo) %>%
    arrange(span, lo)
  y <- numeric(nrow(sig))
  for (i in seq_len(nrow(sig))) {
    y[i] <- max(group_max[sig$lo[i]:sig$hi[i]]) + pad
    earlier_overlap <- which(seq_len(nrow(sig)) < i & sig$lo <= sig$hi[i] & sig$hi >= sig$lo[i])
    if (length(earlier_overlap) > 0) y[i] <- max(y[i], max(y[earlier_overlap]) + step)
  }
  mutate(sig, y = y, x_mid = (x1 + x2) / 2)
}

group_max_of <- function(animals, y_var) {
  animals <- mutate(animals, Treatment = droplevels(Treatment))
  vapply(levels(animals$Treatment), function(g) max(animals[[y_var]][animals$Treatment == g]), numeric(1))
}

sig_bracket_layers <- function(sig, groups, group_max, pad, step, tip) {
  sig <- bracket_table(sig, groups, group_max, pad, step)
  if (nrow(sig) == 0) return(list())
  sig <- mutate(sig, y_tip = y - tip)
  list(
    geom_segment(data = sig, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE),
    geom_segment(data = sig, aes(x = x1, xend = x1, y = y, yend = y_tip), inherit.aes = FALSE),
    geom_segment(data = sig, aes(x = x2, xend = x2, y = y, yend = y_tip), inherit.aes = FALSE),
    geom_text(data = sig, aes(x = x_mid, y = y, label = stars), vjust = 0.01, size = 3.5, inherit.aes = FALSE)
  )
}

# Group tests use the original measurements (one row per animal).
box_data_reservoir <- function(study) {
  reservoir_covar %>% filter(cell_type == "PBMC", Study == study, !is.na(SIV_DNA)) %>% mutate(Treatment = droplevels(Treatment))
}
box_data_immune <- function(study, timing) {
  immune_covar %>% filter(Study == study, Timing == timing, !is.na(breadth_Total)) %>% mutate(Treatment = droplevels(Treatment))
}
# The permutation tests take about a minute, so they are cached; the cache is tagged with the method so a change
# of method (or RERUN) recomputes them.
tests_file <- here(out_dir, "treatment_group_tests.RData")
tests_version <- "v2: pairwise_survdiff / pairwise.wilcox.test reused inside the permutation (label shuffle within pair); permutation when smaller group < 20; 10000 permutations; Holm"
cache_ok <- !RERUN && file.exists(tests_file)
if (cache_ok) {
  load(tests_file)
  cache_ok <- exists("tests_version_saved") && identical(tests_version_saved, tests_version)
}
if (!cache_ok) {
  treatment_tests <- list(
    reservoir_1 = pairwise_group_tests(box_data_reservoir(1), "SIV_DNA", "Censored"),
    reservoir_2 = pairwise_group_tests(box_data_reservoir(2), "SIV_DNA", "Censored"),
    immune_2    = pairwise_group_tests(box_data_immune(2, "ART interruption"), "breadth_Total"),
    immune_4    = pairwise_group_tests(box_data_immune(4, "Post MVA vaccination"), "breadth_Total"))
  tests_version_saved <- tests_version
  save(treatment_tests, tests_version_saved, file = tests_file)
}
print(bind_rows(treatment_tests, .id = "panel") %>%
        select(panel, group1, group2, n1, n2, method, p_asymptotic, p, p_adj, stars) %>%
        mutate(across(c(p_asymptotic, p, p_adj), ~ signif(.x, 3))))
write_csv(bind_rows(treatment_tests, .id = "panel"), here(out_dir, "treatment_group_tests.csv"))

paper_studies <- reservoir_animals$Study %in% 1:2
paper_treatments <- names(colors)[names(colors) %in% c(as.character(reservoir_animals$Treatment[paper_studies]),
                                                       as.character(immune_animals$Treatment))]

# Shared y-limits per row: lower edge covers the jittered below-limit / zero points, upper edge the highest
# bracket in the row plus headroom for its stars. pad = clearance above the data, step = spacing of stacked brackets.
box_specs <- list(
  reservoir_1 = list(animals = filter(reservoir_animals, Study == 1), y_var = "log_SIV_DNA"),
  reservoir_2 = list(animals = filter(reservoir_animals, Study == 2), y_var = "log_SIV_DNA"),
  immune_2    = list(animals = animal_sets(immune_animals, immune_sets)[["Study 2"]], y_var = "observation"),
  immune_4    = list(animals = animal_sets(immune_animals, immune_sets)[["Study 4"]], y_var = "observation"))
panel_top <- function(name, pad, step) {
  a <- mutate(box_specs[[name]]$animals, Treatment = droplevels(Treatment))
  gm <- group_max_of(a, box_specs[[name]]$y_var)
  max(c(gm, bracket_table(treatment_tests[[name]], levels(a$Treatment), gm, pad, step)$y))
}
res_pad <- 0.3; res_step <- 0.3; res_head <- 0.2
imm_pad <- 2;   imm_step <- 2;   imm_head <- 1.5
res_ylim <- c(min(0, min(reservoir_animals$y[paper_studies]) - 0.1),
              max(panel_top("reservoir_1", res_pad, res_step), panel_top("reservoir_2", res_pad, res_step)) + res_head)
imm_ylim <- c(min(immune_animals$y_count) - 0.3,
              max(panel_top("immune_2", imm_pad, imm_step), panel_top("immune_4", imm_pad, imm_step)) + imm_head)

# Without brackets the rows need no headroom: limits come from the data alone.
res_ylim_nosig <- c(res_ylim[1], max(reservoir_animals$log_SIV_DNA[paper_studies]) + 0.3)
imm_ylim_nosig <- c(imm_ylim[1], max(immune_animals$observation) + 2)

paper_box <- function(animals, y_var, y_plot_var, shapes, tests, pad, step, tip, ylim, ylab, title, limit_line = FALSE, ylog = FALSE,
                      show_sig = TRUE) {
  animals <- mutate(animals, Treatment = droplevels(Treatment))
  p <- ggplot(animals, aes(Treatment, .data[[y_var]]))
  if (limit_line) p <- p + geom_hline(yintercept = log10(3), lty = "dashed", color = "gray75")
  p <- p + if (ylog) {
    scale_y_continuous(limits = ylim, oob = scales::oob_keep, breaks = integer_breaks(ylim),
                       labels = scales::label_math(expr = 10^.x))
  } else {
    scale_y_continuous(limits = ylim, oob = scales::oob_keep)
  }
  p +
    geom_boxplot(aes(color = Treatment), outlier.shape = NA) +
    geom_beeswarm(aes(y = .data[[y_plot_var]], color = Treatment, shape = status),
                  size = point_size, alpha = 0.7, cex = 1.8, priority = "density") +
    (if (show_sig) sig_bracket_layers(tests, levels(animals$Treatment), group_max_of(animals, y_var), pad, step, tip)) +
    scale_shape_manual(values = shapes, name = NULL) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_x_discrete(labels = scales::label_wrap(10)) +   # horizontal labels, at most two lines
    guides(color = "none", shape = "none") +
    coord_cartesian(ylim = ylim) +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw()
}

res_ylab <- "SIV DNA (per million cells)"
reservoir_row <- function(study, desc, tests, approach, show_sig = TRUE) {
  nm <- paste("Study", study)
  ylim_use <- if (show_sig) res_ylim else res_ylim_nosig
  fit_panel <- if (approach == "parametric") {
    tobit_panel(nm, all_treatments = paper_treatments, label = paste0(nm, ": ", desc), ylim = ylim_use,
                model_in_title = FALSE, ylog = TRUE, show_n = FALSE, show_stats = show_sig) + ylab(res_ylab)
  } else {
    kendall_panel(animal_sets(reservoir_animals, reservoir_sets)[[nm]], reservoir_kendall[[nm]],
                  paste0(nm, ": ", desc), bquote("Individual"~t['a']~"estimate"),
                  res_ylab, shapes = reservoir_shapes, limit_line = TRUE, all_treatments = paper_treatments,
                  ylim = ylim_use, ylog = TRUE, show_n = FALSE, show_stats = show_sig)  + ylab(res_ylab)
  }
  list(
    fit_panel,
    paper_box(filter(reservoir_animals, Study == study), "log_SIV_DNA", "log_SIV_DNA", reservoir_shapes, tests,
              res_pad, res_step, tip = 0.08, ylim = ylim_use, ylab = res_ylab, title = paste0(nm, ": ", desc),
              limit_line = TRUE, ylog = TRUE, show_sig = show_sig) + ylab(res_ylab)) 
}

immune_ylab <- "Number of peptide subpools\ngiving IFN-g+ ELISPOT"
immune_row <- function(study, desc, tests, approach, show_sig = TRUE) {
  nm <- paste("Study", study)
  ylim_use <- if (show_sig) imm_ylim else imm_ylim_nosig
  animals <- animal_sets(immune_animals, immune_sets)[[nm]]
  v <- immune_view(animals, nm)
  fit_panel <- if (approach == "parametric") {
    betabinom_panel(nm, all_treatments = paper_treatments, label = paste0(nm, ": ", desc), ylim = ylim_use,
                    model_in_title = FALSE, show_n = FALSE, show_stats = show_sig) + ylab(immune_ylab)
  } else {
    kendall_panel(v$animals, immune_kendall[[nm]], paste0(nm, ": ", desc),
                  bquote("Individual"~N['p']~"estimate"), v$ylab,
                  shapes = immune_shapes, scale = v$scale, all_treatments = paper_treatments, nonneg_line = TRUE,
                  ylim = ylim_use, show_n = FALSE, show_stats = show_sig) + ylab(immune_ylab)
  }
  list(
    fit_panel,
    paper_box(animals, "observation", "observation", immune_shapes, tests, imm_pad, imm_step, tip = 0.5,
              ylim = ylim_use, ylab = v$ylab, title = paste0(nm, ": ", desc), show_sig = show_sig) + ylab(immune_ylab))
}

immune_shapes <- c("Positive breadth" = 16, "Zero breadth" = 16)
point_size <- 3
# approach = "nonparametric": Kendall's tau + robust line in both rows (reservoir: censored Kendall /
# Akritas-Theil-Sen; immune: Kendall / Theil-Sen). approach = "parametric": regression in both rows
# (reservoir: Tobit; immune: beta-binomial). approach = "mixed": nonparametric reservoir, parametric immune.
make_paper_figure <- function(approach = c("nonparametric", "parametric", "mixed"), show_sig = TRUE) {
  approach <- match.arg(approach)
  reservoir_approach <- if (approach == "parametric") "parametric" else "nonparametric"
  immune_approach <- if (approach == "nonparametric") "nonparametric" else "parametric"
  panels <- c(reservoir_row(1, "SIV-TLR7", treatment_tests$reservoir_1, reservoir_approach, show_sig),
              reservoir_row(2, "SIV-TLR7-Vac", treatment_tests$reservoir_2, reservoir_approach, show_sig),
              immune_row(2, "SIV-TLR7-Vac", treatment_tests$immune_2, immune_approach, show_sig),
              immune_row(4, "SHIV-TLR7-Vac-Ab", treatment_tests$immune_4, immune_approach, show_sig))
  # Each study is a pair: group-comparison boxplot first, then the association panel.
  panels <- panels[c(2, 1, 4, 3, 6, 5, 8, 7)]
  # One Treatment legend for the whole figure (from the first association panel, b); the association
  # panels keep their own shape keys (below-limit / zero breadth), which merge across panels that share them.
  fit_panels <- c(4, 6, 8)
  panels[fit_panels] <- map(panels[fit_panels], ~ .x + guides(color = "none"))
  # the Treatment and shape keys (taken from panel b) are a single row; both guides need nrow = 1, otherwise
  # patchwork moves the collected legend out of the guide area
  panels[[2]] <- panels[[2]] + guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = point_size)),
                                      shape = guide_legend(nrow = 1))
  # the immune panels' zero-breadth open circles are explained in the figure legend text, not in the key
  panels[c(6, 8)] <- map(panels[c(6, 8)], ~ .x + guides(shape = "none"))
  # y-axis titles only on the first panel of each row (the axes are matched across the row)
  panels[-c(1, 5)] <- map(panels[-c(1, 5)], ~ .x + theme(axis.title.y = element_blank()))
  # layout: reservoir row (A-D), shared legend (E), immune row (F-I)
  wrap_plots(c(panels[1:4], list(guide_area()), panels[5:8]), design = "ABCD\nEEEE\nFGHI",
             heights = c(1, 0.14, 1)) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom", plot.tag = element_text(face = 2),
          plot.title = element_text(size = 14, hjust = 0.5), axis.title = element_text(size = 12))
}

ggsave(file.path(fig_dir, "SIV_paper_figure_nonparametric.pdf"), make_paper_figure("nonparametric") & theme(aspect.ratio = 1),
       width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_parametric.pdf"), make_paper_figure("parametric") & theme(aspect.ratio = 1),
       width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_mixed.pdf"), make_paper_figure("mixed") & theme(aspect.ratio = 1),
       width = 12, height = 8, bg = "white")



# Three-panel paper figure ####
# Per row (reservoir; immune): the association panel with the studies combined, then one boxplot per study.
# Nothing is annotated for significance (show_stats / show_sig are off); the numbers behind the figure go to
# three_panel_figure_results.csv so they can be quoted in the caption.
#   - Reservoir row: the scatter combines Studies 1-3; the boxplots show Studies 1 and 2 only (Study 3 has no group
#     comparison to show: 41 of 44 animals are below detection). Study 3 has the highest t_a and almost no detectable
#     virus, so it adds a between-study contrast to the combined tau (within Study 3 tau is about 0).
#   - Immune row: Studies 2 and 4. The two studies use different subpool panels (61 vs 88), so the whole row is
#     plotted as the PROPORTION of subpools positive, which lets the scatter and both boxplots share one y-axis.
#   - The combined line/tau is the pooled Kendall / ATS (reservoir) or Theil-Sen (immune) fit to all animals of the
#     two studies together. It mixes between-study differences with within-study association; the per-study tau
#     values are in the CSV next to it.

study_desc <- c("1" = "SIV-TLR7", "2" = "SIV-TLR7-Vac", "4" = "SHIV-TLR7-Vac-Ab")
res_comb <- "Combined (Studies 1-3)"
imm_comb <- "Combined (Studies 2 + 4)"
imm_excl <- "Combined, excl. TLR7 + Ab"
immune_ylab_prop <- "Proportion of peptide subpools\ngiving IFN-g+ ELISPOT"

# y-limits from the data alone (no brackets), shared by the three panels of a row; treatments cover every study shown
three_res_ylim <- c(min(0, min(reservoir_animals$y) - 0.1), max(reservoir_animals$log_SIV_DNA) + 0.3)
three_imm_ylim_counts <- c(min(immune_animals$y_count) - 0.3, max(immune_animals$observation) + 2)
three_treatments <- names(colors)[names(colors) %in% c(as.character(reservoir_animals$Treatment), as.character(immune_animals$Treatment))]
three_imm_ylim <- c(min(immune_animals$y) - 0.01, max(immune_animals$prop_positive) + 0.03)

# reservoir_fit: "kendall" (censored Kendall / ATS line), "tobit" or "zinb" (censored ZINB with pi depending on t_a).
# immune_fit: "kendall" (Kendall / Theil-Sen line) or "betabinom". The parametric rows fit the combined studies with one
# curve and no study term, exactly as the combined Kendall line does; study-adjusted slopes are in the results CSV.
# exclude_tlr7ab = TRUE drops the TLR7 + Ab animals (Study 4) from the immune scatter and its fit; the boxplots still show
# every group.
# Immune ZINB panel: EBE points on the raw count scale (zeros jittered below 0) with the pooled expected count and band.
immune_zinb_panel <- function(nm, all_treatments, label, ylim = NULL, show_n = TRUE, show_stats = TRUE) {
  animals <- animal_sets(immune_animals, immune_sets)[[nm]] %>% mutate(y = y_count)
  fit <- immune_zinb[[nm]]
  tab <- pool_summary(fit$slope)
  animal_base(animals, immune_shapes, all_treatments = all_treatments, ylim = ylim,
              underlay = geom_ribbon(data = fit$grid, aes(x, ymin = lower, ymax = upper), inherit.aes = FALSE,
                                     fill = "grey40", alpha = 0.25)) +
    geom_line(data = fit$grid, aes(x, estimate), inherit.aes = FALSE, linewidth = 1) +
    (if (show_stats) geom_text(data = tibble(label = sprintf("RR = %.2f per log10\n%s", exp(tab$estimate), p_label(tab$p_value))),
                               aes(Inf, Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 4)) +
    labs(x = bquote("Individual"~N['p']~"estimate"), y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT", color = "Treatment",
         title = paste0(label, if (show_n) paste0(" (n = ", fit$n_animals, ")"))) +
    theme_bw()
}

make_three_panel_figure <- function(reservoir_fit = c("kendall", "tobit", "zinb"), immune_fit = c("kendall", "betabinom", "zinb"),
                                    exclude_tlr7ab = FALSE) {
  reservoir_fit <- match.arg(reservoir_fit); immune_fit <- match.arg(immune_fit)
  # the no-offset ZINB models raw counts, so that option plots the whole immune row (scatter and boxplots) as counts
  imm_counts <- immune_fit == "zinb"
  imm_y_var <- if (imm_counts) "observation" else "prop_positive"
  imm_ylim <- if (imm_counts) three_imm_ylim_counts else three_imm_ylim
  imm_ylab <- if (imm_counts) immune_ylab else immune_ylab_prop
  imm_set <- if (exclude_tlr7ab) imm_excl else imm_comb
  imm_label <- if (exclude_tlr7ab) "Studies 2 and 4\n(excl. TLR7 + Ab)" else "Studies 2 and 4"
  res_animals <- semi_join(reservoir_animals, reservoir_sets[[res_comb]], by = c("Study", "id"))
  imm_animals <- immune_view(semi_join(immune_animals, immune_sets[[imm_set]], by = c("Study", "id")), imm_set)

  res_scatter <- switch(reservoir_fit,
    kendall = kendall_panel(res_animals, reservoir_kendall[[res_comb]], "Studies 1-3",
                            bquote("Individual"~t['a']~"estimate"), res_ylab, shapes = reservoir_shapes,
                            limit_line = TRUE, all_treatments = three_treatments, ylim = three_res_ylim,
                            ylog = TRUE, show_n = FALSE, show_stats = FALSE),
    tobit = tobit_panel(res_comb, all_treatments = three_treatments, label = "Studies 1-3", ylim = three_res_ylim,
                        model_in_title = FALSE, ylog = TRUE, show_n = FALSE, show_stats = FALSE),
    zinb = zinb_panel(res_comb, "zinb_x", all_treatments = three_treatments, label = "Studies 1-3", ylim = three_res_ylim,
                      show_n = FALSE, show_stats = FALSE, model_in_title = FALSE)) + ylab(res_ylab)
  imm_scatter <- (switch(immune_fit,
    kendall = kendall_panel(imm_animals$animals, immune_kendall[[imm_set]], imm_label,
                            bquote("Individual"~N['p']~"estimate"), imm_ylab, shapes = immune_shapes,
                            scale = imm_animals$scale, all_treatments = three_treatments, nonneg_line = TRUE,
                            ylim = imm_ylim, show_n = FALSE, show_stats = FALSE),
    betabinom = betabinom_panel(imm_set, all_treatments = three_treatments, label = imm_label,
                                ylim = imm_ylim, model_in_title = FALSE, show_n = FALSE, show_stats = FALSE),
    zinb = immune_zinb_panel(imm_set, all_treatments = three_treatments, label = imm_label, ylim = imm_ylim,
                             show_n = FALSE, show_stats = FALSE))) +
    ylab(imm_ylab)

  res_box <- function(study) {
    paper_box(filter(reservoir_animals, Study == study), "log_SIV_DNA", "log_SIV_DNA", reservoir_shapes, NULL,
              pad = 0, step = 0, tip = 0, ylim = three_res_ylim, ylab = res_ylab,
              title = paste0("Study ", study, ": ", study_desc[[as.character(study)]]),
              limit_line = TRUE, ylog = TRUE, show_sig = FALSE)
  }
  imm_box <- function(study) {
    paper_box(filter(immune_animals, Study == study), imm_y_var, imm_y_var, immune_shapes, NULL,
              pad = 0, step = 0, tip = 0, ylim = imm_ylim, ylab = imm_ylab,
              title = paste0("Study ", study, ": ", study_desc[[as.character(study)]]), show_sig = FALSE)
  }

  panels <- list(res_scatter, res_box(1), res_box(2), imm_scatter, imm_box(2), imm_box(4))
  # one Treatment legend (from the first panel), a single row; the other panels suppress theirs
  panels[[1]] <- panels[[1]] + guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = point_size)),
                                      shape = guide_legend(nrow = 1))
  panels[[4]] <- panels[[4]] + guides(color = "none", shape = "none")
  # y-axis titles only on the first panel of each row
  panels[-c(1, 4)] <- map(panels[-c(1, 4)], ~ .x + theme(axis.title.y = element_blank()))
  wrap_plots(c(panels[1:3], list(guide_area()), panels[4:6]), design = "ABC\nDDD\nEFG", heights = c(1, 0.14, 1)) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom", plot.tag = element_text(face = 2),
          plot.title = element_text(size = 14, hjust = 0.5), axis.title = element_text(size = 12))
}

ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel.pdf"), make_three_panel_figure("kendall", "kendall") & theme(aspect.ratio = 1),
       width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_tobit.pdf"), make_three_panel_figure("tobit", "betabinom") & theme(aspect.ratio = 1),
       width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_zinb.pdf"), make_three_panel_figure("zinb", "betabinom") & theme(aspect.ratio = 1),
       width = 12, height = 8, bg = "white")
# the same three figures with the TLR7 + Ab animals removed from the immune scatter
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_exclTLR7Ab.pdf"),
       make_three_panel_figure("kendall", "kendall", exclude_tlr7ab = TRUE) & theme(aspect.ratio = 1), width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_tobit_exclTLR7Ab.pdf"),
       make_three_panel_figure("tobit", "betabinom", exclude_tlr7ab = TRUE) & theme(aspect.ratio = 1), width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_zinb_exclTLR7Ab.pdf"),
       make_three_panel_figure("zinb", "betabinom", exclude_tlr7ab = TRUE) & theme(aspect.ratio = 1), width = 12, height = 8, bg = "white")
# reservoir ZINB with the immune ZINB (raw counts, no offset, no study term) instead of the beta-binomial
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_zinb_zinb.pdf"),
       make_three_panel_figure("zinb", "zinb") & theme(aspect.ratio = 1), width = 12, height = 8, bg = "white")
ggsave(file.path(fig_dir, "SIV_paper_figure_three_panel_zinb_zinb_exclTLR7Ab.pdf"),
       make_three_panel_figure("zinb", "zinb", exclude_tlr7ab = TRUE) & theme(aspect.ratio = 1), width = 12, height = 8, bg = "white")

# Numbers for the caption: association statistics, per-group summaries and the pairwise tests behind the figure.
association_rows <- function(fits, outcome, datasets, slope_units) {
  map_dfr(datasets, function(nm) {
    f <- fits[[nm]]; tau <- pool_summary(f$tau); sl <- pool_summary(f$slope)
    bind_rows(
      tibble(term = "Kendall tau", estimate = tau$estimate, lower = tau$lower, upper = tau$upper, p_value = tau$p_value),
      tibble(term = paste0("line slope (", slope_units, ")"), estimate = sl$estimate, lower = sl$lower, upper = sl$upper,
             p_value = sl$p_value)) %>%
      mutate(section = "association", outcome = outcome, dataset = nm, n_animals = f$n_animals, .before = 1)
  })
}

# Study-stratified summary: inverse-variance average of the per-study taus (each already pooled over draws). It removes
# between-study offsets, which the combined fit above mixes in; the SEs are the bootstrap/asymptotic ones from the pooling.
stratified_tau_row <- function(fits, outcome, datasets, label) {
  t <- map_dfr(datasets, ~ pool_summary(fits[[.x]]$tau)); w <- 1 / t$se^2
  est <- sum(w * t$estimate) / sum(w); se <- 1 / sqrt(sum(w))
  tibble(section = "association", outcome = outcome, dataset = label,
         term = "Kendall tau (inverse-variance average of within-study taus)",
         n_animals = sum(map_dbl(datasets, ~ fits[[.x]]$n_animals)),
         estimate = est, lower = est - 1.96 * se, upper = est + 1.96 * se, p_value = 2 * pnorm(-abs(est / se)))
}

group_rows <- bind_rows(
  map_dfr(1:2, function(s) {
    box_data_reservoir(s) %>% left_join(ebe %>% select(Study, id, ta_mode), by = c("Study", "id")) %>%
      group_by(Treatment) %>%
      summarise(n_animals = n(), n_flagged = sum(Censored == 1),
                median = if (any(Censored == 0)) median(SIV_DNA[Censored == 0]) else NA_real_,
                q25 = if (any(Censored == 0)) quantile(SIV_DNA[Censored == 0], 0.25, names = FALSE) else NA_real_,
                q75 = if (any(Censored == 0)) quantile(SIV_DNA[Censored == 0], 0.75, names = FALSE) else NA_real_,
                median_param = median(ta_mode), .groups = "drop") %>%
      mutate(outcome = "Reservoir (SIV DNA, copies per million cells; median/IQR among detectable animals)",
             dataset = paste("Study", s), median_prop = NA_real_,
             note = "n_flagged = animals below detection limit; median_param = median EBE t_a")
  }),
  map2_dfr(c(2, 4), c("ART interruption", "Post MVA vaccination"), function(s, timing) {
    box_data_immune(s, timing) %>% left_join(ebe %>% select(Study, id, NPs_mode), by = c("Study", "id")) %>%
      mutate(n_trials = immune_breadth$N_trials[match(Study, immune_breadth$Study)]) %>%
      group_by(Treatment) %>%
      summarise(n_animals = n(), n_flagged = sum(breadth_Total == 0), median = median(breadth_Total),
                q25 = quantile(breadth_Total, 0.25, names = FALSE), q75 = quantile(breadth_Total, 0.75, names = FALSE),
                median_prop = median(breadth_Total / first(n_trials)), median_param = median(NPs_mode), .groups = "drop") %>%
      mutate(outcome = "Immune (positive subpools; median/IQR among all animals)", dataset = paste("Study", s),
             note = "n_flagged = animals with zero breadth; median_param = median EBE N_p")
  })
) %>% mutate(section = "group summary", term = "group summary", group1 = as.character(Treatment), .keep = "unused")

pairwise_rows <- bind_rows(treatment_tests[c("reservoir_1", "reservoir_2", "immune_2", "immune_4")], .id = "panel") %>%
  mutate(section = "pairwise test",
         outcome = ifelse(grepl("^reservoir", panel), "Reservoir (SIV DNA)", "Immune (positive subpools)"),
         dataset = paste("Study", sub(".*_", "", panel)), term = "pairwise group comparison",
         n_animals = n1 + n2, p_value = p,
         note = paste0("p_value = permutation p; p_adj = Holm-adjusted within study and outcome; asymptotic p = ",
                       signif(p_asymptotic, 2))) %>%
  select(section, outcome, dataset, term, group1, group2, n_animals, p_value, p_adj, method, note)

# Parametric fits behind the Tobit / ZINB / beta-binomial versions of the figure (combined and per study). The study-adjusted
# slopes add study intercepts (Tobit) or a study effect (beta-binomial) to the combined fit, the analogue of the stratified
# tau. A study-adjusted ZINB is not fit: the zero-inflation part is not identifiable within Studies 1 and 3.
tobit_adjusted <- pool_tobit_draws(reservoir_sets[[res_comb]], adjust_study = TRUE)
betabinom_adjusted <- pool_betabinom_draws(immune_sets[[imm_comb]], adjust_study = TRUE)
betabinom_adjusted_excl <- pool_betabinom_draws(immune_sets[[imm_excl]], adjust_study = TRUE)
parametric_rows <- bind_rows(
  bind_rows(tobit_table, zinb_table, betabinom_table, immune_zinb_table) %>%
    filter(dataset %in% c(res_comb, imm_comb, imm_excl, "Study 1", "Study 2", "Study 3", "Study 4")) %>%
    transmute(section = "association (parametric)", outcome, dataset, term, n_animals, estimate, lower, upper, p_value,
              method = analysis, note),
  pool_summary(tobit_adjusted$slope) %>%
    mutate(section = "association (parametric)", outcome = "log10 SIV DNA vs log10 t_a", dataset = "Studies 1-3, adjusted for study",
           term = "slope per log10 t_a, study intercepts", n_animals = tobit_adjusted$n_animals, method = "Reservoir Tobit regression",
           note = "latent log10 SIV DNA scale") %>% select(-se, -df),
  pool_summary(betabinom_adjusted$slope) %>%
    mutate(section = "association (parametric)", outcome = "Immune breadth (positive subpools / N_trials)",
           dataset = "Studies 2 and 4, adjusted for study", term = "log-odds slope per log10 N_p, study effect",
           n_animals = betabinom_adjusted$n_animals, method = "Beta-binomial regression", note = "odds ratio = exp(estimate)") %>%
    select(-se, -df),
  pool_summary(betabinom_adjusted_excl$slope) %>%
    mutate(section = "association (parametric)", outcome = "Immune breadth (positive subpools / N_trials)",
           dataset = "Studies 2 and 4 excl. TLR7 + Ab, adjusted for study", term = "log-odds slope per log10 N_p, study effect",
           n_animals = betabinom_adjusted_excl$n_animals, method = "Beta-binomial regression", note = "odds ratio = exp(estimate)") %>%
    select(-se, -df))

three_panel_results <- bind_rows(
  association_rows(reservoir_kendall, "Reservoir (SIV DNA, censored)", c(res_comb, "Study 1", "Study 2", "Study 3"), "log10 SIV DNA per log10 t_a"),
  association_rows(immune_kendall, "Immune (proportion of subpools positive)", c(imm_comb, imm_excl, "Study 2", "Study 4"), "proportion per log10 N_p"),
  stratified_tau_row(reservoir_kendall, "Reservoir (SIV DNA, censored)", c("Study 1", "Study 2", "Study 3"), "Studies 1-3, stratified"),
  stratified_tau_row(reservoir_kendall, "Reservoir (SIV DNA, censored)", c("Study 1", "Study 2"), "Studies 1 and 2, stratified"),
  stratified_tau_row(immune_kendall, "Immune (proportion of subpools positive)", c("Study 2", "Study 4"), "Studies 2 and 4, stratified"),
  parametric_rows, group_rows, pairwise_rows) %>%
  select(any_of(c("section", "outcome", "dataset", "term", "group1", "group2", "n_animals", "n_flagged", "median", "q25", "q75",
                  "median_prop", "median_param", "estimate", "lower", "upper", "p_value", "p_adj", "method", "note")))
write_csv(three_panel_results, file.path(fig_dir, "three_panel_figure_results.csv"), na = "")
message("Wrote ", file.path(fig_dir, "three_panel_figure_results.csv"), " (", nrow(three_panel_results), " rows)")
