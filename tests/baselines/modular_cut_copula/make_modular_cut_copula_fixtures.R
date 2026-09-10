#!/usr/bin/env Rscript
##
## Regenerates the frozen external-oracle fixtures in this directory.
##
## WHY THESE FIXTURES ARE DIFFERENT IN KIND. PsyRAT's accuracy suite is
## otherwise regression-only: its oracle duplicates the production formulas, so a
## green run proves "unchanged", not "correct". The CSVs this script writes come
## from the owner's INDEPENDENT R + CmdStan reference implementation
## (person_specific_dynamic_concurrent_modular), so a MATLAB test that matches
## them is externally grounded rather than self-referential.
##
## Usage (from anywhere):
##
##   Rscript --vanilla make_modular_cut_copula_fixtures.R [output_dir]
##
## Point MPSDC_REFERENCE_DIR at the reference bundle if it is not at the default
## path below. The script is fully deterministic: rerunning it on an unchanged
## reference bundle reproduces every file byte for byte, so an unexpected diff
## means the reference changed and the MATLAB tolerances should be re-examined
## rather than the fixtures silently re-baselined.
##
## Regenerate ONLY from the reference bundle. Never hand-edit a fixture, and
## never regenerate them from the MATLAB implementation they exist to check -
## that would convert an external oracle into a regression snapshot.

options(stringsAsFactors = FALSE)

## Pin the RNG kinds so regeneration is byte-for-byte stable across R releases.
## R 3.6.0 changed sample()'s default; every release since uses these kinds, so
## this call is a no-op today and a guard tomorrow. It must precede every
## set.seed() below.
RNGversion("3.6.0")

code_dir <- Sys.getenv(
  "MPSDC_REFERENCE_DIR",
  "~/Desktop/stan_chisquare/person_specific_dynamic_concurrent_modular"
)
code_dir <- path.expand(code_dir)
if (!file.exists(file.path(code_dir, "R", "copula_stage.R"))) {
  stop("Reference bundle not found at: ", code_dir,
       "\nSet MPSDC_REFERENCE_DIR to its location.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else {
  file_arg <- sub("^--file=", "",
                  commandArgs(trailingOnly = FALSE)[
                    grepl("^--file=", commandArgs(trailingOnly = FALSE))])
  if (length(file_arg)) dirname(normalizePath(file_arg[[1]])) else getwd()
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c("common.R", "model_generation.R", "copula_stage.R", "reliability.R")) {
  source(file.path(code_dir, "R", f))
}
source(file.path(code_dir, "tests", "simulation_helpers.R"))

wcsv <- function(d, name) {
  ## A binary connection with an explicit eol pins the line endings to LF on
  ## every platform; a plain path would go through a text-mode connection and
  ## produce CRLF fixtures on Windows, which the frozen git bytes are not.
  con <- file(file.path(out_dir, name), open = "wb")
  on.exit(close(con), add = TRUE)
  utils::write.csv(d, con, row.names = FALSE, na = "", eol = "\n")
  cat("wrote", name, "-", nrow(d), "rows\n")
}

## ---------------------------------------------------------------------------
## 1. Gauss-Hermite copula -> standardized observed-scale covariance.
##    Isolates psyrat_gamma_copula_rescov across small/large nu and both signs
##    of rho, including the rho = 0 short circuit.
## ---------------------------------------------------------------------------
grid1 <- expand.grid(
  nu_theta = c(2.5, 18, 140),
  nu_alpha = c(4, 30, 90),
  rho      = c(-0.85, -0.30, 0, 0.087, 0.55, 0.93),
  KEEP.OUT.ATTRS = FALSE
)
wcsv(do.call(rbind, lapply(seq_len(nrow(grid1)), function(i) {
  r <- grid1[i, ]
  q <- copula_standard_covariance(r$nu_theta, r$nu_alpha, r$rho)
  data.frame(nu_theta = r$nu_theta, nu_alpha = r$nu_alpha, rho = r$rho,
             value = q$value, nodes = q$nodes, error = q$error,
             bound = sqrt((2 / r$nu_theta) * (2 / r$nu_alpha)))
})), "modular_cut_copula_quadrature_expected.csv")

## ---------------------------------------------------------------------------
## 2. Conditional rho grid posterior from sufficient statistics.
##    The uniform draw is supplied explicitly, which is what makes the
##    inverse-CDF draw comparable at all: MATLAB's random stream is not R's, so
##    only a SUPPLIED uniform can be matched. Cases 6 and 7 pile the posterior
##    against the +/-0.95 bound - the only regime where the trapezoid end
##    weights carry real mass, and the only one that drives rho_boundary_mass
##    into the range that raises the numerical-invalidity flag.
## ---------------------------------------------------------------------------
suff_cases <- list(
  list(N = 2000,  szt = 1996.2,   sza = 2011.7,  cross =  520.3, u = 0.5),
  list(N = 2000,  szt = 1996.2,   sza = 2011.7,  cross =  520.3, u = 0.113),
  list(N = 49858, szt = 49820.0,  sza = 49903.0, cross = 4322.0, u = 0.5),
  list(N = 500,   szt = 498.0,    sza = 502.0,   cross =  -95.0, u = 0.87),
  list(N = 5000,  szt = 4998.0,   sza = 5002.0,  cross =    3.0, u = 0.5),
  list(N = 500,   szt = 500.0,    sza = 500.0,   cross =  492.0, u = 0.5),
  list(N = 800,   szt = 800.0,    sza = 800.0,   cross = -790.0, u = 0.31)
)
wcsv(do.call(rbind, lapply(seq_along(suff_cases), function(i) {
  s <- suff_cases[[i]]
  suff <- c(N = s$N, sum_z_theta_sq = s$szt, sum_z_alpha_sq = s$sza,
            sum_cross = s$cross)
  post <- conditional_rho_posterior(suff, s$u, 20001L, 0.5)
  data.frame(case = i, N = s$N, sum_z_theta_sq = s$szt,
             sum_z_alpha_sq = s$sza, sum_cross = s$cross, uniform_draw = s$u,
             loglik_at_0_10 = gaussian_copula_log_likelihood(0.10, suff),
             loglik_at_neg_0_40 = gaussian_copula_log_likelihood(-0.40, suff),
             as.list(post))
})), "modular_cut_copula_rhogrid_expected.csv")

## ---------------------------------------------------------------------------
## 3. Stable Gamma PIT, with planted extremes so BOTH clip branches fire.
## ---------------------------------------------------------------------------
set.seed(20260809)
n3 <- 400L
mu3 <- exp(rnorm(n3, log(5.5), 0.30))
sigma3 <- exp(rnorm(n3, log(1.8333), 0.25))
y3 <- qgamma(pnorm(rnorm(n3)), shape = (mu3 / sigma3)^2, rate = mu3 / sigma3^2)
y3[1:6] <- c(1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2)   # far lower tail
y3[7:12] <- mu3[7:12] + 60 * sigma3[7:12]           # far upper tail
pit3 <- stable_gamma_pit(y3, mu3, sigma3, 1e-12)
wcsv(data.frame(y = y3, mu = mu3, sigma = sigma3, z = pit3$z),
     "modular_cut_copula_pit_expected.csv")
wcsv(data.frame(lower_clipped = pit3$lower_clipped,
                upper_clipped = pit3$upper_clipped, clip = 1e-12, n = n3),
     "modular_cut_copula_pit_clipcounts.csv")

## ---------------------------------------------------------------------------
## 4/5. FULL stage-2 run, one-facet and two-facet.
##    Random effects are deliberately NON-ZERO and distinct at every level, so
##    an indexing error in the reconstruction cannot pass by symmetry. The
##    two-facet run also dumps its draw-1 margins, which pin the reconstruction
##    to machine precision independently of the incomplete-gamma routines.
## ---------------------------------------------------------------------------
settings <- mpsdc_settings()
settings$copula_cores <- 1L
settings$reliability_cores <- 1L
settings$max_copula_draws <- 3L
settings$copula_grid_points <- 20001L
settings$seed <- 4242L

emit_stage2 <- function(facet, dimension, P, I, O, tag, dump_margins = FALSE) {
  st <- settings
  st$facet <- facet
  st$dimension <- dimension
  tmp <- tempfile(paste0("mpsdc_", tag, "_")); dir.create(tmp)
  path <- file.path(tmp, "d.rds")
  saveRDS(simulate_person_scale_concurrent_data(facet, P = P, I = I, O = O,
                                                seed = 31337L), path)
  prep <- prepare_modular_person_specific_data(path, st)
  sd_ <- prep$standata
  K <- sd_$K

  set.seed(90210)
  nd <- 3L
  draw <- data.frame(.chain = 1L, .iteration = seq_len(nd), .draw = seq_len(nd))
  for (m in 1:2) {
    draw[[sprintf("alpha_mu[%d]", m)]] <- c(0.05, -0.11, 0.02) * m
    draw[[sprintf("alpha_sigma[%d]", m)]] <- c(-0.07, 0.13, 0.04) * m
    for (k in seq_len(K)) {
      draw[[sprintf("beta_mu[%d,%d]", m, k)]] <- rnorm(nd, 0, 0.10)
      draw[[sprintf("beta_sigma[%d,%d]", m, k)]] <- rnorm(nd, 0, 0.08)
    }
  }
  for (p in seq_len(sd_$P)) for (m in 1:4) {
    draw[[sprintf("u_p_joint[%d,%d]", p, m)]] <- rnorm(nd, 0, 0.25)
  }
  for (i in seq_len(sd_$I)) for (m in 1:2) {
    draw[[sprintf("u_i[%d,%d]", i, m)]] <- rnorm(nd, 0, 0.10)
  }
  if (facet == "trial_occasion") {
    for (nm in c("o", "pi", "po", "io")) {
      lv <- sd_[[toupper(nm)]]
      for (l in seq_len(lv)) for (m in 1:2) {
        draw[[sprintf("u_%s[%d,%d]", nm, l, m)]] <- rnorm(nd, 0, 0.07)
      }
    }
  }

  stage <- run_modular_copula_stage(draw, sd_, st)

  ## Long format so MATLAB reassembles without a column-order convention.
  long <- do.call(rbind, lapply(
    setdiff(names(draw), c(".chain", ".iteration", ".draw")), function(nm) {
      base <- sub("\\[.*$", "", nm)
      idx <- gsub("\\[|\\]", "", regmatches(nm, regexpr("\\[.*\\]", nm)))
      parts <- as.integer(strsplit(idx, ",")[[1]])
      data.frame(draw = seq_len(nd), param = base, idx1 = parts[1],
                 idx2 = if (length(parts) > 1) parts[2] else NA_integer_,
                 value = draw[[nm]])
    }))
  wcsv(long, sprintf("modular_cut_copula_%s_draws.csv", tag))

  obs <- data.frame(y_theta = sd_$y_theta, y_alpha = sd_$y_alpha,
                    p_id = sd_$p_id, i_id = sd_$i_id)
  for (k in seq_len(K)) obs[[sprintf("Z%d", k)]] <- sd_$Z[, k]
  if (facet == "trial_occasion") {
    obs$o_id <- sd_$o_id; obs$pi_id <- sd_$pi_id
    obs$po_id <- sd_$po_id; obs$io_id <- sd_$io_id
  }
  wcsv(obs, sprintf("modular_cut_copula_%s_standata.csv", tag))

  wcsv(data.frame(
    facet = facet, dimension = dimension, N = sd_$N, K = K, P = sd_$P, I = sd_$I,
    O = if (is.null(sd_$O)) NA_integer_ else sd_$O,
    PI = if (is.null(sd_$PI)) NA_integer_ else sd_$PI,
    PO = if (is.null(sd_$PO)) NA_integer_ else sd_$PO,
    IO = if (is.null(sd_$IO)) NA_integer_ else sd_$IO,
    mu_offset_1 = sd_$mu_offset[1], mu_offset_2 = sd_$mu_offset[2],
    log_sigma_offset_1 = sd_$log_sigma_offset[1],
    log_sigma_offset_2 = sd_$log_sigma_offset[2],
    n_draws = nd, grid_points = st$copula_grid_points,
    prior_sd = st$copula_prior_sd, pit_clip = st$pit_clip
  ), sprintf("modular_cut_copula_%s_meta.csv", tag))

  wcsv(stage, sprintf("modular_cut_copula_%s_stage2_expected.csv", tag))

  if (dump_margins) {
    mar <- conditional_margin_parameters(draw[1, , drop = FALSE], sd_, facet)
    p1 <- stable_gamma_pit(sd_$y_theta, mar$mu_theta, mar$sigma_theta, 1e-12)
    p2 <- stable_gamma_pit(sd_$y_alpha, mar$mu_alpha, mar$sigma_alpha, 1e-12)
    wcsv(data.frame(mu_1 = mar$mu_theta, sigma_1 = mar$sigma_theta, z_1 = p1$z,
                    mu_2 = mar$mu_alpha, sigma_2 = mar$sigma_alpha, z_2 = p2$z,
                    shape_2 = (mar$mu_alpha / mar$sigma_alpha)^2),
         sprintf("modular_cut_copula_%s_draw1_margins.csv", tag))
  }
  unlink(tmp, recursive = TRUE)
  ## Returned so the stage-3 section below can reuse this exact draw set and
  ## standata. Nothing above is recomputed there, so the stage-2 fixtures stay
  ## byte-identical to what this script produced before stage 3 existed.
  invisible(list(draw = draw, stage = stage, standata = sd_,
                 metadata = prep$metadata, settings = st))
}

onefacet <- emit_stage2("trial", "one_dimension", 6L, 5L, 1L, "onefacet")
twofacet <- emit_stage2("trial_occasion", "two_dimensions", 5L, 4L, 2L, "twofacet",
                        dump_margins = TRUE)

## ---------------------------------------------------------------------------
## 6. FULL stage-3 run, TWO-FACET only.
##    This is the externally grounded oracle for
##    psyrat_ssrel_diffdynrel_trt_rescor_gamma_ls and
##    psyrat_rel_diffdynrel_trt_gamma_ls. It matters more than the one-facet
##    equivalent would: the reference bundle's own two-facet fits were prepared
##    against a SIMULATED second occasion and were never sampled, so there is no
##    completed external run to replay the two-facet code path against. These
##    rows are the only external check it has.
##
##    WHAT THIS PINS AND WHAT IT DOES NOT. It pins the FORMULAS - the six
##    observed-scale components, their cross-event covariances, the pooled
##    pio_e residual channel, the mean-matched quadrature, and the CE/CS/CES
##    coefficients. It CANNOT pin the Stan builder's facet-index assignment
##    (the reference's u_pi/u_po/u_io against the toolbox's tid/oid/to), because
##    stage 3 consumes SDs by name and never sees an index. That mapping needs a
##    live round-trip and remains open.
##
##    Both D-study designs are emitted. 'actual' uses each person's own
##    (n_i, n_o) - note the concurrent design measures both events on the same
##    trials, so n_i is shared across events rather than per-event.
##    'standardized' uses the settings values. THEY MUST DIFFER FROM THE ACTUAL
##    COUNTS AND FROM EACH OTHER: the original fixture set n' = (4, 2), which
##    coincided exactly with the simulated design's (I, O) = (4, 2), so the
##    'actual' and 'standardized' blocks came out byte-identical and the
##    design-selector divisors (pi/n_i', po/n_o', pio/(n_i'*n_o'), i/n_i',
##    o/n_o', io/(n_i'*n_o')) were unpinned - a swapped or mis-assigned divisor
##    would have passed. Distinct primes (7, 3), unequal to I = 4, O = 2 and to
##    each other, give every divisor slot a unique value.
## ---------------------------------------------------------------------------
rel_settings <- twofacet$settings
rel_settings$max_reliability_draws <- 0L   # keep every draw; no thinning
rel_settings$n_i_prime <- 7
rel_settings$n_o_prime <- 3

## Stage 3 additionally reads the POPULATION SDs and cross-event correlations,
## which stage 2 never touches (it works from the realized random effects alone).
## They are added HERE rather than in emit_stage2 on purpose: the stage-2 draws
## fixture enumerates every column of that frame, so adding them upstream would
## rewrite a fixture this section has no business changing.
##
## Every value is DISTINCT across events and across facets, and cor_po is
## negative. That is deliberate - with symmetric values an index swap between,
## say, pi and po would reproduce the right answer and the fixture would pin
## nothing. Magnitudes match the generating SDs of the realized effects above so
## the draw stays internally coherent.
rel_draws <- twofacet$draw
nd_rel <- nrow(rel_draws)
jitter <- seq(0, by = 0.004, length.out = nd_rel)   # per-draw variation
pop_pars <- list(
  p_joint = c(0.260, 0.230), i = c(0.110, 0.130), o = c(0.075, 0.062),
  pi = c(0.068, 0.081), po = c(0.059, 0.047), io = c(0.052, 0.066)
)
pop_cors <- c(p_mu = 0.42, i = 0.31, o = 0.24, pi = 0.19, po = -0.14, io = 0.27)
for (nm in names(pop_pars)) {
  for (m in 1:2) {
    rel_draws[[sprintf("sd_%s[%d]", nm, m)]] <- pop_pars[[nm]][m] + jitter
  }
}
for (nm in names(pop_cors)) {
  rel_draws[[sprintf("cor_%s", nm)]] <- pop_cors[[nm]] + jitter
}

rel_draws$rho_e <- twofacet$stage$rho_e
rel_draws$copula_stage_invalidity_flag <-
  twofacet$stage$copula_stage_invalidity_flag
rel_draws$pit_calibration_review_flag <-
  twofacet$stage$pit_calibration_review_flag

wcsv(twofacet$metadata$person_table,
     "modular_cut_copula_twofacet_persontable.csv")

## The stage-3-only parameters, emitted so the MATLAB test reads them rather
## than restating them as constants. Facet names here are the REFERENCE's
## (i = trial, o = occasion, pi = person x trial, po = person x occasion,
## io = trial x occasion); the test is responsible for mapping them onto the
## toolbox's trl/occ/tid/oid/to, and the deliberately distinct values above are
## what make a mis-mapping fail rather than pass by symmetry.
wcsv(do.call(rbind, c(
  lapply(names(pop_pars), function(nm) data.frame(
    param = sprintf("sd_%s", nm), event = rep(1:2, each = nd_rel),
    draw = rep(seq_len(nd_rel), 2),
    value = c(rel_draws[[sprintf("sd_%s[1]", nm)]],
              rel_draws[[sprintf("sd_%s[2]", nm)]]))),
  lapply(names(pop_cors), function(nm) data.frame(
    param = sprintf("cor_%s", nm), event = NA_integer_,
    draw = seq_len(nd_rel), value = rel_draws[[sprintf("cor_%s", nm)]])),
  list(data.frame(param = "rho_e", event = NA_integer_,
                  draw = seq_len(nd_rel), value = rel_draws$rho_e)))),
  "modular_cut_copula_twofacet_stage3_params.csv")

wcsv(compute_person_specific_reliability(rel_draws, rel_settings,
                                         twofacet$metadata, twofacet$standata),
     "modular_cut_copula_twofacet_stage3_expected.csv")

## ---------------------------------------------------------------------------
## 7. FULL stage-3 run, ONE-FACET.
##    The externally grounded UNIT-LANE oracle for
##    psyrat_ssrel_diffdynrel_rescor_gamma_ls (analysis 19). The one-facet path
##    already has the external replay against the owner's completed fits, but
##    that check is BUNDLE-GATED and skips wherever the reference bundle is
##    absent (CI included); these rows put an external pin inside the repo so
##    the unit lane never runs oracle-less. Same structure as section 6; the
##    population SDs are again added here, all-distinct, so an event or facet
##    swap fails rather than passing by symmetry.
##
##    n_i' = 7: prime and unequal to the simulated design's I = 5, so the
##    'standardized' rows differ from 'actual' and the n_i' divisor is pinned.
## ---------------------------------------------------------------------------
rel1_settings <- onefacet$settings
rel1_settings$max_reliability_draws <- 0L   # keep every draw; no thinning
rel1_settings$n_i_prime <- 7
rel1_settings$n_o_prime <- 1

rel1_draws <- onefacet$draw
nd_rel1 <- nrow(rel1_draws)
jitter1 <- seq(0, by = 0.004, length.out = nd_rel1)   # per-draw variation
pop_pars_1f <- list(p_joint = c(0.265, 0.235), i = c(0.112, 0.128))
pop_cors_1f <- c(p_mu = 0.41, i = 0.28)
for (nm in names(pop_pars_1f)) {
  for (m in 1:2) {
    rel1_draws[[sprintf("sd_%s[%d]", nm, m)]] <- pop_pars_1f[[nm]][m] + jitter1
  }
}
for (nm in names(pop_cors_1f)) {
  rel1_draws[[sprintf("cor_%s", nm)]] <- pop_cors_1f[[nm]] + jitter1
}

rel1_draws$rho_e <- onefacet$stage$rho_e
rel1_draws$copula_stage_invalidity_flag <-
  onefacet$stage$copula_stage_invalidity_flag
rel1_draws$pit_calibration_review_flag <-
  onefacet$stage$pit_calibration_review_flag

wcsv(onefacet$metadata$person_table,
     "modular_cut_copula_onefacet_persontable.csv")

wcsv(do.call(rbind, c(
  lapply(names(pop_pars_1f), function(nm) data.frame(
    param = sprintf("sd_%s", nm), event = rep(1:2, each = nd_rel1),
    draw = rep(seq_len(nd_rel1), 2),
    value = c(rel1_draws[[sprintf("sd_%s[1]", nm)]],
              rel1_draws[[sprintf("sd_%s[2]", nm)]]))),
  lapply(names(pop_cors_1f), function(nm) data.frame(
    param = sprintf("cor_%s", nm), event = NA_integer_,
    draw = seq_len(nd_rel1), value = rel1_draws[[sprintf("cor_%s", nm)]])),
  list(data.frame(param = "rho_e", event = NA_integer_,
                  draw = seq_len(nd_rel1), value = rel1_draws$rho_e)))),
  "modular_cut_copula_onefacet_stage3_params.csv")

wcsv(compute_person_specific_reliability(rel1_draws, rel1_settings,
                                         onefacet$metadata, onefacet$standata),
     "modular_cut_copula_onefacet_stage3_expected.csv")

## ---------------------------------------------------------------------------
## Provenance: exactly what produced these numbers. Every reference file this
## script sources is hashed - the sourcing loop above plus simulation_helpers.R
## - so a changed bundle is detectable from the fixture directory alone.
## ---------------------------------------------------------------------------
ref_files <- c(file.path("R", c("common.R", "model_generation.R",
                                "copula_stage.R", "reliability.R")),
               file.path("tests", "simulation_helpers.R"))
wcsv(data.frame(
  file = c(ref_files, "R_version", "RNG_kinds"),
  md5 = c(unname(tools::md5sum(file.path(code_dir, ref_files))),
          NA_character_, NA_character_),
  detail = c(rep("reference bundle source", length(ref_files)),
             R.version.string, paste(RNGkind(), collapse = " "))
), "modular_cut_copula_provenance.csv")

cat("\nAll modular cut-copula fixtures regenerated.\n")
