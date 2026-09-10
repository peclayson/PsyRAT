#!/usr/bin/env Rscript
##
## Regenerates the frozen external-oracle fixtures for the person-specific
## dynamic NONCONCURRENT difference-of-differences designs (PsyRAT analyses
## 28 one-facet / 29 two-facet).
##
## WHY THESE FIXTURES ARE DIFFERENT IN KIND. PsyRAT's accuracy suite is
## otherwise regression-only: its oracle duplicates the production formulas, so
## a green run proves "unchanged", not "correct". The CSVs this script writes
## come from the owner's INDEPENDENT R reference implementation
## (person_specific_dynamic_nonconcurrent_dod_modular), so a MATLAB test that
## matches them is externally grounded rather than self-referential.
##
## Usage (from anywhere):
##
##   Rscript --vanilla make_nonconcurrent_dod_fixtures.R [output_dir]
##
## Point MPSDOD_REFERENCE_DIR at the reference bundle if it is not at the
## default path below. The script is fully deterministic - no RNG is consumed
## after the data simulation, and every posterior value is an explicit literal
## plus an arithmetic jitter - so rerunning it on an unchanged reference bundle
## reproduces every file byte for byte. An unexpected diff means the reference
## changed and the MATLAB tolerances should be re-examined rather than the
## fixtures silently re-baselined.
##
## Regenerate ONLY from the reference bundle. Never hand-edit a fixture, and
## never regenerate them from the MATLAB implementation they exist to check -
## that would convert an external oracle into a regression snapshot.
##
## SCOPE NOTE (owner ruling 2026-08-12): the toolbox implements the HARMONIC
## count rule only; the reference also emits exact-shared-count sensitivity
## columns (*_exact_shared). Those columns are frozen here anyway - they come
## for free from the reference and let a later increment add the sensitivity
## without re-freezing - but the MATLAB tests compare harmonic-primary columns
## only.

options(stringsAsFactors = FALSE)

## Pin the RNG kinds so the data-simulation step is byte-for-byte stable across
## R releases. Must precede every set.seed() (inside the simulator).
RNGversion("3.6.0")

code_dir <- Sys.getenv(
  "MPSDOD_REFERENCE_DIR",
  "~/Desktop/stan_chisquare/person_specific_dynamic_nonconcurrent_dod_modular"
)
code_dir <- path.expand(code_dir)
if (!file.exists(file.path(code_dir, "R", "independence_stage.R"))) {
  stop("Reference bundle not found at: ", code_dir,
       "\nSet MPSDOD_REFERENCE_DIR to its location.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else {
  file_arg <- sub("^--file=", "",
                  commandArgs(trailingOnly = FALSE)[
                    grepl("^--file=", commandArgs(trailingOnly = FALSE))])
  if (length(file_arg)) dirname(normalizePath(file_arg[[1]])) else getwd()
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref_files <- c(file.path("R", c("common.R", "model_generation.R",
                                "independence_stage.R", "reliability.R")),
               file.path("tests", "simulation_helpers.R"))
for (f in ref_files) source(file.path(code_dir, f))

wcsv <- function(d, name) {
  ## A binary connection with an explicit eol pins the line endings to LF on
  ## every platform; a plain path would go through a text-mode connection and
  ## produce CRLF fixtures on Windows, which the frozen git bytes are not.
  con <- file(file.path(out_dir, name), open = "wb")
  on.exit(close(con), add = TRUE)
  utils::write.csv(d, con, row.names = FALSE, na = "", eol = "\n")
  cat("wrote", name, "-", nrow(d), "rows\n")
}

## Deterministic correlation matrix from explicit factor loadings. L %*% t(L)
## plus a positive diagonal is positive definite by construction, so cov2cor of
## it is a valid correlation matrix whatever loadings are chosen; mixed-sign
## loadings put NEGATIVE correlations in the result, which is what makes an
## index or sign error fail rather than pass by symmetry. The small per-draw
## diagonal inflation keeps the matrices distinct across draws (a reader that
## reused draw 1's matrix for every draw would otherwise pass).
det_corr <- function(loadings, dvec, inflate) {
  stats::cov2cor(loadings %*% t(loadings) + diag(dvec + inflate))
}

## ---------------------------------------------------------------------------
## Hand-constructed posterior draws. Three draws; every parameter is an
## explicit literal plus jitter = 0, 0.004, 0.008 per draw, so all values are
## distinct across cells, facets, persons AND draws. All-distinct is the load-
## bearing property (the copula fixture script's doctrine): symmetric values
## would let a swapped index reproduce the right answer.
## ---------------------------------------------------------------------------
ND <- 3L
jit <- c(0, 0.004, 0.008)

build_draws <- function(standata, facet) {
  K <- standata$K
  P <- standata$P
  draw <- data.frame(.draw = seq_len(ND))

  alpha_mu_base    <- c(0.055, -0.110, 0.020, 0.085)
  alpha_sigma_base <- c(-0.070, 0.130, 0.045, -0.025)
  ## Explicit per-cell dimension slopes on BOTH submodels (columns = predictor
  ## terms: K = 1 -> Dim1_z; K = 3 -> Dim1_z, Dim2_z, Dim1_z*Dim2_z).
  beta_mu_base <- matrix(c( 0.115, -0.062,  0.041, -0.088,
                            -0.053,  0.097, -0.076,  0.034,
                             0.029, -0.045,  0.083, -0.017), nrow = 4)[, seq_len(K), drop = FALSE]
  beta_sigma_base <- matrix(c(-0.048,  0.071, -0.033,  0.056,
                               0.062, -0.027,  0.049, -0.081,
                              -0.019,  0.038, -0.064,  0.023), nrow = 4)[, seq_len(K), drop = FALSE]

  for (q in 1:4) {
    draw[[sprintf("alpha_mu[%d]", q)]] <- alpha_mu_base[q] + jit
    draw[[sprintf("alpha_sigma[%d]", q)]] <- alpha_sigma_base[q] + jit
    for (k in seq_len(K)) {
      draw[[sprintf("beta_mu[%d,%d]", q, k)]] <- beta_mu_base[q, k] + jit
      draw[[sprintf("beta_sigma[%d,%d]", q, k)]] <- beta_sigma_base[q, k] + jit
    }
  }

  ## Person block: 4 location SDs then 4 log-sigma SDs, all distinct.
  sd_p_base <- c(0.260, 0.235, 0.210, 0.185, 0.115, 0.095, 0.130, 0.105)
  for (j in 1:8) draw[[sprintf("sd_p_joint[%d]", j)]] <- sd_p_base[j] + jit / 2

  ## 8 x 2 loadings with mixed signs -> negative correlations present in the
  ## joint person correlation matrix (e.g. between location 2 and scale 6).
  L8 <- matrix(c( 0.55,  0.30,
                  0.48, -0.35,
                 -0.42,  0.28,
                  0.36,  0.22,
                  0.25, -0.35,
                 -0.35, -0.25,
                  0.28,  0.32,
                 -0.22,  0.38), nrow = 8, byrow = TRUE)
  d8 <- seq(0.55, 0.90, length.out = 8)
  for (d in seq_len(ND)) {
    R8 <- det_corr(L8, d8, 0.05 * (d - 1))
    for (i in 1:8) for (j in 1:8) {
      nm <- sprintf("R_p_joint[%d,%d]", i, j)
      if (d == 1L) draw[[nm]] <- numeric(ND)
      draw[[nm]][d] <- R8[i, j]
    }
  }

  ## Person realized effects: explicit, all-distinct, nonzero in BOTH the
  ## location (1:4) and log-sigma (5:8) halves, varying across draws.
  base_p <- seq(-0.25, 0.25, length.out = P)
  col_scale <- c(0.9, -1.1, 1.25, -0.8, 1.05, -1.3, 0.7, 1.15) * 0.2
  for (p in seq_len(P)) for (j in 1:8) {
    draw[[sprintf("u_p_joint[%d,%d]", p, j)]] <-
      base_p[p] * col_scale[j] + 0.011 * j * (-1)^p + jit * (-1)^j
  }

  ## Facet blocks: per-cell SDs distinct within and across facets; each
  ## correlation matrix from its own loading set with at least one negative
  ## off-diagonal.
  facet_sds <- list(
    i  = c(0.112, 0.128, 0.098, 0.121),
    o  = c(0.075, 0.062, 0.088, 0.070),
    pi = c(0.068, 0.081, 0.059, 0.074),
    po = c(0.059, 0.047, 0.066, 0.052),
    io = c(0.052, 0.066, 0.044, 0.058)
  )
  facet_loadings <- list(
    i  = matrix(c( 0.50,  0.25, -0.38,  0.30,  0.27, -0.33,  0.36,  0.22), nrow = 4),
    o  = matrix(c( 0.44, -0.29,  0.35,  0.26, -0.31,  0.28,  0.24,  0.37), nrow = 4),
    pi = matrix(c( 0.41,  0.33,  0.29, -0.36,  0.30,  0.26, -0.28,  0.24), nrow = 4),
    po = matrix(c(-0.39,  0.31,  0.27,  0.34,  0.25, -0.29,  0.32,  0.21), nrow = 4),
    io = matrix(c( 0.37,  0.28, -0.34,  0.25,  0.29,  0.35,  0.23, -0.27), nrow = 4)
  )
  d4 <- seq(0.60, 0.85, length.out = 4)
  effects <- if (facet == "trial") "i" else c("i", "o", "pi", "po", "io")
  for (e in effects) {
    for (q in 1:4) {
      draw[[sprintf("sd_%s[%d]", e, q)]] <- facet_sds[[e]][q] + jit / 2
    }
    for (d in seq_len(ND)) {
      R4 <- det_corr(facet_loadings[[e]], d4, 0.05 * (d - 1))
      for (i in 1:4) for (j in 1:4) {
        nm <- sprintf("R_%s[%d,%d]", e, i, j)
        if (d == 1L) draw[[nm]] <- numeric(ND)
        draw[[nm]][d] <- R4[i, j]
      }
    }
  }
  draw
}

## Long format so MATLAB reassembles without a column-order convention.
draws_long <- function(draw) {
  do.call(rbind, lapply(setdiff(names(draw), ".draw"), function(nm) {
    base <- sub("\\[.*$", "", nm)
    idx <- regmatches(nm, regexpr("\\[.*\\]", nm))
    parts <- if (length(idx)) as.integer(strsplit(gsub("\\[|\\]", "", idx), ",")[[1]]) else integer(0)
    data.frame(draw = seq_len(ND), param = base,
               idx1 = if (length(parts) >= 1) parts[1] else NA_integer_,
               idx2 = if (length(parts) >= 2) parts[2] else NA_integer_,
               value = draw[[nm]])
  }))
}

## ---------------------------------------------------------------------------
## One lane per facet structure. Standardized n-primes are DISTINCT PRIMES,
## unequal to the simulated design's counts (I_error = 5, I_correct = 7,
## O = 2) and to each other, so every D-study divisor slot gets a unique value
## and a swapped or mis-assigned divisor fails rather than passing by
## coincidence (the lesson recorded in the copula fixture script at its
## stage-3 section).
## ---------------------------------------------------------------------------
emit_lane <- function(facet, dimension, P, O, seed, tag) {
  st <- mpsdod_settings()
  st$facet <- facet
  st$dimension <- dimension
  st$reliability_cores <- 1L
  st$max_reliability_draws <- 0L      # keep every draw; no thinning
  st$allow_synthetic_occasions <- (O > 1L)
  st$n_i_error <- 11; st$n_i_correct <- 13
  st$n_o_error <- 3;  st$n_o_correct <- 5
  ## Exact-shared settings must still VALIDATE (the reference computes those
  ## columns); they are frozen but not compared by the toolbox (harmonic-only
  ## scope). Each is strictly below its bound so the exact-shared columns
  ## differ from the harmonic ones in the frozen output.
  st$n_common_i <- 4; st$n_common_o <- 2; st$n_common_io <- 30

  tmp <- tempfile(paste0("mpsdod_", tag, "_")); dir.create(tmp)
  path <- file.path(tmp, "d.rds")
  saveRDS(simulate_nonconcurrent_dod_source(P = P, I_error = 5L, I_correct = 7L,
                                            O = O, seed = seed), path)
  prep <- prepare_modular_dod_data(path, st)
  sd_ <- prep$standata

  draw <- build_draws(sd_, facet)
  wcsv(draws_long(draw), sprintf("nonconcurrent_dod_%s_draws.csv", tag))

  obs <- data.frame(y_theta = sd_$y_theta, y_alpha = sd_$y_alpha,
                    condition_id = sd_$condition_id,
                    p_id = sd_$p_id, i_id = sd_$i_id)
  for (k in seq_len(sd_$K)) obs[[sprintf("Z%d", k)]] <- sd_$Z[, k]
  if (facet == "trial_occasion") {
    obs$o_id <- sd_$o_id; obs$pi_id <- sd_$pi_id
    obs$po_id <- sd_$po_id; obs$io_id <- sd_$io_id
  }
  wcsv(obs, sprintf("nonconcurrent_dod_%s_standata.csv", tag))

  wcsv(data.frame(
    facet = facet, dimension = dimension, N = sd_$N, K = sd_$K, P = sd_$P,
    I = sd_$I,
    O = if (is.null(sd_$O)) NA_integer_ else sd_$O,
    PI = if (is.null(sd_$PI)) NA_integer_ else sd_$PI,
    PO = if (is.null(sd_$PO)) NA_integer_ else sd_$PO,
    IO = if (is.null(sd_$IO)) NA_integer_ else sd_$IO,
    mu_offset_1 = sd_$mu_offset[1], mu_offset_2 = sd_$mu_offset[2],
    mu_offset_3 = sd_$mu_offset[3], mu_offset_4 = sd_$mu_offset[4],
    log_sigma_offset_1 = sd_$log_sigma_offset[1],
    log_sigma_offset_2 = sd_$log_sigma_offset[2],
    log_sigma_offset_3 = sd_$log_sigma_offset[3],
    log_sigma_offset_4 = sd_$log_sigma_offset[4],
    n_draws = ND,
    n_i_error_std = st$n_i_error, n_i_correct_std = st$n_i_correct,
    n_o_error_std = st$n_o_error, n_o_correct_std = st$n_o_correct,
    n_common_i_std = st$n_common_i, n_common_o_std = st$n_common_o,
    n_common_io_std = st$n_common_io,
    numerical_tolerance = st$numerical_tolerance,
    w1 = st$contrast_weights[1], w2 = st$contrast_weights[2],
    w3 = st$contrast_weights[3], w4 = st$contrast_weights[4],
    outcome_1 = st$outcome_labels[1], outcome_2 = st$outcome_labels[2],
    outcome_3 = st$outcome_labels[3], outcome_4 = st$outcome_labels[4],
    inference_framework = st$inference_framework,
    estimand_version = st$estimand_version
  ), sprintf("nonconcurrent_dod_%s_meta.csv", tag))

  wcsv(prep$metadata$person_table,
       sprintf("nonconcurrent_dod_%s_persontable.csv", tag))

  audit <- validate_residual_independence_boundary(st, prep$metadata)
  wcsv(audit$manifest, sprintf("nonconcurrent_dod_%s_residual_audit.csv", tag))

  wcsv(compute_person_specific_reliability(draw, st, prep$metadata, sd_),
       sprintf("nonconcurrent_dod_%s_stage3_expected.csv", tag))

  unlink(tmp, recursive = TRUE)
  invisible(NULL)
}

emit_lane("trial",          "one_dimension",  6L, 1L, 90210L, "onefacet")
emit_lane("trial_occasion", "two_dimensions", 5L, 2L, 13579L, "twofacet")

## ---------------------------------------------------------------------------
## Provenance: exactly what produced these numbers. Every reference file this
## script sources is hashed, so a changed bundle is detectable from the
## fixture directory alone.
## ---------------------------------------------------------------------------
wcsv(data.frame(
  file = c(ref_files, "R_version", "RNG_kinds"),
  md5 = c(unname(tools::md5sum(file.path(code_dir, ref_files))),
          NA_character_, NA_character_),
  detail = c(rep("reference bundle source", length(ref_files)),
             R.version.string, paste(RNGkind(), collapse = " "))
), "nonconcurrent_dod_provenance.csv")

cat("\nAll nonconcurrent DoD fixtures regenerated.\n")
