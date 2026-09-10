functions {
real partial_log_lik(array[] int seq_slice, int start, int end,
  vector meas, vector mu, vector sigma) {
  return normal_lpdf(meas[start:end] | mu[start:end], sigma[start:end]);
}
}
data {
int<lower=1> NOBS;  //number of observations
int<lower=1> grainsize;  // grainsize for reduce_sum
array[NOBS] int<lower=1, upper=NOBS> seq;
int<lower=1> NSUB;  //number of subjects
int<lower=1> NOCC;  //number of occasions
int<lower=1> NTRL;  //number of trials
int<lower=1> NTID;  // total unique for trial*id
int<lower=1> NOID;  // total unique for occ*id
int<lower=1> NTO;  // total unique for trial*occ
int<lower=1> KDIM;  //number of dimension predictors (1 or 3)
array[NOBS] int<lower=1, upper=NSUB> id;  //id variable
array[NOBS] int<lower=1, upper=NOCC> occ;  //occ variable
array[NOBS] int<lower=1, upper=NTRL> trl;  //trl variable
array[NOBS] int<lower=1, upper=NTID> trlxid;  //trlxid variable
array[NOBS] int<lower=1, upper=NOID> occxid;  //occxid variable
array[NOBS] int<lower=1, upper=NTO> trlxocc;  //trlxocc variable
vector[NOBS * KDIM] Xdim_vec;  //flattened (column-major) dimension design
vector[NOBS] meas;  // response variable
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
real Intercept;  // population intercept (mean at z = 0)
real Intercept_sigma;  // population sigma intercept (log residual at z = 0)
vector[KDIM] b;  // mean dimension slopes
vector[KDIM] b_sigma;  // log-residual (scale) dimension slopes
vector<lower=0>[2] gro_sds;  // group-level standard deviations (sigma_p, sigma_delta_p)
matrix[2, NSUB] gro_effs_stndzd;  // standardized group-level effects
cholesky_factor_corr[2] chol_corrmat;  // cholesky factor of correlation matrix
real<lower=0> sig_occ;  // occ-level std dev
real<lower=0> sig_trl;  // trl-level std dev
real<lower=0> sig_trlxid;  // trlxid std dev
real<lower=0> sig_occxid;  // occxid std dev
real<lower=0> sig_trlxocc;  // trlxocc std dev
vector[NOCC] occ_raw;  // occ means
vector[NTRL] trl_raw;  // trl means
vector[NTID] trlxid_raw;  // trlxid means
vector[NOID] occxid_raw;  // occxid means
vector[NTO] trlxocc_raw;  // trlxocc means
}
transformed parameters {
matrix[NSUB, 2] gro_effs_actual;  // actual group-level effects
// using vectors speeds up indexing in loops
vector[NSUB] ind_bs;
vector[NSUB] ind_sd;
vector[NOCC] occ_terms;  // occ-level terms
vector[NTRL] trl_terms;  // trl-level terms
vector[NTID] trlxid_terms;  // trlxid terms
vector[NOID] occxid_terms;  // occxid terms
vector[NTO] trlxocc_terms;  // trlxocc terms
// compute actual group-level effects
gro_effs_actual = (diag_pre_multiply(gro_sds, chol_corrmat) * gro_effs_stndzd)';
ind_bs = gro_effs_actual[, 1];
ind_sd = gro_effs_actual[, 2];
occ_terms = sig_occ * occ_raw;
trl_terms = sig_trl * trl_raw;
trlxid_terms = sig_trlxid * trlxid_raw;
occxid_terms = sig_occxid * occxid_raw;
trlxocc_terms = sig_trlxocc * trlxocc_raw;
}
model {
// dimension-conditioned linear predictors (fixed slopes on mean and scale)
vector[NOBS] mu = Intercept + Xdim * b;
vector[NOBS] sigma = Intercept_sigma + Xdim * b_sigma;
mu += ind_bs[id];
mu += occ_terms[occ] + trl_terms[trl] + trlxid_terms[trlxid]
  + occxid_terms[occxid] + trlxocc_terms[trlxocc];
sigma += ind_sd[id];
sigma = exp(sigma);
// priors including all constants
target += student_t_lpdf(Intercept | 3, 15.7, 12.4);
target += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
target += student_t_lpdf(gro_sds | 3, 0, 12.4)
  - 2 * student_t_lccdf(0 | 3, 0, 12.4);
target += std_normal_lpdf(to_vector(gro_effs_stndzd));
target += lkj_corr_cholesky_lpdf(chol_corrmat | 1);
target += std_normal_lpdf(occ_raw);
target += std_normal_lpdf(trl_raw);
target += std_normal_lpdf(trlxid_raw);
target += std_normal_lpdf(occxid_raw);
target += std_normal_lpdf(trlxocc_raw);
target += normal_lpdf(b | 0, 10);
target += normal_lpdf(b_sigma | 0, 1);
target += cauchy_lpdf(sig_occ | 0, 5);
target += cauchy_lpdf(sig_trl | 0, 10);
target += cauchy_lpdf(sig_trlxid | 0, 5);
target += cauchy_lpdf(sig_occxid | 0, 5);
target += cauchy_lpdf(sig_trlxocc | 0, 1);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize, meas, mu, sigma);
}
