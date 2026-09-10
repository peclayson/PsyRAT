functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
  real partial_log_lik(array[] int seq_slice, int start, int end,
    vector meas, vector mu, vector nu) {
    real lp = 0;
    for (n in start:end) {
      lp += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
    }
    return lp;
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
real Intercept;  // log expected-score grand intercept (alpha, mean at z = 0)
real Intercept_sigma;  // log residual-SD grand intercept (at z = 0)
vector[KDIM] b;  // log-mean dimension slopes
vector[KDIM] b_sigma;  // log residual-SD dimension slopes
vector<lower=0>[2] gro_sds;  // person log-scale SDs: [mean (s_p), log-sigma]
matrix[2, NSUB] gro_effs_stndzd;  // standardized person effects
cholesky_factor_corr[2] chol_corrmat;  // person (mean, log-sigma) correlation
real<lower=0> sig_occ;  // occasion main-effect SD on the log mean
real<lower=0> sig_trl;  // trial main-effect SD on the log mean
real<lower=0> sig_trlxid;  // trial x person SD on the log mean
real<lower=0> sig_occxid;  // occasion x person SD on the log mean
real<lower=0> sig_trlxocc;  // trial x occasion SD on the log mean
vector[NOCC] occ_raw;  // standardized occasion main effects
vector[NTRL] trl_raw;  // standardized trial main effects
vector[NTID] trlxid_raw;  // standardized trial x person effects
vector[NOID] occxid_raw;  // standardized occasion x person effects
vector[NTO] trlxocc_raw;  // standardized trial x occasion effects
}
transformed parameters {
matrix[NSUB, 2] gro_effs_actual;  // actual person effects
vector[NSUB] ind_bs;  // person log-mean effect (u_p)
vector[NSUB] ind_sd;  // person log-residual-SD effect (w_p)
vector[NOCC] occ_terms;  // occasion main effects
vector[NTRL] trl_terms;  // trial main effects
vector[NTID] trlxid_terms;  // trial x person effects
vector[NOID] occxid_terms;  // occasion x person effects
vector[NTO] trlxocc_terms;  // trial x occasion effects
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
// dimension-conditioned, log-linked mean: log mu = alpha + Xdim*b + person
// + occasion + trial + interactions
vector[NOBS] mu = Intercept + Xdim * b;
// dimension-conditioned, log-linked, person-varying residual SD:
// log sigma = Intercept_sigma + Xdim*b_sigma + person (NO facet terms, so the
// residual stays a single sigma_pto,e(z) per participant)
vector[NOBS] sigma = Intercept_sigma + Xdim * b_sigma;
vector[NOBS] nu;  // derived below; declared here so all declarations precede statements
mu += ind_bs[id];
mu += occ_terms[occ] + trl_terms[trl] + trlxid_terms[trlxid]
  + occxid_terms[occxid] + trlxocc_terms[trlxocc];
mu = exp(mu);
sigma += ind_sd[id];
sigma = exp(sigma);
// moment-match the scaled chi-square to (mean, SD): Var(Y | mu, sigma) = sigma^2
nu = 2 * square(mu ./ sigma);
// priors including all constants
target += normal_lpdf(Intercept | 1.70475, 0.5);
target += normal_lpdf(Intercept_sigma | 0.606136, 0.5);
target += normal_lpdf(b | 0, 0.5);
target += normal_lpdf(b_sigma | 0, 0.25);
target += student_t_lpdf(gro_sds[1] | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += student_t_lpdf(gro_sds[2] | 3, 0, 0.5)
  - 1 * student_t_lccdf(0 | 3, 0, 0.5);
target += student_t_lpdf(sig_occ | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += student_t_lpdf(sig_trl | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += student_t_lpdf(sig_trlxid | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += student_t_lpdf(sig_occxid | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += student_t_lpdf(sig_trlxocc | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += std_normal_lpdf(to_vector(gro_effs_stndzd));
target += lkj_corr_cholesky_lpdf(chol_corrmat | 1);
target += std_normal_lpdf(occ_raw);
target += std_normal_lpdf(trl_raw);
target += std_normal_lpdf(trlxid_raw);
target += std_normal_lpdf(occxid_raw);
target += std_normal_lpdf(trlxocc_raw);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize, meas, mu, nu);
}
