data {
int<lower=1> NOBS;  //number of observations
int<lower=1> NSUB;  //number of subjects
int<lower=1> NTRL;  //number of trials
int<lower=1> KDIM;  //number of dimension predictors (1 or 3)
array[NOBS] int<lower=1> id;  //id variable
array[NOBS] int<lower=1, upper=NTRL> trl;  //trial variable
vector[NOBS * KDIM] Xdim_vec;  //flattened (column-major) dimension design
vector[NOBS] meas;  // response variable
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
real Intercept;  // population intercept (residual scale at z = 0)
real Intercept_sigma;  // population sigma intercept (log residual at z = 0)
vector[KDIM] b;  // mean dimension slopes
vector[KDIM] b_sigma;  // log-residual (scale) dimension slopes
vector<lower=0>[2] gro_sds;  // group-level standard deviations (sigma_p, sigma_delta_p)
matrix[2, NSUB] gro_effs_stndzd;  // standardized group-level effects
cholesky_factor_corr[2] chol_corrmat;  // cholesky factor of correlation matrix
real<lower=0> sig_trl;  // trial main-effect SD (sigma_i)
vector[NTRL] trl_raw;  // standardized trial main effects
}
transformed parameters {
matrix[NSUB, 2] gro_effs_actual;  // actual group-level effects
// using vectors speeds up indexing in loops
vector[NSUB] ind_bs;
vector[NSUB] ind_sd;
vector[NTRL] trl_b;  // non-centered trial main effects
// compute actual group-level effects
gro_effs_actual = (diag_pre_multiply(gro_sds, chol_corrmat) * gro_effs_stndzd)';
ind_bs = gro_effs_actual[, 1];
ind_sd = gro_effs_actual[, 2];
trl_b = sig_trl * trl_raw;
}
model {
// dimension-conditioned linear predictors (fixed slopes on mean and scale)
vector[NOBS] mu = Intercept + Xdim * b;
vector[NOBS] sigma = Intercept_sigma + Xdim * b_sigma;
mu += ind_bs[id];
mu += trl_b[trl];
sigma += ind_sd[id];
sigma = exp(sigma);
// priors including all constants
target += student_t_lpdf(Intercept | 3, 15.7, 12.4);
target += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
target += student_t_lpdf(gro_sds | 3, 0, 12.4)
  - 2 * student_t_lccdf(0 | 3, 0, 12.4);
target += std_normal_lpdf(to_vector(gro_effs_stndzd));
target += lkj_corr_cholesky_lpdf(chol_corrmat | 1);
target += std_normal_lpdf(trl_raw);
target += normal_lpdf(b | 0, 10);
target += normal_lpdf(b_sigma | 0, 1);
target += cauchy_lpdf(sig_trl | 0, 10);
// likelihood including all constants
target += normal_lpdf(meas | mu, sigma);
}
