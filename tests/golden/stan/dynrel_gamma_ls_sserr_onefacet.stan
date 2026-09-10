functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
}
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
real Intercept;  // log expected-score grand intercept (alpha, mean at z = 0)
real Intercept_sigma;  // log residual-SD grand intercept (at z = 0)
vector[KDIM] b;  // log-mean dimension slopes
vector[KDIM] b_sigma;  // log residual-SD dimension slopes
vector<lower=0>[2] gro_sds;  // person log-scale SDs: [mean (s_p), log-sigma]
matrix[2, NSUB] gro_effs_stndzd;  // standardized person effects
cholesky_factor_corr[2] chol_corrmat;  // person (mean, log-sigma) correlation
real<lower=0> sig_trl;  // trial main-effect SD on the log mean (sigma_i,log)
vector[NTRL] trl_raw;  // standardized trial main effects
}
transformed parameters {
matrix[NSUB, 2] gro_effs_actual;  // actual person effects
vector[NSUB] ind_bs;  // person log-mean effect (u_p)
vector[NSUB] ind_sd;  // person log-residual-SD effect (w_p)
vector[NTRL] trl_b;  // non-centered trial main effects
gro_effs_actual = (diag_pre_multiply(gro_sds, chol_corrmat) * gro_effs_stndzd)';
ind_bs = gro_effs_actual[, 1];
ind_sd = gro_effs_actual[, 2];
trl_b = sig_trl * trl_raw;
}
model {
// dimension-conditioned, log-linked mean: log mu = alpha + Xdim*b + person + trial
vector[NOBS] mu = Intercept + Xdim * b;
// dimension-conditioned, log-linked, person-varying residual SD:
// log sigma = Intercept_sigma + Xdim*b_sigma + person (NO trial term, so the
// residual stays a single sigma_pi,e(z) per participant)
vector[NOBS] sigma = Intercept_sigma + Xdim * b_sigma;
vector[NOBS] nu;  // derived below; declared here so all declarations precede statements
mu += ind_bs[id];
mu += trl_b[trl];
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
target += student_t_lpdf(sig_trl | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += std_normal_lpdf(to_vector(gro_effs_stndzd));
target += lkj_corr_cholesky_lpdf(chol_corrmat | 1);
target += std_normal_lpdf(trl_raw);
// likelihood including all constants
for (n in 1:NOBS) {
  target += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
}
}
