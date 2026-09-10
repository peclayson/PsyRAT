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
array[NOBS] int<lower=1> id;  //id variable
array[NOBS] int<lower=1, upper=NTRL> trl;  //trial variable
vector[NOBS] meas;  // response variable
}
parameters {
real Intercept;  // log expected-score grand intercept (alpha)
real Intercept_nu;  // log degrees-of-freedom grand intercept (beta)
vector<lower=0>[2] gro_sds;  // person log-scale SDs: [mean (s_p), log-nu]
matrix[2, NSUB] gro_effs_stndzd;  // standardized person effects
cholesky_factor_corr[2] chol_corrmat;  // person (mean, log-nu) correlation
real<lower=0> sig_trl;  // trial main-effect SD on the log mean (sigma_i,log)
vector[NTRL] trl_raw;  // standardized trial main effects
}
transformed parameters {
matrix[NSUB, 2] gro_effs_actual;  // actual person effects
vector[NSUB] ind_bs;  // person log-mean effect (u_p)
vector[NSUB] ind_nu;  // person log-nu effect (v_p)
vector[NTRL] trl_b;  // non-centered trial main effects
gro_effs_actual = (diag_pre_multiply(gro_sds, chol_corrmat) * gro_effs_stndzd)';
ind_bs = gro_effs_actual[, 1];
ind_nu = gro_effs_actual[, 2];
trl_b = sig_trl * trl_raw;
}
model {
// log-linked mean: log mu = alpha + person + trial
vector[NOBS] mu = Intercept + rep_vector(0, NOBS);
// log-linked, person-varying dispersion: log nu = beta + person
vector[NOBS] nu = Intercept_nu + rep_vector(0, NOBS);
mu += ind_bs[id];
mu += trl_b[trl];
mu = exp(mu);
nu += ind_nu[id];
nu = exp(nu);
// priors including all constants
target += normal_lpdf(Intercept | 1.70475, 0.5);
target += normal_lpdf(Intercept_nu | 2.89037, 0.75);
target += student_t_lpdf(gro_sds[1] | 3, 0, 0.35)
  - 1 * student_t_lccdf(0 | 3, 0, 0.35);
target += student_t_lpdf(gro_sds[2] | 3, 0, 1)
  - 1 * student_t_lccdf(0 | 3, 0, 1);
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
