functions {
  real biv_scaled_chisq_copula_log_density(real y1, real y2,
    real mu1, real mu2, real nu1, real nu2, real rho) {
    if (y1 <= 0 || y2 <= 0) {
      return negative_infinity();
    }
    real shape1 = 0.5 * nu1;
    real shape2 = 0.5 * nu2;
    real rate1 = 0.5 * nu1 / mu1;
    real rate2 = 0.5 * nu2 / mu2;
    real u1 = gamma_cdf(y1 | shape1, rate1);
    real u2 = gamma_cdf(y2 | shape2, rate2);
    real z1;
    real z2;
    real one_minus_r2 = 1 - square(rho);
    real log_c;
    u1 = fmin(1 - 1e-10, fmax(1e-10, u1));
    u2 = fmin(1 - 1e-10, fmax(1e-10, u2));
    z1 = inv_Phi(u1);
    z2 = inv_Phi(u2);
    log_c = -0.5 * log(one_minus_r2)
      - 0.5 * (square(rho) * (square(z1) + square(z2)) - 2 * rho * z1 * z2) / one_minus_r2;
    return gamma_lpdf(y1 | shape1, rate1) + gamma_lpdf(y2 | shape2, rate2) + log_c;
  }
  real partial_log_lik(array[] int seq_slice, int start, int end,
    vector meas1, vector meas2, vector mu1, vector mu2,
    vector nu1, vector nu2, real rho_e) {
    real lp = 0;
    for (n in start:end) {
      lp += biv_scaled_chisq_copula_log_density(meas1[n], meas2[n],
        mu1[n], mu2[n], nu1[n], nu2[n], rho_e);
    }
    return lp;
  }
}
data {
int<lower=1> NOBS;  // number of paired observations
int<lower=1> grainsize;  // grainsize for reduce_sum
array[NOBS] int<lower=1,upper=NOBS> seq;
vector[NOBS] meas1;  // response variable for event 1
vector[NOBS] meas2;  // response variable for event 2
int<lower=1> NSUB;  // number of grouping levels
array[NOBS] int<lower=1> ID;  // grouping indicator per observation
int<lower=1> NTRL;  // number of grouping levels
array[NOBS] int<lower=1> TRL;  // grouping indicator per observation
}
parameters {
vector[2] b;  // per-event log-mean intercepts
vector[2] b_nu;  // per-event GLOBAL log-nu intercepts (fixed dispersion)
vector<lower=0>[2] sd_id;  // person mean SDs: [mean_1, mean_2]
matrix[2, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[2] L_id;  // cross-event person-mean correlation
vector<lower=0>[2] sd_trl;  // per-event trial main-effect log-SDs
matrix[2, NTRL] z_trl;  // standardized trial effects
cholesky_factor_corr[2] L_trl;  // cross-event trial correlation
real<lower=-0.95, upper=0.95> rho_e;  // Gaussian-copula residual correlation
}
transformed parameters {
matrix[NSUB, 2] r_id;
matrix[NTRL, 2] r_trl;
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
}
model {
vector[NOBS] mu1;
vector[NOBS] mu2;
vector[NOBS] nu1;
vector[NOBS] nu2;
for (n in 1:NOBS) {
  mu1[n] = b[1] + r_id[ID[n], 1] + r_trl[TRL[n], 1];
  mu2[n] = b[2] + r_id[ID[n], 2] + r_trl[TRL[n], 2];
  nu1[n] = b_nu[1];
  nu2[n] = b_nu[2];
}
mu1 = exp(mu1);
mu2 = exp(mu2);
nu1 = exp(nu1);
nu2 = exp(nu2);
// priors including all constants
to_vector(z_id) ~ std_normal();
to_vector(z_trl) ~ std_normal();
b ~ normal(1.70475, 0.5);
b_nu ~ normal(2.89037, 0.75);
sd_id ~ student_t(3, 0, 0.35);
sd_trl ~ student_t(3, 0, 0.35);
L_id ~ lkj_corr_cholesky(1);
L_trl ~ lkj_corr_cholesky(1);
rho_e ~ normal(0, 0.5);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize, meas1, meas2, mu1, mu2, nu1, nu2, rho_e);
}
generated quantities {
corr_matrix[2] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
real cor_p = cor_id[1, 2];  // cross-event person-mean correlation
real cor_i = cor_trl[1, 2];  // cross-event trial correlation
}
