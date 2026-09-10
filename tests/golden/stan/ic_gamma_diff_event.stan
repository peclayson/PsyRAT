functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
}
data {
int<lower=1> NOBS;  // number of observations
vector[NOBS] meas;  // response variable
matrix[NOBS, 2] X;  // population-level design matrix
int<lower=1> NSUB;  // number of grouping levels
array[NOBS] int<lower=1> ID;  // grouping indicator per observation
int<lower=1> NTRL;  // number of grouping levels
array[NOBS] int<lower=1> TRL;  // grouping indicator per observation
vector[NOBS] meas1;
vector[NOBS] meas2;
}
parameters {
vector[2] b;  // per-event log-mean intercepts
vector[2] b_nu;  // per-event log-nu intercepts
vector<lower=0>[4] sd_id;  // person SDs: [mean_1, mean_2, lognu_1, lognu_2]
matrix[4, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[4] L_id;  // person correlation (means + log-nus)
vector<lower=0>[2] sd_trl;  // per-event trial main-effect log-SDs
matrix[2, NTRL] z_trl;  // standardized trial effects
cholesky_factor_corr[2] L_trl;  // cross-event trial correlation
}
transformed parameters {
matrix[NSUB, 4] r_id;
matrix[NTRL, 2] r_trl;
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
}
model {
vector[NOBS] mu = X * b;
vector[NOBS] nu = X * b_nu;
for (n in 1:NOBS) {
  mu[n] += r_id[ID[n], 1] * meas1[n] + r_id[ID[n], 2] * meas2[n]
        + r_trl[TRL[n], 1] * meas1[n] + r_trl[TRL[n], 2] * meas2[n];
  nu[n] += r_id[ID[n], 3] * meas1[n] + r_id[ID[n], 4] * meas2[n];
}
mu = exp(mu);
nu = exp(nu);
// priors including all constants
to_vector(z_id) ~ std_normal();
to_vector(z_trl) ~ std_normal();
b ~ normal(1.70475, 0.5);
b_nu ~ normal(2.89037, 0.75);
sd_id[1] ~ student_t(3, 0, 0.35);
sd_id[2] ~ student_t(3, 0, 0.35);
sd_id[3] ~ student_t(3, 0, 1);
sd_id[4] ~ student_t(3, 0, 1);
sd_trl ~ student_t(3, 0, 0.35);
L_id ~ lkj_corr_cholesky(1);
L_trl ~ lkj_corr_cholesky(1);
// likelihood including all constants
for (n in 1:NOBS) {
  target += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
}
}
generated quantities {
corr_matrix[4] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
real cor_p = cor_id[1, 2];  // cross-event person-mean correlation
real cor_i = cor_trl[1, 2];  // cross-event trial correlation
}
