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
array[NOBS] int<lower=1,upper=4> cell;  // cell index (1-4)
int<lower=1> NSUB;  // number of persons
array[NOBS] int<lower=1> ID;  // person indicator per observation
int<lower=1> NTRL;  // number of trials
array[NOBS] int<lower=1> TRL;  // trial indicator per observation
}
parameters {
vector[4] b;  // per-cell log-mean intercepts
vector[4] b_nu;  // per-cell log-nu intercepts
vector<lower=0>[8] sd_id;  // person SDs: [mean_1..4, lognu_1..4]
matrix[8, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[8] L_id;  // person correlation (means + log-nus)
vector<lower=0>[4] sd_trl;  // per-cell trial main-effect log-SDs
matrix[4, NTRL] z_trl;  // standardized trial effects
cholesky_factor_corr[4] L_trl;  // cross-cell trial correlation
}
transformed parameters {
matrix[NSUB, 8] r_id;
matrix[NTRL, 4] r_trl;
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
}
model {
vector[NOBS] mu;
vector[NOBS] nu;
for (n in 1:NOBS) {
  int c = cell[n];
  mu[n] = b[c] + r_id[ID[n], c] + r_trl[TRL[n], c];
  nu[n] = b_nu[c] + r_id[ID[n], c + 4];
}
mu = exp(mu);
nu = exp(nu);
// priors including all constants
to_vector(z_id) ~ std_normal();
to_vector(z_trl) ~ std_normal();
b ~ normal(1.70475, 0.5);
b_nu ~ normal(2.89037, 0.75);
sd_id[1:4] ~ student_t(3, 0, 0.35);
sd_id[5:8] ~ student_t(3, 0, 1);
sd_trl ~ student_t(3, 0, 0.35);
L_id ~ lkj_corr_cholesky(1);
L_trl ~ lkj_corr_cholesky(1);
// likelihood including all constants
for (n in 1:NOBS) {
  target += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
}
}
generated quantities {
corr_matrix[8] cor_id_full = multiply_lower_tri_self_transpose(L_id);
corr_matrix[4] cor_trl = multiply_lower_tri_self_transpose(L_trl);
matrix[4, 4] cor_id = block(cor_id_full, 1, 1, 4, 4);  // person-mean cross-cell corr
}
