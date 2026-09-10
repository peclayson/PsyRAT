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
vector[2] b;  // population-level effects
vector[2] b_sigma;  // population-level effects
vector<lower=0>[2] sd_id;  // group-level standard deviations
matrix[2, NSUB] z_id;  // standardized group-level effects
cholesky_factor_corr[2] L_id;  // cholesky factor of correlation matrix
vector<lower=0>[2] sd_trl;  // group-level standard deviations
matrix[2, NTRL] z_trl;  // standardized group-level effects
cholesky_factor_corr[2] L_trl;  // cholesky factor of correlation matrix
}
transformed parameters {
matrix[NSUB, 2] r_1;
vector[NSUB] r_1_1;
vector[NSUB] r_1_2;
matrix[NTRL, 2] r_2;
vector[NTRL] r_2_1;
vector[NTRL] r_2_2;
r_1 = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_1_1 = r_1[, 1];
r_1_2 = r_1[, 2];
r_2 = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
r_2_1 = r_2[, 1];
r_2_2 = r_2[, 2];
}
model {
// initialize linear predictor term
vector[NOBS] mu = X * b;
vector[NOBS] sig_err = X * b_sigma;
for (n in 1:NOBS) {
  // add more terms to the linear predictor
  mu[n] += r_1_1[ID[n]] * meas1[n] + r_1_2[ID[n]] * meas2[n] + 
  r_2_1[TRL[n]] * meas1[n] + r_2_2[TRL[n]] * meas2[n];
}
for (n in 1:NOBS) {
  // apply the inverse link function
  sig_err[n] = exp(sig_err[n]);
}
// priors including all constants
target += normal_lpdf(b | 0,7);
target += student_t_lpdf(b_sigma | 3, 0, 8);
target += student_t_lpdf(sd_id | 3, 0, 9)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_id));
target += lkj_corr_cholesky_lpdf(L_id | 2);
target += student_t_lpdf(sd_trl | 3, 0, 11)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_trl));
target += lkj_corr_cholesky_lpdf(L_trl | 2);
// likelihood including all constants
target += normal_lpdf(meas | mu, sig_err);
}
generated quantities {
matrix[2,2] id_varcov;
matrix[2,2] trl_varcov;
// compute group-level correlations
corr_matrix[2] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
id_varcov = quad_form_diag(cor_id,sd_id);
trl_varcov = quad_form_diag(cor_trl,sd_trl);
}
