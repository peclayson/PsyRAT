functions {
real partial_log_lik(array[] int seq_slice, int start, int end,
  vector meas, vector mu, vector sig_err) {
  return normal_lpdf(meas[start:end] | mu[start:end], sig_err[start:end]);
}
}
data {
int<lower=1> NOBS;  // number of observations
int<lower=1> grainsize;  // grainsize for reduce_sum
array[NOBS] int<lower=1,upper=NOBS> seq;
vector[NOBS] meas;  // response variable
matrix[NOBS, 2] X;  // population-level design matrix
int<lower=1> KDIM;  // number of event-by-dimension predictors (2 or 6)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) event-by-dimension design
int<lower=1> NSUB;  // number of grouping levels
array[NOBS] int<lower=1> ID;  // grouping indicator per observation
int<lower=1> NTRL;  // number of grouping levels
array[NOBS] int<lower=1> TRL;  // grouping indicator per observation
vector[NOBS] meas1;
vector[NOBS] meas2;
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
vector[2] b;  // population-level effects
vector[2] b_sigma;  // population-level effects
vector[KDIM] b_dim;  // event-specific mean dimension slopes
vector[KDIM] b_sigma_dim;  // event-specific log-residual dimension slopes
vector<lower=0>[4] sd_id;  // joint person block SDs (1-2 location, 3-4 scale)
matrix[4, NSUB] z_id;  // standardized group-level effects
cholesky_factor_corr[4] L_id;  // cholesky factor of correlation matrix
vector<lower=0>[2] sd_trl;  // group-level standard deviations
matrix[2, NTRL] z_trl;  // standardized group-level effects
cholesky_factor_corr[2] L_trl;  // cholesky factor of correlation matrix
}
transformed parameters {
matrix[NSUB, 4] r_1;
vector[NSUB] r_1_1;
vector[NSUB] r_1_2;
vector[NSUB] r_1_s1;
vector[NSUB] r_1_s2;
matrix[NTRL, 2] r_2;
vector[NTRL] r_2_1;
vector[NTRL] r_2_2;
r_1 = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_1_1 = r_1[, 1];
r_1_2 = r_1[, 2];
r_1_s1 = r_1[, 3];
r_1_s2 = r_1[, 4];
r_2 = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
r_2_1 = r_2[, 1];
r_2_2 = r_2[, 2];
}
model {
// dimension-conditioned linear predictors (event-specific fixed slopes)
vector[NOBS] mu = X * b + Xdim * b_dim;
vector[NOBS] sig_err = X * b_sigma + Xdim * b_sigma_dim;
for (n in 1:NOBS) {
  // add more terms to the linear predictor
  mu[n] += r_1_1[ID[n]] * meas1[n] + r_1_2[ID[n]] * meas2[n] + 
  r_2_1[TRL[n]] * meas1[n] + r_2_2[TRL[n]] * meas2[n];
}
for (n in 1:NOBS) {
  // person event-specific scale random effect (subject-level residual)
  sig_err[n] += r_1_s1[ID[n]] * meas1[n] + r_1_s2[ID[n]] * meas2[n];
}
for (n in 1:NOBS) {
  // apply the inverse link function
  sig_err[n] = exp(sig_err[n]);
}
// priors including all constants
target += normal_lpdf(b | 0,5);
target += student_t_lpdf(b_sigma | 3, 0, 10);
target += student_t_lpdf(sd_id[1:2] | 3, 0, 10)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += student_t_lpdf(sd_id[3:4] | 3, 0, 2.5)
- 2 * student_t_lccdf(0 | 3, 0, 2.5);
target += std_normal_lpdf(to_vector(z_id));
target += lkj_corr_cholesky_lpdf(L_id | 2);
target += student_t_lpdf(sd_trl | 3, 0, 10)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_trl));
target += lkj_corr_cholesky_lpdf(L_trl | 2);
target += normal_lpdf(b_dim | 0, 10);
target += normal_lpdf(b_sigma_dim | 0, 1);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize,
  meas, mu, sig_err);
}
generated quantities {
matrix[2,2] id_varcov;
matrix[2,2] trl_varcov;
matrix[NSUB,2] er_var_ss;
// compute group-level correlations
corr_matrix[4] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
matrix[4,4] id_varcov_full = quad_form_diag(cor_id,sd_id);
id_varcov = id_varcov_full[1:2, 1:2];
trl_varcov = quad_form_diag(cor_trl,sd_trl);
for (s in 1:NSUB) {
  er_var_ss[s,1] = b_sigma[1] + r_1_s1[s];
  er_var_ss[s,2] = b_sigma[2] + r_1_s2[s];
}
}
