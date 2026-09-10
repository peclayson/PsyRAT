data {
int<lower=1> NOBS;  // number of paired observations
vector[NOBS] meas1;  // response variable for event 1
vector[NOBS] meas2;  // response variable for event 2
int<lower=1> KDIM;  // number of dimension predictors (1 or 3)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) dimension design
int<lower=1> NSUB;  // number of grouping levels
array[NOBS] int<lower=1> ID;  // grouping indicator per observation
int<lower=1> NTRL;  // number of grouping levels
array[NOBS] int<lower=1> TRL;  // grouping indicator per observation
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
vector[2] b;  // population-level effects
vector[2] b_sigma;  // population-level effects
vector[KDIM] b_dim1;  // event-1 mean dimension slopes
vector[KDIM] b_dim2;  // event-2 mean dimension slopes
vector[KDIM] b_sigma_dim1;  // event-1 log-residual dimension slopes
vector[KDIM] b_sigma_dim2;  // event-2 log-residual dimension slopes
vector<lower=0>[4] sd_id;  // joint person block SDs (1-2 location, 3-4 scale)
matrix[4, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[4] L_id;  // joint person cholesky (location + scale)
vector<lower=0>[2] sd_trl;  // group-level standard deviations
matrix[2, NTRL] z_trl;  // standardized group-level effects
cholesky_factor_corr[2] L_trl;  // cholesky factor of correlation matrix
cholesky_factor_corr[2] Lrescor;  // residual correlation matrix
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
vector[NOBS] mu1 = rep_vector(b[1], NOBS) + Xdim * b_dim1;
vector[NOBS] mu2 = rep_vector(b[2], NOBS) + Xdim * b_dim2;
vector[NOBS] logsig1 = rep_vector(b_sigma[1], NOBS) + Xdim * b_sigma_dim1;
vector[NOBS] logsig2 = rep_vector(b_sigma[2], NOBS) + Xdim * b_sigma_dim2;
vector[NOBS] sigma1;
vector[NOBS] sigma2;
for (n in 1:NOBS) {
  // person location + trial REs on the mean; person scale RE on the residual
  mu1[n] += r_1_1[ID[n]] + r_2_1[TRL[n]];
  mu2[n] += r_1_2[ID[n]] + r_2_2[TRL[n]];
  sigma1[n] = exp(logsig1[n] + r_1_s1[ID[n]]);
  sigma2[n] = exp(logsig2[n] + r_1_s2[ID[n]]);
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
target += lkj_corr_cholesky_lpdf(Lrescor | 2);
target += normal_lpdf(b_dim1 | 0, 10);
target += normal_lpdf(b_dim2 | 0, 10);
target += normal_lpdf(b_sigma_dim1 | 0, 1);
target += normal_lpdf(b_sigma_dim2 | 0, 1);
// likelihood including all constants
for (n in 1:NOBS) {
  vector[2] Y = transpose([meas1[n], meas2[n]]);
  vector[2] Mu = transpose([mu1[n], mu2[n]]);
  matrix[2,2] LSigma = diag_pre_multiply(transpose([sigma1[n], sigma2[n]]), Lrescor);
  target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
}
}
generated quantities {
matrix[2,2] id_varcov;
matrix[2,2] trl_varcov;
matrix[2,2] err_varcov;
vector[1] rescor;
matrix[NSUB,2] er_var_ss;
// joint 4x4 person block: location (1:2) + scale (3:4) cholesky
corr_matrix[4] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
corr_matrix[2] cor_res = multiply_lower_tri_self_transpose(Lrescor);
matrix[4,4] id_varcov_full = quad_form_diag(cor_id,sd_id);
id_varcov = id_varcov_full[1:2, 1:2];
trl_varcov = quad_form_diag(cor_trl,sd_trl);
// residual covariance at z = 0 (typical person, scale RE = 0)
err_varcov = quad_form_diag(cor_res,exp(b_sigma));
rescor[1] = cor_res[1,2];
for (s in 1:NSUB) {
  er_var_ss[s,1] = b_sigma[1] + r_1_s1[s];
  er_var_ss[s,2] = b_sigma[2] + r_1_s2[s];
}
}
