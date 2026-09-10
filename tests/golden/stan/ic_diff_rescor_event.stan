data {
int<lower=1> NOBS;  // number of paired observations
vector[NOBS] meas1;  // response variable for event 1
vector[NOBS] meas2;  // response variable for event 2
int<lower=1> NSUB;  // number of grouping levels
array[NOBS] int<lower=1> ID;  // grouping indicator per observation
int<lower=1> NTRL;  // number of grouping levels
array[NOBS] int<lower=1> TRL;  // grouping indicator per observation
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
cholesky_factor_corr[2] Lrescor;  // residual correlation matrix
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
vector[NOBS] mu1 = rep_vector(b[1], NOBS);
vector[NOBS] mu2 = rep_vector(b[2], NOBS);
vector[2] sigma = exp(b_sigma);
matrix[2,2] LSigma = diag_pre_multiply(sigma, Lrescor);
for (n in 1:NOBS) {
  mu1[n] += r_1_1[ID[n]] + r_2_1[TRL[n]];
  mu2[n] += r_1_2[ID[n]] + r_2_2[TRL[n]];
}
// priors including all constants
target += normal_lpdf(b | 0,5);
target += student_t_lpdf(b_sigma | 3, 0, 10);
target += student_t_lpdf(sd_id | 3, 0, 10)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_id));
target += lkj_corr_cholesky_lpdf(L_id | 2);
target += student_t_lpdf(sd_trl | 3, 0, 10)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_trl));
target += lkj_corr_cholesky_lpdf(L_trl | 2);
target += lkj_corr_cholesky_lpdf(Lrescor | 2);
// likelihood including all constants
for (n in 1:NOBS) {
  vector[2] Y;
  vector[2] Mu;
  Y = transpose([meas1[n], meas2[n]]);
  Mu = transpose([mu1[n], mu2[n]]);
  target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
}
}
generated quantities {
matrix[2,2] id_varcov;
matrix[2,2] trl_varcov;
matrix[2,2] err_varcov;
vector[1] rescor;
// compute group-level correlations
corr_matrix[2] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
corr_matrix[2] cor_res = multiply_lower_tri_self_transpose(Lrescor);
id_varcov = quad_form_diag(cor_id,sd_id);
trl_varcov = quad_form_diag(cor_trl,sd_trl);
err_varcov = quad_form_diag(cor_res,exp(b_sigma));
rescor[1] = cor_res[1,2];
}
