functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
}
data {
int<lower=1> NOBS;  // number of paired observations
vector<lower=0>[NOBS] meas1;  // event-1 response (positive support)
vector<lower=0>[NOBS] meas2;  // event-2 response (positive support)
int<lower=1> KDIM;  // number of dimension predictors (1 or 3)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) dimension design
int<lower=1> NSUB;  // number of person levels
array[NOBS] int<lower=1> ID;  // person index per observation
int<lower=1> NTRL;  // number of trial levels
array[NOBS] int<lower=1> TRL;  // trial index per observation
int<lower=1> NOCC;  // number of occasion levels
array[NOBS] int<lower=1> OCC;  // occasion index per observation
int<lower=1> NTID;  // number of trial x person levels
array[NOBS] int<lower=1> TID;  // trial x person index per observation
int<lower=1> NOID;  // number of occasion x person levels
array[NOBS] int<lower=1> OID;  // occasion x person index per observation
int<lower=1> NTO;  // number of trial x occasion levels
array[NOBS] int<lower=1> TO;  // trial x occasion index per observation
vector[2] mu_offset;  // data-derived per-event log-mean centre
vector[2] log_sigma_offset;  // data-derived per-event log residual-SD centre
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
vector[2] b;  // per-event log-mean intercept DEVIATIONS from mu_offset
vector[2] b_sigma;  // per-event log residual-SD DEVIATIONS from log_sigma_offset
vector[KDIM] b_dim1;  // event-1 log-mean dimension slopes
vector[KDIM] b_dim2;  // event-2 log-mean dimension slopes
vector[KDIM] b_sigma_dim1;  // event-1 log residual-SD dimension slopes
vector[KDIM] b_sigma_dim2;  // event-2 log residual-SD dimension slopes
vector<lower=0>[4] sd_id;  // person SDs: [mean_1, mean_2, logsigma_1, logsigma_2]
matrix[4, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[4] L_id;  // person correlation (means + log-sigmas)
vector<lower=0>[2] sd_trl;  // trial log-SD (per event)
matrix[2, NTRL] z_trl;  // trial standardized effects
cholesky_factor_corr[2] L_trl;  // trial cross-event corr
vector<lower=0>[2] sd_occ;  // occasion log-SD (per event)
matrix[2, NOCC] z_occ;  // occasion standardized effects
cholesky_factor_corr[2] L_occ;  // occasion cross-event corr
vector<lower=0>[2] sd_tid;  // trial x person log-SD (per event)
matrix[2, NTID] z_tid;  // trial x person standardized effects
cholesky_factor_corr[2] L_tid;  // trial x person cross-event corr
vector<lower=0>[2] sd_oid;  // occasion x person log-SD (per event)
matrix[2, NOID] z_oid;  // occasion x person standardized effects
cholesky_factor_corr[2] L_oid;  // occasion x person cross-event corr
vector<lower=0>[2] sd_to;  // trial x occasion log-SD (per event)
matrix[2, NTO] z_to;  // trial x occasion standardized effects
cholesky_factor_corr[2] L_to;  // trial x occasion cross-event corr
}
transformed parameters {
matrix[NSUB, 4] r_id;
matrix[NTRL, 2] r_trl;
matrix[NOCC, 2] r_occ;
matrix[NTID, 2] r_tid;
matrix[NOID, 2] r_oid;
matrix[NTO, 2] r_to;
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
r_occ = (diag_pre_multiply(sd_occ, L_occ) * z_occ)';
r_tid = (diag_pre_multiply(sd_tid, L_tid) * z_tid)';
r_oid = (diag_pre_multiply(sd_oid, L_oid) * z_oid)';
r_to = (diag_pre_multiply(sd_to, L_to) * z_to)';
}
model {
// dimension-conditioned, log-linked per-event mean and residual SD.
// The offsets centre both submodels on the observed response scale.
vector[NOBS] mu1 = rep_vector(mu_offset[1] + b[1], NOBS) + Xdim * b_dim1;
vector[NOBS] mu2 = rep_vector(mu_offset[2] + b[2], NOBS) + Xdim * b_dim2;
vector[NOBS] sigma1 = rep_vector(log_sigma_offset[1] + b_sigma[1], NOBS) + Xdim * b_sigma_dim1;
vector[NOBS] sigma2 = rep_vector(log_sigma_offset[2] + b_sigma[2], NOBS) + Xdim * b_sigma_dim2;
vector[NOBS] nu1;  // derived below; declared here so all declarations precede statements
vector[NOBS] nu2;
for (n in 1:NOBS) {
  // person location effect on the mean (the five crossed facets follow)
  mu1[n] += r_id[ID[n], 1];
  mu2[n] += r_id[ID[n], 2];
  mu1[n] += r_trl[TRL[n], 1];
  mu2[n] += r_trl[TRL[n], 2];
  mu1[n] += r_occ[OCC[n], 1];
  mu2[n] += r_occ[OCC[n], 2];
  mu1[n] += r_tid[TID[n], 1];
  mu2[n] += r_tid[TID[n], 2];
  mu1[n] += r_oid[OID[n], 1];
  mu2[n] += r_oid[OID[n], 2];
  mu1[n] += r_to[TO[n], 1];
  mu2[n] += r_to[TO[n], 2];
  // person scale effect on the per-event residual; no facet touches it
  sigma1[n] += r_id[ID[n], 3];
  sigma2[n] += r_id[ID[n], 4];
}
mu1 = exp(mu1);
mu2 = exp(mu2);
sigma1 = exp(sigma1);
sigma2 = exp(sigma2);
// moment-match the scaled chi-square to (mean, SD): Var(Y | mu, sigma) = sigma^2
nu1 = 2 * square(mu1 ./ sigma1);
nu2 = 2 * square(mu2 ./ sigma2);
// priors including all constants
target += normal_lpdf(b | 0, 0.75);
target += normal_lpdf(b_sigma | 0, 0.5);
target += student_t_lpdf(sd_id[1:2] | 3, 0, 0.5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += student_t_lpdf(sd_id[3:4] | 3, 0, 0.35)
- 2 * student_t_lccdf(0 | 3, 0, 2.5);
target += std_normal_lpdf(to_vector(z_id));
target += lkj_corr_cholesky_lpdf(L_id | 2);
target += student_t_lpdf(sd_trl | 3, 0, 0.5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_trl));
target += lkj_corr_cholesky_lpdf(L_trl | 2);
target += student_t_lpdf(sd_occ | 3, 0, 0.5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_occ));
target += lkj_corr_cholesky_lpdf(L_occ | 2);
target += student_t_lpdf(sd_tid | 3, 0, 0.5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_tid));
target += lkj_corr_cholesky_lpdf(L_tid | 2);
target += student_t_lpdf(sd_oid | 3, 0, 0.5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_oid));
target += lkj_corr_cholesky_lpdf(L_oid | 2);
target += student_t_lpdf(sd_to | 3, 0, 0.5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_to));
target += lkj_corr_cholesky_lpdf(L_to | 2);
target += normal_lpdf(b_dim1 | 0, 0.5);
target += normal_lpdf(b_dim2 | 0, 0.5);
target += normal_lpdf(b_sigma_dim1 | 0, 0.35);
target += normal_lpdf(b_sigma_dim2 | 0, 0.35);
// likelihood including all constants
for (n in 1:NOBS) {
  target += scaled_chi_square_lpdf(meas1[n] | mu1[n], nu1[n]);
  target += scaled_chi_square_lpdf(meas2[n] | mu2[n], nu2[n]);
}
}
generated quantities {
corr_matrix[4] cor_id = multiply_lower_tri_self_transpose(L_id);
real cor_p = cor_id[1, 2];  // cross-event person-mean correlation
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
real cor_t = cor_trl[1, 2];
corr_matrix[2] cor_occ = multiply_lower_tri_self_transpose(L_occ);
real cor_o = cor_occ[1, 2];
corr_matrix[2] cor_tid = multiply_lower_tri_self_transpose(L_tid);
real cor_pt = cor_tid[1, 2];
corr_matrix[2] cor_oid = multiply_lower_tri_self_transpose(L_oid);
real cor_po = cor_oid[1, 2];
corr_matrix[2] cor_to = multiply_lower_tri_self_transpose(L_to);
real cor_ot = cor_to[1, 2];
}
