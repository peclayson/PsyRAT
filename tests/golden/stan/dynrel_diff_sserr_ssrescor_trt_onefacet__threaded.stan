functions {
real partial_log_lik(array[] int seq_slice, int start, int end,
  vector meas1, vector meas2, vector mu1, vector mu2,
  vector sigma1, vector sigma2, vector rho, array[] int ID) {
  real lp = 0;
  for (n in start:end) {
    vector[2] Y = transpose([meas1[n], meas2[n]]);
    vector[2] Mu = transpose([mu1[n], mu2[n]]);
    matrix[2,2] Lcorr;
    matrix[2,2] LSigma;
    Lcorr[1,1] = 1;
    Lcorr[1,2] = 0;
    Lcorr[2,1] = rho[ID[n]];
    Lcorr[2,2] = sqrt(1 - square(rho[ID[n]]));
    LSigma = diag_pre_multiply(transpose([sigma1[n], sigma2[n]]), Lcorr);
    lp += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
  }
  return lp;
}
}
data {
int<lower=1> NOBS;  // number of paired observations
int<lower=1> grainsize;  // grainsize for reduce_sum
array[NOBS] int<lower=1,upper=NOBS> seq;
vector[NOBS] meas1;  // response variable for condition 1
vector[NOBS] meas2;  // response variable for condition 2
int<lower=1> KDIM;  // number of dimension predictors (1 or 3)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) dimension design
int<lower=1> NSUB;  // number of person levels
array[NOBS] int<lower=1> ID;  // person index per observation
int<lower=1> NTRL;  // number of trial sigma(i) levels
array[NOBS] int<lower=1> TRL;  // trial sigma(i) index per observation
int<lower=1> NOCC;  // number of occasion sigma(o) levels
array[NOBS] int<lower=1> OCC;  // occasion sigma(o) index per observation
int<lower=1> NTID;  // number of person x trial sigma(pi) levels
array[NOBS] int<lower=1> TID;  // person x trial sigma(pi) index per observation
int<lower=1> NOID;  // number of person x occasion sigma(po) levels
array[NOBS] int<lower=1> OID;  // person x occasion sigma(po) index per observation
int<lower=1> NTO;  // number of trial x occasion sigma(io) levels
array[NOBS] int<lower=1> TO;  // trial x occasion sigma(io) index per observation
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
vector[2] b;  // population-level condition means
vector[2] b_sigma;  // population-level residual log-SD
vector[KDIM] b_dim1;  // event-1 mean dimension slopes
vector[KDIM] b_dim2;  // event-2 mean dimension slopes
vector[KDIM] b_sigma_dim1;  // event-1 log-residual dimension slopes
vector[KDIM] b_sigma_dim2;  // event-2 log-residual dimension slopes
vector<lower=0>[4] sd_id;  // joint person block SDs (1-2 location, 3-4 scale)
matrix[4, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[4] L_id;  // joint person cholesky (location + scale)
vector<lower=0>[2] sd_trl;  // trial sigma(i) SD (per condition)
matrix[2, NTRL] z_trl;  // trial sigma(i) standardized effects
cholesky_factor_corr[2] L_trl;  // trial sigma(i) cross-condition corr
vector<lower=0>[2] sd_occ;  // occasion sigma(o) SD (per condition)
matrix[2, NOCC] z_occ;  // occasion sigma(o) standardized effects
cholesky_factor_corr[2] L_occ;  // occasion sigma(o) cross-condition corr
vector<lower=0>[2] sd_tid;  // person x trial sigma(pi) SD (per condition)
matrix[2, NTID] z_tid;  // person x trial sigma(pi) standardized effects
cholesky_factor_corr[2] L_tid;  // person x trial sigma(pi) cross-condition corr
vector<lower=0>[2] sd_oid;  // person x occasion sigma(po) SD (per condition)
matrix[2, NOID] z_oid;  // person x occasion sigma(po) standardized effects
cholesky_factor_corr[2] L_oid;  // person x occasion sigma(po) cross-condition corr
vector<lower=0>[2] sd_to;  // trial x occasion sigma(io) SD (per condition)
matrix[2, NTO] z_to;  // trial x occasion sigma(io) standardized effects
cholesky_factor_corr[2] L_to;  // trial x occasion sigma(io) cross-condition corr
real rescor_mu;  // population mean of the (atanh) residual correlation
real<lower=0> sd_rescor;  // between-person SD of the (atanh) residual correlation
vector[NSUB] z_rescor;  // standardized per-person residual-correlation effects
}
transformed parameters {
matrix[NSUB, 4] r_1;
vector[NSUB] r_1_1;
vector[NSUB] r_1_2;
vector[NSUB] r_1_s1;
vector[NSUB] r_1_s2;
vector[NSUB] rho;  // per-person residual correlation
r_1 = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_1_1 = r_1[, 1];
r_1_2 = r_1[, 2];
r_1_s1 = r_1[, 3];
r_1_s2 = r_1[, 4];
for (s in 1:NSUB) rho[s] = tanh(rescor_mu + sd_rescor * z_rescor[s]);
matrix[NTRL, 2] r_trl;
vector[NTRL] r_trl_1;
vector[NTRL] r_trl_2;
matrix[NOCC, 2] r_occ;
vector[NOCC] r_occ_1;
vector[NOCC] r_occ_2;
matrix[NTID, 2] r_tid;
vector[NTID] r_tid_1;
vector[NTID] r_tid_2;
matrix[NOID, 2] r_oid;
vector[NOID] r_oid_1;
vector[NOID] r_oid_2;
matrix[NTO, 2] r_to;
vector[NTO] r_to_1;
vector[NTO] r_to_2;
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
r_trl_1 = r_trl[, 1];
r_trl_2 = r_trl[, 2];
r_occ = (diag_pre_multiply(sd_occ, L_occ) * z_occ)';
r_occ_1 = r_occ[, 1];
r_occ_2 = r_occ[, 2];
r_tid = (diag_pre_multiply(sd_tid, L_tid) * z_tid)';
r_tid_1 = r_tid[, 1];
r_tid_2 = r_tid[, 2];
r_oid = (diag_pre_multiply(sd_oid, L_oid) * z_oid)';
r_oid_1 = r_oid[, 1];
r_oid_2 = r_oid[, 2];
r_to = (diag_pre_multiply(sd_to, L_to) * z_to)';
r_to_1 = r_to[, 1];
r_to_2 = r_to[, 2];
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
  // person location effect on the mean (the five crossed facets follow)
  mu1[n] += r_1_1[ID[n]];
  mu2[n] += r_1_2[ID[n]];
  mu1[n] += r_trl_1[TRL[n]];
  mu1[n] += r_occ_1[OCC[n]];
  mu1[n] += r_tid_1[TID[n]];
  mu1[n] += r_oid_1[OID[n]];
  mu1[n] += r_to_1[TO[n]];
  mu2[n] += r_trl_2[TRL[n]];
  mu2[n] += r_occ_2[OCC[n]];
  mu2[n] += r_tid_2[TID[n]];
  mu2[n] += r_oid_2[OID[n]];
  mu2[n] += r_to_2[TO[n]];
  // person scale random effect on the per-event residual
  sigma1[n] = exp(logsig1[n] + r_1_s1[ID[n]]);
  sigma2[n] = exp(logsig2[n] + r_1_s2[ID[n]]);
}
// priors including all constants
target += normal_lpdf(b | 0,5);
target += student_t_lpdf(b_sigma | 3, 0, 10);
target += student_t_lpdf(sd_trl | 3, 0, 10)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_trl));
target += lkj_corr_cholesky_lpdf(L_trl | 2);
target += student_t_lpdf(sd_occ | 3, 0, 5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_occ));
target += lkj_corr_cholesky_lpdf(L_occ | 2);
target += student_t_lpdf(sd_tid | 3, 0, 5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_tid));
target += lkj_corr_cholesky_lpdf(L_tid | 2);
target += student_t_lpdf(sd_oid | 3, 0, 5)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_oid));
target += lkj_corr_cholesky_lpdf(L_oid | 2);
target += student_t_lpdf(sd_to | 3, 0, 1)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_to));
target += lkj_corr_cholesky_lpdf(L_to | 2);
target += student_t_lpdf(sd_id[1:2] | 3, 0, 10)
- 2 * student_t_lccdf(0 | 3, 0, 10);
target += student_t_lpdf(sd_id[3:4] | 3, 0, 2.5)
- 2 * student_t_lccdf(0 | 3, 0, 2.5);
target += std_normal_lpdf(to_vector(z_id));
target += lkj_corr_cholesky_lpdf(L_id | 2);
target += normal_lpdf(rescor_mu | 0, 1);
target += student_t_lpdf(sd_rescor | 3, 0, 1)
- student_t_lccdf(0 | 3, 0, 1);
target += std_normal_lpdf(z_rescor);
target += normal_lpdf(b_dim1 | 0, 10);
target += normal_lpdf(b_dim2 | 0, 10);
target += normal_lpdf(b_sigma_dim1 | 0, 1);
target += normal_lpdf(b_sigma_dim2 | 0, 1);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize,
  meas1, meas2, mu1, mu2, sigma1, sigma2, rho, ID);
}
generated quantities {
matrix[2,2] id_varcov;
matrix[NSUB,2] er_var_ss;
matrix[2,2] err_varcov;
vector[1] rescor;
matrix[2,2] trl_varcov;
matrix[2,2] occ_varcov;
matrix[2,2] tid_varcov;
matrix[2,2] oid_varcov;
matrix[2,2] to_varcov;
// compute group-level (co)variances
corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
trl_varcov = quad_form_diag(cor_trl,sd_trl);
corr_matrix[2] cor_occ = multiply_lower_tri_self_transpose(L_occ);
occ_varcov = quad_form_diag(cor_occ,sd_occ);
corr_matrix[2] cor_tid = multiply_lower_tri_self_transpose(L_tid);
tid_varcov = quad_form_diag(cor_tid,sd_tid);
corr_matrix[2] cor_oid = multiply_lower_tri_self_transpose(L_oid);
oid_varcov = quad_form_diag(cor_oid,sd_oid);
corr_matrix[2] cor_to = multiply_lower_tri_self_transpose(L_to);
to_varcov = quad_form_diag(cor_to,sd_to);
// typical-person residual covariance at z = 0 using the population-average rho
err_varcov[1,1] = exp(2 * b_sigma[1]);
err_varcov[2,2] = exp(2 * b_sigma[2]);
err_varcov[1,2] = mean(rho) * exp(b_sigma[1]) * exp(b_sigma[2]);
err_varcov[2,1] = err_varcov[1,2];
rescor[1] = mean(rho);
// joint 4x4 person block: location (1:2) + scale (3:4) cholesky
corr_matrix[4] cor_id = multiply_lower_tri_self_transpose(L_id);
matrix[4,4] id_varcov_full = quad_form_diag(cor_id,sd_id);
id_varcov = id_varcov_full[1:2, 1:2];
for (s in 1:NSUB) {
  er_var_ss[s,1] = b_sigma[1] + r_1_s1[s];
  er_var_ss[s,2] = b_sigma[2] + r_1_s2[s];
}
}
