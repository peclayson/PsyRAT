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
matrix[NOBS, 2] X;  // population-level (cell) design matrix
int<lower=1> KDIM;  // number of event-by-dimension predictors (2 or 6)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) event-by-dimension design
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
vector[NOBS] meas1;  // condition-1 indicator (1,0)
vector[NOBS] meas2;  // condition-2 indicator (1,0)
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
vector[2] b;  // population-level cell means
vector[2] b_sigma;  // population-level residual log-SD
vector[KDIM] b_dim;  // event-specific mean dimension slopes
vector[KDIM] b_sigma_dim;  // event-specific log-residual dimension slopes
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
}
transformed parameters {
matrix[NSUB, 4] r_1;
vector[NSUB] r_1_1;
vector[NSUB] r_1_2;
vector[NSUB] r_1_s1;
vector[NSUB] r_1_s2;
r_1 = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_1_1 = r_1[, 1];
r_1_2 = r_1[, 2];
r_1_s1 = r_1[, 3];
r_1_s2 = r_1[, 4];
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
vector[NOBS] mu = X * b + Xdim * b_dim;
vector[NOBS] sig_err = X * b_sigma + Xdim * b_sigma_dim;
for (n in 1:NOBS) {
  // add the person + crossed cross-condition location effects
  mu[n] += r_1_1[ID[n]] * meas1[n] + r_1_2[ID[n]] * meas2[n];
  mu[n] += r_trl_1[TRL[n]] * meas1[n] + r_trl_2[TRL[n]] * meas2[n];
  mu[n] += r_occ_1[OCC[n]] * meas1[n] + r_occ_2[OCC[n]] * meas2[n];
  mu[n] += r_tid_1[TID[n]] * meas1[n] + r_tid_2[TID[n]] * meas2[n];
  mu[n] += r_oid_1[OID[n]] * meas1[n] + r_oid_2[OID[n]] * meas2[n];
  mu[n] += r_to_1[TO[n]] * meas1[n] + r_to_2[TO[n]] * meas2[n];
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
target += normal_lpdf(b_dim | 0, 10);
target += normal_lpdf(b_sigma_dim | 0, 1);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize,
  meas, mu, sig_err);
}
generated quantities {
matrix[2,2] id_varcov;
matrix[NSUB,2] er_var_ss;
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
// joint 4x4 person block: location (1:2) + scale (3:4) cholesky
corr_matrix[4] cor_id = multiply_lower_tri_self_transpose(L_id);
matrix[4,4] id_varcov_full = quad_form_diag(cor_id,sd_id);
id_varcov = id_varcov_full[1:2, 1:2];
for (s in 1:NSUB) {
  er_var_ss[s,1] = b_sigma[1] + r_1_s1[s];
  er_var_ss[s,2] = b_sigma[2] + r_1_s2[s];
}
}
