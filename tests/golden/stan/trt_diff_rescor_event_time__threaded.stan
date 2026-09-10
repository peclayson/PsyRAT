functions {
real partial_log_lik(array[] int seq_slice, int start, int end,
  vector meas1, vector meas2, vector mu1, vector mu2, matrix LSigma) {
  real lp = 0;
  for (n in start:end) {
    vector[2] Y;
    vector[2] Mu;
    Y = transpose([meas1[n], meas2[n]]);
    Mu = transpose([mu1[n], mu2[n]]);
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
int<lower=1> NSUB;  // number of person sigma(p) levels
array[NOBS] int<lower=1> ID;  // person sigma(p) index per observation
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
parameters {
vector[2] b;  // population-level condition means
vector[2] b_sigma;  // population-level residual log-SD
vector<lower=0>[2] sd_id;  // person sigma(p) SD (per condition)
matrix[2, NSUB] z_id;  // person sigma(p) standardized effects
cholesky_factor_corr[2] L_id;  // person sigma(p) cross-condition corr
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
cholesky_factor_corr[2] Lrescor;  // residual correlation matrix
}
transformed parameters {
matrix[NSUB, 2] r_id;
vector[NSUB] r_id_1;
vector[NSUB] r_id_2;
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
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_id_1 = r_id[, 1];
r_id_2 = r_id[, 2];
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
vector[NOBS] mu1 = rep_vector(b[1], NOBS);
vector[NOBS] mu2 = rep_vector(b[2], NOBS);
vector[2] sigma = exp(b_sigma);
matrix[2,2] LSigma = diag_pre_multiply(sigma, Lrescor);
for (n in 1:NOBS) {
  mu1[n] += r_id_1[ID[n]];
  mu1[n] += r_trl_1[TRL[n]];
  mu1[n] += r_occ_1[OCC[n]];
  mu1[n] += r_tid_1[TID[n]];
  mu1[n] += r_oid_1[OID[n]];
  mu1[n] += r_to_1[TO[n]];
  mu2[n] += r_id_2[ID[n]];
  mu2[n] += r_trl_2[TRL[n]];
  mu2[n] += r_occ_2[OCC[n]];
  mu2[n] += r_tid_2[TID[n]];
  mu2[n] += r_oid_2[OID[n]];
  mu2[n] += r_to_2[TO[n]];
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
target += lkj_corr_cholesky_lpdf(Lrescor | 2);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize,
  meas1, meas2, mu1, mu2, LSigma);
}
generated quantities {
matrix[2,2] id_varcov;
matrix[2,2] trl_varcov;
matrix[2,2] occ_varcov;
matrix[2,2] tid_varcov;
matrix[2,2] oid_varcov;
matrix[2,2] to_varcov;
matrix[2,2] err_varcov;
vector[1] rescor;
// compute group-level (co)variances
corr_matrix[2] cor_id = multiply_lower_tri_self_transpose(L_id);
id_varcov = quad_form_diag(cor_id,sd_id);
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
corr_matrix[2] cor_res = multiply_lower_tri_self_transpose(Lrescor);
err_varcov = quad_form_diag(cor_res,exp(b_sigma));
rescor[1] = cor_res[1,2];
}
