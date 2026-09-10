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
matrix[NOBS, 2] X;  // population-level (cell) design matrix
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
vector[NOBS] meas1;  // condition-1 indicator (1,0)
vector[NOBS] meas2;  // condition-2 indicator (1,0)
}
parameters {
vector[2] b;  // per-event log-mean intercepts
vector[2] b_nu;  // per-event log-nu intercepts
vector<lower=0>[4] sd_id;  // person SDs: [mean_1, mean_2, lognu_1, lognu_2]
matrix[4, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[4] L_id;  // person correlation (means + log-nus)
vector<lower=0>[2] sd_trl;  // trial SD (per condition)
matrix[2, NTRL] z_trl;  // trial standardized effects
cholesky_factor_corr[2] L_trl;  // trial cross-condition corr
vector<lower=0>[2] sd_occ;  // occasion SD (per condition)
matrix[2, NOCC] z_occ;  // occasion standardized effects
cholesky_factor_corr[2] L_occ;  // occasion cross-condition corr
vector<lower=0>[2] sd_tid;  // trial x person SD (per condition)
matrix[2, NTID] z_tid;  // trial x person standardized effects
cholesky_factor_corr[2] L_tid;  // trial x person cross-condition corr
vector<lower=0>[2] sd_oid;  // occasion x person SD (per condition)
matrix[2, NOID] z_oid;  // occasion x person standardized effects
cholesky_factor_corr[2] L_oid;  // occasion x person cross-condition corr
vector<lower=0>[2] sd_to;  // trial x occasion SD (per condition)
matrix[2, NTO] z_to;  // trial x occasion standardized effects
cholesky_factor_corr[2] L_to;  // trial x occasion cross-condition corr
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
vector[NOBS] mu = X * b;
vector[NOBS] nu = X * b_nu;
for (n in 1:NOBS) {
  mu[n] += r_id[ID[n], 1] * meas1[n] + r_id[ID[n], 2] * meas2[n];
  mu[n] += r_trl[TRL[n], 1] * meas1[n] + r_trl[TRL[n], 2] * meas2[n];
  mu[n] += r_occ[OCC[n], 1] * meas1[n] + r_occ[OCC[n], 2] * meas2[n];
  mu[n] += r_tid[TID[n], 1] * meas1[n] + r_tid[TID[n], 2] * meas2[n];
  mu[n] += r_oid[OID[n], 1] * meas1[n] + r_oid[OID[n], 2] * meas2[n];
  mu[n] += r_to[TO[n], 1] * meas1[n] + r_to[TO[n], 2] * meas2[n];
  nu[n] += r_id[ID[n], 3] * meas1[n] + r_id[ID[n], 4] * meas2[n];
}
mu = exp(mu);
nu = exp(nu);
// priors including all constants
to_vector(z_id) ~ std_normal();
to_vector(z_trl) ~ std_normal();
to_vector(z_occ) ~ std_normal();
to_vector(z_tid) ~ std_normal();
to_vector(z_oid) ~ std_normal();
to_vector(z_to) ~ std_normal();
b ~ normal(1.70475, 0.5);
b_nu ~ normal(2.89037, 0.75);
sd_id[1] ~ student_t(3, 0, 0.35);
sd_id[2] ~ student_t(3, 0, 0.35);
sd_id[3] ~ student_t(3, 0, 1);
sd_id[4] ~ student_t(3, 0, 1);
L_id ~ lkj_corr_cholesky(1);
sd_trl ~ student_t(3, 0, 0.35);
L_trl ~ lkj_corr_cholesky(1);
sd_occ ~ student_t(3, 0, 0.35);
L_occ ~ lkj_corr_cholesky(1);
sd_tid ~ student_t(3, 0, 0.35);
L_tid ~ lkj_corr_cholesky(1);
sd_oid ~ student_t(3, 0, 0.35);
L_oid ~ lkj_corr_cholesky(1);
sd_to ~ student_t(3, 0, 0.35);
L_to ~ lkj_corr_cholesky(1);
// likelihood including all constants
for (n in 1:NOBS) {
  target += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
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
