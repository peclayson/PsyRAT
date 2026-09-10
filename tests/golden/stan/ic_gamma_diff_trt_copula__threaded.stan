functions {
  real biv_scaled_chisq_copula_log_density(real y1, real y2,
    real mu1, real mu2, real nu1, real nu2, real rho) {
    if (y1 <= 0 || y2 <= 0) {
      return negative_infinity();
    }
    real shape1 = 0.5 * nu1;
    real shape2 = 0.5 * nu2;
    real rate1 = 0.5 * nu1 / mu1;
    real rate2 = 0.5 * nu2 / mu2;
    real u1 = gamma_cdf(y1 | shape1, rate1);
    real u2 = gamma_cdf(y2 | shape2, rate2);
    real z1;
    real z2;
    real one_minus_r2 = 1 - square(rho);
    real log_c;
    u1 = fmin(1 - 1e-10, fmax(1e-10, u1));
    u2 = fmin(1 - 1e-10, fmax(1e-10, u2));
    z1 = inv_Phi(u1);
    z2 = inv_Phi(u2);
    log_c = -0.5 * log(one_minus_r2)
      - 0.5 * (square(rho) * (square(z1) + square(z2)) - 2 * rho * z1 * z2) / one_minus_r2;
    return gamma_lpdf(y1 | shape1, rate1) + gamma_lpdf(y2 | shape2, rate2) + log_c;
  }
  real partial_log_lik(array[] int seq_slice, int start, int end,
    vector meas1, vector meas2, vector mu1, vector mu2,
    vector nu1, vector nu2, real rho_e) {
    real lp = 0;
    for (n in start:end) {
      lp += biv_scaled_chisq_copula_log_density(meas1[n], meas2[n],
        mu1[n], mu2[n], nu1[n], nu2[n], rho_e);
    }
    return lp;
  }
}
data {
int<lower=1> NOBS;  // number of paired observations
int<lower=1> grainsize;  // grainsize for reduce_sum
array[NOBS] int<lower=1,upper=NOBS> seq;
vector[NOBS] meas1;  // response variable for event 1
vector[NOBS] meas2;  // response variable for event 2
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
real<lower=-0.95, upper=0.95> rho_e;  // Gaussian-copula residual correlation
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
vector[NOBS] mu1;
vector[NOBS] mu2;
vector[NOBS] nu1;
vector[NOBS] nu2;
for (n in 1:NOBS) {
  mu1[n] = b[1] + r_id[ID[n], 1];
  mu2[n] = b[2] + r_id[ID[n], 2];
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
  nu1[n] = b_nu[1] + r_id[ID[n], 3];
  nu2[n] = b_nu[2] + r_id[ID[n], 4];
}
mu1 = exp(mu1);
mu2 = exp(mu2);
nu1 = exp(nu1);
nu2 = exp(nu2);
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
sd_id[3] ~ student_t(3, 0, 0.5);
sd_id[4] ~ student_t(3, 0, 0.5);
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
rho_e ~ normal(0, 0.5);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize, meas1, meas2, mu1, mu2, nu1, nu2, rho_e);
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
