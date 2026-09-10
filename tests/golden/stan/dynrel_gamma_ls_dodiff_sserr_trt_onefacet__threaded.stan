functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
  real partial_log_lik(array[] int seq_slice, int start, int end,
    vector meas, vector mu, vector nu) {
    real lp = 0;
    for (n in start:end) {
      lp += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
    }
    return lp;
  }
}
data {
int<lower=1> NOBS;  // number of observations
int<lower=1> grainsize;  // grainsize for reduce_sum
array[NOBS] int<lower=1,upper=NOBS> seq;
vector<lower=0>[NOBS] meas;  // response (positive support)
array[NOBS] int<lower=1,upper=4> cell;  // constituent index (dodmap order)
int<lower=1> KDIM;  // number of dimension predictors (1 or 3)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) dimension design
int<lower=1> NSUB;  // number of person levels
array[NOBS] int<lower=1> ID;  // person index per observation
int<lower=1> NTRL;  // number of trial levels
array[NOBS] int<lower=1> TRL;  // ordinal trial index per observation
int<lower=1> NOCC;  // number of occasion levels
array[NOBS] int<lower=1> OCC;  // occasion index per observation
int<lower=1> NTID;  // number of person x trial levels
array[NOBS] int<lower=1> TID;  // person x trial index per observation
int<lower=1> NOID;  // number of person x occasion levels
array[NOBS] int<lower=1> OID;  // person x occasion index per observation
int<lower=1> NTO;  // number of trial x occasion levels
array[NOBS] int<lower=1> TO;  // trial x occasion index per observation
vector[4] mu_offset;  // data-derived per-cell log-mean centre
vector[4] log_sigma_offset;  // data-derived per-cell log residual-SD centre
}
transformed data {
matrix[NOBS, KDIM] Xdim = to_matrix(Xdim_vec, NOBS, KDIM);  //column-major
}
parameters {
vector[4] b;  // per-cell log-mean intercept DEVIATIONS from mu_offset
vector[4] b_sigma;  // per-cell log residual-SD DEVIATIONS from log_sigma_offset
matrix[4, KDIM] b_dim;  // per-cell log-mean dimension slopes
matrix[4, KDIM] b_sigma_dim;  // per-cell log residual-SD dimension slopes
vector<lower=0>[8] sd_id;  // person SDs: [loc_1..4, logsd_1..4]
matrix[8, NSUB] z_id;  // standardized person effects
cholesky_factor_corr[8] L_id;  // joint person correlation (locations + log-sigmas)
vector<lower=0>[4] sd_trl;  // per-cell trial main-effect log-SDs
matrix[4, NTRL] z_trl;  // standardized trial effects
cholesky_factor_corr[4] L_trl;  // cross-cell trial correlation
vector<lower=0>[4] sd_occ;  // per-cell occasion main-effect log-SDs
matrix[4, NOCC] z_occ;  // standardized occasion effects
cholesky_factor_corr[4] L_occ;  // cross-cell occasion correlation
vector<lower=0>[4] sd_tid;  // per-cell person x trial log-SDs
matrix[4, NTID] z_tid;  // standardized person x trial effects
cholesky_factor_corr[4] L_tid;  // cross-cell person x trial correlation
vector<lower=0>[4] sd_oid;  // per-cell person x occasion log-SDs
matrix[4, NOID] z_oid;  // standardized person x occasion effects
cholesky_factor_corr[4] L_oid;  // cross-cell person x occasion correlation
vector<lower=0>[4] sd_to;  // per-cell trial x occasion log-SDs
matrix[4, NTO] z_to;  // standardized trial x occasion effects
cholesky_factor_corr[4] L_to;  // cross-cell trial x occasion correlation
}
transformed parameters {
matrix[NSUB, 8] r_id;
matrix[NTRL, 4] r_trl;
matrix[NOCC, 4] r_occ;
matrix[NTID, 4] r_tid;
matrix[NOID, 4] r_oid;
matrix[NTO, 4] r_to;
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
r_occ = (diag_pre_multiply(sd_occ, L_occ) * z_occ)';
r_tid = (diag_pre_multiply(sd_tid, L_tid) * z_tid)';
r_oid = (diag_pre_multiply(sd_oid, L_oid) * z_oid)';
r_to = (diag_pre_multiply(sd_to, L_to) * z_to)';
}
model {
// dimension-conditioned, log-linked per-cell mean and residual SD.
// The offsets centre both submodels on the observed response scale.
vector[NOBS] mu;
vector[NOBS] sigma;
vector[NOBS] nu;  // derived below; declared here so all declarations precede statements
for (n in 1:NOBS) {
  int c = cell[n];
  // person + every crossed facet on the mean; person scale effect ONLY on
  // the residual. No facet ever touches the residual SD.
  mu[n] = mu_offset[c] + b[c] + Xdim[n] * b_dim[c]' + r_id[ID[n], c] + r_trl[TRL[n], c]
    + r_occ[OCC[n], c] + r_tid[TID[n], c] + r_oid[OID[n], c] + r_to[TO[n], c];
  sigma[n] = log_sigma_offset[c] + b_sigma[c] + Xdim[n] * b_sigma_dim[c]' + r_id[ID[n], 4 + c];
}
mu = exp(mu);
sigma = exp(sigma);
// moment-match the scaled chi-square to (mean, SD): Var(Y | mu, sigma) = sigma^2
nu = 2 * square(mu ./ sigma);
// priors including all constants
target += normal_lpdf(b | 0, 0.75);
target += normal_lpdf(b_sigma | 0, 0.5);
target += normal_lpdf(to_vector(b_dim) | 0, 0.5);
target += normal_lpdf(to_vector(b_sigma_dim) | 0, 0.35);
target += student_t_lpdf(sd_id[1:4] | 3, 0, 0.5)
- 4 * student_t_lccdf(0 | 3, 0, 10);
target += student_t_lpdf(sd_id[5:8] | 3, 0, 0.35)
- 4 * student_t_lccdf(0 | 3, 0, 2.5);
target += std_normal_lpdf(to_vector(z_id));
target += lkj_corr_cholesky_lpdf(L_id | 2);
target += student_t_lpdf(sd_trl | 3, 0, 0.5)
- 4 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_trl));
target += lkj_corr_cholesky_lpdf(L_trl | 2);
target += student_t_lpdf(sd_occ | 3, 0, 0.5)
- 4 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_occ));
target += lkj_corr_cholesky_lpdf(L_occ | 2);
target += student_t_lpdf(sd_tid | 3, 0, 0.5)
- 4 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_tid));
target += lkj_corr_cholesky_lpdf(L_tid | 2);
target += student_t_lpdf(sd_oid | 3, 0, 0.5)
- 4 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_oid));
target += lkj_corr_cholesky_lpdf(L_oid | 2);
target += student_t_lpdf(sd_to | 3, 0, 0.5)
- 4 * student_t_lccdf(0 | 3, 0, 10);
target += std_normal_lpdf(to_vector(z_to));
target += lkj_corr_cholesky_lpdf(L_to | 2);
// likelihood including all constants
target += reduce_sum(partial_log_lik, seq, grainsize,
  meas, mu, nu);
}
generated quantities {
corr_matrix[8] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[4] cor_trl = multiply_lower_tri_self_transpose(L_trl);
corr_matrix[4] cor_occ = multiply_lower_tri_self_transpose(L_occ);
corr_matrix[4] cor_tid = multiply_lower_tri_self_transpose(L_tid);
corr_matrix[4] cor_oid = multiply_lower_tri_self_transpose(L_oid);
corr_matrix[4] cor_to = multiply_lower_tri_self_transpose(L_to);
}
