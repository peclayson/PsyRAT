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
vector<lower=0>[NOBS] meas;  // response (positive support)
array[NOBS] int<lower=1,upper=4> cell;  // constituent index (dodmap order)
int<lower=1> KDIM;  // number of dimension predictors (1 or 3)
vector[NOBS * KDIM] Xdim_vec;  // flattened (column-major) dimension design
int<lower=1> NSUB;  // number of person levels
array[NOBS] int<lower=1> ID;  // person index per observation
int<lower=1> NTRL;  // number of trial levels
array[NOBS] int<lower=1> TRL;  // ordinal trial index per observation
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
}
transformed parameters {
matrix[NSUB, 8] r_id;
matrix[NTRL, 4] r_trl;
r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
}
model {
// dimension-conditioned, log-linked per-cell mean and residual SD.
// The offsets centre both submodels on the observed response scale.
vector[NOBS] mu;
vector[NOBS] sigma;
vector[NOBS] nu;  // derived below; declared here so all declarations precede statements
for (n in 1:NOBS) {
  int c = cell[n];
  // person location + trial effects on the mean; person scale effect on the
  // residual. The trial factor never touches the residual SD.
  mu[n] = mu_offset[c] + b[c] + Xdim[n] * b_dim[c]' + r_id[ID[n], c] + r_trl[TRL[n], c];
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
// likelihood including all constants
for (n in 1:NOBS) {
  target += scaled_chi_square_lpdf(meas[n] | mu[n], nu[n]);
}
}
generated quantities {
corr_matrix[8] cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[4] cor_trl = multiply_lower_tri_self_transpose(L_trl);
}
