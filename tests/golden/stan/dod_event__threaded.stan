functions {
real partial_log_lik(array[] int seq_slice, int start, int end,
  vector Y, vector mu, vector log_sigma) {
  return normal_lpdf(Y[start:end] | mu[start:end], exp(log_sigma[start:end]));
}
}
data {
int<lower=1> N;
int<lower=1> grainsize;
array[N] int<lower=1, upper=N> seq;
vector[N] Y;
int<lower=4, upper=4> C;
array[N] int<lower=1, upper=C> cell;
int<lower=1> N_id;
array[N] int<lower=1, upper=N_id> id;
int<lower=1> N_trial;
array[N] int<lower=1, upper=N_trial> trial;
int<lower=0, upper=1> prior_only;
}
parameters {
vector[C] b_cell;
vector[C] b_sigma_cell;
vector<lower=0>[2 * C] sd_id;
cholesky_factor_corr[2 * C] L_id;
matrix[2 * C, N_id] z_id;
vector<lower=0>[2 * C] sd_trial;
cholesky_factor_corr[2 * C] L_trial;
matrix[2 * C, N_trial] z_trial;
}
transformed parameters {
matrix[N_id, 2 * C] r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
matrix[N_trial, 2 * C] r_trial = (diag_pre_multiply(sd_trial, L_trial) * z_trial)';
}
model {
b_cell ~ normal(0, 10);
b_sigma_cell ~ student_t(10, 0, 2);
sd_id ~ student_t(10, 0, 2);
sd_trial ~ student_t(10, 0, 2);
L_id ~ lkj_corr_cholesky(2);
L_trial ~ lkj_corr_cholesky(2);
to_vector(z_id) ~ std_normal();
to_vector(z_trial) ~ std_normal();
if (!prior_only) {
  vector[N] mu;
  vector[N] log_sigma;
  for (n in 1:N) {
    int c = cell[n];
    int c_sigma = c + C;
    mu[n] = b_cell[c] + r_id[id[n], c] + r_trial[trial[n], c];
    log_sigma[n] = b_sigma_cell[c] + r_id[id[n], c_sigma] + r_trial[trial[n], c_sigma];
  }
  target += reduce_sum(partial_log_lik, seq, grainsize, Y, mu, log_sigma);
}
}
generated quantities {
corr_matrix[2 * C] Cor_id = multiply_lower_tri_self_transpose(L_id);
corr_matrix[2 * C] Cor_trial = multiply_lower_tri_self_transpose(L_trial);
vector[C] sd_ID_Cell = sd_id[1:C];
vector[C] sd_TrialNumber_Cell = sd_trial[1:C];
matrix[C, C] cor_ID_Cell = block(Cor_id, 1, 1, C, C);
matrix[C, C] cor_TrialNumber_Cell = block(Cor_trial, 1, 1, C, C);
vector[C] err_var_Cell = exp(2 * b_sigma_cell);
matrix[C, C] id_varcov = quad_form_diag(cor_ID_Cell, sd_ID_Cell);
matrix[C, C] trl_varcov = quad_form_diag(cor_TrialNumber_Cell, sd_TrialNumber_Cell);
}
