data {
  int<lower=0> NPAIR;
  int<lower=1> NPAIR_SAFE;
  int<lower=0> NUNP;
  int<lower=1> NUNP_SAFE;
  int<lower=1> NSUB;
  int<lower=1> NTRL;
  vector[NPAIR_SAFE] meas1;
  vector[NPAIR_SAFE] meas2;
  array[NPAIR_SAFE] int<lower=1,upper=NSUB> ID_PAIR;
  array[NPAIR_SAFE] int<lower=1,upper=NTRL> TRL_PAIR;
  vector[NUNP_SAFE] meas_unp;
  array[NUNP_SAFE] int<lower=1,upper=2> EVENT_UNP;
  array[NUNP_SAFE] int<lower=1,upper=NSUB> ID_UNP;
  array[NUNP_SAFE] int<lower=1,upper=NTRL> TRL_UNP;
  int<lower=0,upper=1> estimate_rescor;
}
parameters {
  vector[2] b;
  vector[2] b_sigma;
  vector<lower=0>[2] sd_id;
  matrix[2, NSUB] z_id;
  cholesky_factor_corr[2] L_id;
  vector<lower=0>[2] sd_trl;
  matrix[2, NTRL] z_trl;
  cholesky_factor_corr[2] L_trl;
  vector<lower=0>[2] sd_id_sigma;
  matrix[2, NSUB] z_id_sigma;
  cholesky_factor_corr[2] L_id_sigma;
  real rescor_mu;
  real<lower=0> sd_rescor;
  vector[NSUB] z_rescor;
}
transformed parameters {
  matrix[NSUB, 2] r_id;
  matrix[NTRL, 2] r_trl;
  matrix[NSUB, 2] r_sigma;
  vector[NSUB] rho;
  r_id = (diag_pre_multiply(sd_id, L_id) * z_id)';
  r_trl = (diag_pre_multiply(sd_trl, L_trl) * z_trl)';
  r_sigma = (diag_pre_multiply(sd_id_sigma, L_id_sigma) * z_id_sigma)';
  for (s in 1:NSUB) {
    if (estimate_rescor == 1) {
      rho[s] = tanh(rescor_mu + sd_rescor * z_rescor[s]);
    } else {
      rho[s] = 0;
    }
  }
}
model {
  target += normal_lpdf(b | 0, 7);
  target += student_t_lpdf(b_sigma | 3, 0, 8);
  target += student_t_lpdf(sd_id | 3, 0, 9)
    - 2 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sd_trl | 3, 0, 11)
    - 2 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sd_id_sigma | 3, 0, 3)
    - 2 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(to_vector(z_id));
  target += std_normal_lpdf(to_vector(z_trl));
  target += std_normal_lpdf(to_vector(z_id_sigma));
  target += lkj_corr_cholesky_lpdf(L_id | 2);
  target += lkj_corr_cholesky_lpdf(L_trl | 2);
  target += lkj_corr_cholesky_lpdf(L_id_sigma | 2);
  target += normal_lpdf(rescor_mu | 0, 2);
  target += student_t_lpdf(sd_rescor | 3, 0, 1.5)
    - student_t_lccdf(0 | 3, 0, 1);
  target += std_normal_lpdf(z_rescor);
  if (NPAIR > 0) {
    for (n in 1:NPAIR) {
      int sid = ID_PAIR[n];
      int trl = TRL_PAIR[n];
      vector[2] Y;
      vector[2] Mu;
      vector[2] sigma;
      matrix[2,2] Lcorr;
      matrix[2,2] LSigma;
      Y = transpose([meas1[n], meas2[n]]);
      Mu[1] = b[1] + r_id[sid,1] + r_trl[trl,1];
      Mu[2] = b[2] + r_id[sid,2] + r_trl[trl,2];
      sigma[1] = exp(b_sigma[1] + r_sigma[sid,1]);
      sigma[2] = exp(b_sigma[2] + r_sigma[sid,2]);
      Lcorr[1,1] = 1;
      Lcorr[1,2] = 0;
      Lcorr[2,1] = rho[sid];
      Lcorr[2,2] = sqrt(1 - square(rho[sid]));
      LSigma = diag_pre_multiply(sigma, Lcorr);
      target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
    }
  }
  if (NUNP > 0) {
    for (n in 1:NUNP) {
      int e = EVENT_UNP[n];
      int sid = ID_UNP[n];
      int trl = TRL_UNP[n];
      real mu = b[e] + r_id[sid,e] + r_trl[trl,e];
      real sigma = exp(b_sigma[e] + r_sigma[sid,e]);
      target += normal_lpdf(meas_unp[n] | mu, sigma);
    }
  }
}
generated quantities {
  matrix[2,2] id_varcov;
  matrix[2,2] trl_varcov;
  matrix[2,2] err_varcov;
  matrix[NSUB,2] er_var_ss;
  vector[NSUB] wp_cov_ss;
  vector[1] rescor;
  corr_matrix[2] cor_id = multiply_lower_tri_self_transpose(L_id);
  corr_matrix[2] cor_trl = multiply_lower_tri_self_transpose(L_trl);
  id_varcov = quad_form_diag(cor_id, sd_id);
  trl_varcov = quad_form_diag(cor_trl, sd_trl);
  rescor[1] = mean(rho);
  err_varcov = rep_matrix(0, 2, 2);
  for (s in 1:NSUB) {
    er_var_ss[s,1] = b_sigma[1] + r_sigma[s,1];
    er_var_ss[s,2] = b_sigma[2] + r_sigma[s,2];
    wp_cov_ss[s] = rho[s] * exp(er_var_ss[s,1]) * exp(er_var_ss[s,2]);
    err_varcov[1,1] += exp(2 * er_var_ss[s,1]) / NSUB;
    err_varcov[2,2] += exp(2 * er_var_ss[s,2]) / NSUB;
    err_varcov[1,2] += wp_cov_ss[s] / NSUB;
    err_varcov[2,1] += wp_cov_ss[s] / NSUB;
  }
}
