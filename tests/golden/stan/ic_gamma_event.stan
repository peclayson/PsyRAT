functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
}
data {
  int<lower=0> NE1; //number of obs in E1
  int<lower=0> NE2; //number of obs in E2
  int<lower=0> JE1; //number of subj in E1
  int<lower=0> JE2; //number of subj in E2
  array[NE1] int<lower=0, upper=JE1> id_E1; //subject id E1;
  array[NE2] int<lower=0, upper=JE2> id_E2; //subject id E2;
  int<lower=1> NTRL_E1; //number of trials in E1
  int<lower=1> NTRL_E2; //number of trials in E2
  array[NE1] int<lower=1, upper=NTRL_E1> trl_E1; //trial id E1
  array[NE2] int<lower=1, upper=NTRL_E2> trl_E2; //trial id E2
  vector[NE1] meas_E1;
  vector[NE2] meas_E2;
}
parameters {
  real Intercept_E1;  // log-mean grand intercept
  real Intercept_nu_E1;  // log-nu grand intercept
  vector<lower=0>[2] gro_sds_E1;  // person SDs: [mean, log-nu]
  matrix[2, JE1] gro_effs_stndzd_E1;
  cholesky_factor_corr[2] chol_corrmat_E1;
  real<lower=0> sig_trl_E1;  // trial main-effect log-SD
  vector[NTRL_E1] trl_raw_E1;
  real Intercept_E2;  // log-mean grand intercept
  real Intercept_nu_E2;  // log-nu grand intercept
  vector<lower=0>[2] gro_sds_E2;  // person SDs: [mean, log-nu]
  matrix[2, JE2] gro_effs_stndzd_E2;
  cholesky_factor_corr[2] chol_corrmat_E2;
  real<lower=0> sig_trl_E2;  // trial main-effect log-SD
  vector[NTRL_E2] trl_raw_E2;
}
transformed parameters {
  matrix[JE1, 2] gro_effs_actual_E1;
  vector[JE1] ind_bs_E1;
  vector[JE1] ind_nu_E1;
  vector[NTRL_E1] trl_b_E1;
  matrix[JE2, 2] gro_effs_actual_E2;
  vector[JE2] ind_bs_E2;
  vector[JE2] ind_nu_E2;
  vector[NTRL_E2] trl_b_E2;
  gro_effs_actual_E1 = (diag_pre_multiply(gro_sds_E1, chol_corrmat_E1) * gro_effs_stndzd_E1)';
  ind_bs_E1 = gro_effs_actual_E1[, 1];
  ind_nu_E1 = gro_effs_actual_E1[, 2];
  trl_b_E1 = sig_trl_E1 * trl_raw_E1;
  gro_effs_actual_E2 = (diag_pre_multiply(gro_sds_E2, chol_corrmat_E2) * gro_effs_stndzd_E2)';
  ind_bs_E2 = gro_effs_actual_E2[, 1];
  ind_nu_E2 = gro_effs_actual_E2[, 2];
  trl_b_E2 = sig_trl_E2 * trl_raw_E2;
}
model {
  vector[NE1] mu_E1 = Intercept_E1 + rep_vector(0, NE1);
  vector[NE1] nu_E1 = Intercept_nu_E1 + rep_vector(0, NE1);
  mu_E1 += ind_bs_E1[id_E1];
  mu_E1 += trl_b_E1[trl_E1];
  mu_E1 = exp(mu_E1);
  nu_E1 += ind_nu_E1[id_E1];
  nu_E1 = exp(nu_E1);
  vector[NE2] mu_E2 = Intercept_E2 + rep_vector(0, NE2);
  vector[NE2] nu_E2 = Intercept_nu_E2 + rep_vector(0, NE2);
  mu_E2 += ind_bs_E2[id_E2];
  mu_E2 += trl_b_E2[trl_E2];
  mu_E2 = exp(mu_E2);
  nu_E2 += ind_nu_E2[id_E2];
  nu_E2 = exp(nu_E2);
  gro_effs_stndzd_E1[1] ~ std_normal();
  gro_effs_stndzd_E1[2] ~ std_normal();
  trl_raw_E1 ~ std_normal();
  Intercept_E1 ~ normal(1.70475, 0.5);
  Intercept_nu_E1 ~ normal(2.89037, 0.75);
  gro_sds_E1[1] ~ student_t(3, 0, 0.35);
  gro_sds_E1[2] ~ student_t(3, 0, 1);
  sig_trl_E1 ~ student_t(3, 0, 0.35);
  chol_corrmat_E1 ~ lkj_corr_cholesky(1);
  for (n in 1:NE1) {
    target += scaled_chi_square_lpdf(meas_E1[n] | mu_E1[n], nu_E1[n]);
  }
  gro_effs_stndzd_E2[1] ~ std_normal();
  gro_effs_stndzd_E2[2] ~ std_normal();
  trl_raw_E2 ~ std_normal();
  Intercept_E2 ~ normal(1.70475, 0.5);
  Intercept_nu_E2 ~ normal(2.89037, 0.75);
  gro_sds_E2[1] ~ student_t(3, 0, 0.35);
  gro_sds_E2[2] ~ student_t(3, 0, 1);
  sig_trl_E2 ~ student_t(3, 0, 0.35);
  chol_corrmat_E2 ~ lkj_corr_cholesky(1);
  for (n in 1:NE2) {
    target += scaled_chi_square_lpdf(meas_E2[n] | mu_E2[n], nu_E2[n]);
  }
}
