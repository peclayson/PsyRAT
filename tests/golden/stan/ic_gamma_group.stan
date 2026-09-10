functions {
  real scaled_chi_square_lpdf(real y, real mu, real nu) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(nu) - log(mu) + chi_square_lpdf(y * nu / mu | nu);
  }
}
data {
  int<lower=0> NG1; //number of obs in g1
  int<lower=0> NG2; //number of obs in g2
  int<lower=0> JG1; //number of subj in g1
  int<lower=0> JG2; //number of subj in g2
  array[NG1] int<lower=0, upper=JG1> id_G1; //subject id g1;
  array[NG2] int<lower=0, upper=JG2> id_G2; //subject id g2;
  int<lower=1> NTRL_G1; //number of trials in g1
  int<lower=1> NTRL_G2; //number of trials in g2
  array[NG1] int<lower=1, upper=NTRL_G1> trl_G1; //trial id g1
  array[NG2] int<lower=1, upper=NTRL_G2> trl_G2; //trial id g2
  vector[NG1] meas_G1;
  vector[NG2] meas_G2;
}
parameters {
  real Intercept_G1;  // log-mean grand intercept
  real Intercept_nu_G1;  // log-nu grand intercept
  vector<lower=0>[2] gro_sds_G1;  // person SDs: [mean, log-nu]
  matrix[2, JG1] gro_effs_stndzd_G1;
  cholesky_factor_corr[2] chol_corrmat_G1;
  real<lower=0> sig_trl_G1;  // trial main-effect log-SD
  vector[NTRL_G1] trl_raw_G1;
  real Intercept_G2;  // log-mean grand intercept
  real Intercept_nu_G2;  // log-nu grand intercept
  vector<lower=0>[2] gro_sds_G2;  // person SDs: [mean, log-nu]
  matrix[2, JG2] gro_effs_stndzd_G2;
  cholesky_factor_corr[2] chol_corrmat_G2;
  real<lower=0> sig_trl_G2;  // trial main-effect log-SD
  vector[NTRL_G2] trl_raw_G2;
}
transformed parameters {
  matrix[JG1, 2] gro_effs_actual_G1;
  vector[JG1] ind_bs_G1;
  vector[JG1] ind_nu_G1;
  vector[NTRL_G1] trl_b_G1;
  matrix[JG2, 2] gro_effs_actual_G2;
  vector[JG2] ind_bs_G2;
  vector[JG2] ind_nu_G2;
  vector[NTRL_G2] trl_b_G2;
  gro_effs_actual_G1 = (diag_pre_multiply(gro_sds_G1, chol_corrmat_G1) * gro_effs_stndzd_G1)';
  ind_bs_G1 = gro_effs_actual_G1[, 1];
  ind_nu_G1 = gro_effs_actual_G1[, 2];
  trl_b_G1 = sig_trl_G1 * trl_raw_G1;
  gro_effs_actual_G2 = (diag_pre_multiply(gro_sds_G2, chol_corrmat_G2) * gro_effs_stndzd_G2)';
  ind_bs_G2 = gro_effs_actual_G2[, 1];
  ind_nu_G2 = gro_effs_actual_G2[, 2];
  trl_b_G2 = sig_trl_G2 * trl_raw_G2;
}
model {
  vector[NG1] mu_G1 = Intercept_G1 + rep_vector(0, NG1);
  vector[NG1] nu_G1 = Intercept_nu_G1 + rep_vector(0, NG1);
  mu_G1 += ind_bs_G1[id_G1];
  mu_G1 += trl_b_G1[trl_G1];
  mu_G1 = exp(mu_G1);
  nu_G1 += ind_nu_G1[id_G1];
  nu_G1 = exp(nu_G1);
  vector[NG2] mu_G2 = Intercept_G2 + rep_vector(0, NG2);
  vector[NG2] nu_G2 = Intercept_nu_G2 + rep_vector(0, NG2);
  mu_G2 += ind_bs_G2[id_G2];
  mu_G2 += trl_b_G2[trl_G2];
  mu_G2 = exp(mu_G2);
  nu_G2 += ind_nu_G2[id_G2];
  nu_G2 = exp(nu_G2);
  gro_effs_stndzd_G1[1] ~ std_normal();
  gro_effs_stndzd_G1[2] ~ std_normal();
  trl_raw_G1 ~ std_normal();
  Intercept_G1 ~ normal(1.70475, 0.5);
  Intercept_nu_G1 ~ normal(2.89037, 0.75);
  gro_sds_G1[1] ~ student_t(3, 0, 0.35);
  gro_sds_G1[2] ~ student_t(3, 0, 1);
  sig_trl_G1 ~ student_t(3, 0, 0.35);
  chol_corrmat_G1 ~ lkj_corr_cholesky(1);
  for (n in 1:NG1) {
    target += scaled_chi_square_lpdf(meas_G1[n] | mu_G1[n], nu_G1[n]);
  }
  gro_effs_stndzd_G2[1] ~ std_normal();
  gro_effs_stndzd_G2[2] ~ std_normal();
  trl_raw_G2 ~ std_normal();
  Intercept_G2 ~ normal(1.70475, 0.5);
  Intercept_nu_G2 ~ normal(2.89037, 0.75);
  gro_sds_G2[1] ~ student_t(3, 0, 0.35);
  gro_sds_G2[2] ~ student_t(3, 0, 1);
  sig_trl_G2 ~ student_t(3, 0, 0.35);
  chol_corrmat_G2 ~ lkj_corr_cholesky(1);
  for (n in 1:NG2) {
    target += scaled_chi_square_lpdf(meas_G2[n] | mu_G2[n], nu_G2[n]);
  }
}
