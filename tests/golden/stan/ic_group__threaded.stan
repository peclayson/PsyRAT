functions {
  real partial_log_lik(array[] int seq_slice, int start, int end,
    array[] real meas, vector mu, real sigma) {
    real lp = 0;
    for (n in start:end) {
      lp += normal_lpdf(meas[n] | mu[n], sigma);
    }
    return lp;
  }
}
data {
  int<lower=0> NG1; //number of obs in g1
  int<lower=0> NG2; //number of obs in g2
  int<lower=1> grainsize_G1; //grainsize in g1
  int<lower=1> grainsize_G2; //grainsize in g2
  array[NG1] int<lower=1, upper=NG1> seq_G1;
  array[NG2] int<lower=1, upper=NG2> seq_G2;
  int<lower=0> JG1; //number of subj in g1
  int<lower=0> JG2; //number of subj in g2
  array[NG1] int<lower=0, upper=JG1> id_G1; //subject id g1;
  array[NG2] int<lower=0, upper=JG2> id_G2; //subject id g2;
  int<lower=1> NTRL_G1; //number of trials in g1
  int<lower=1> NTRL_G2; //number of trials in g2
  array[NG1] int<lower=1, upper=NTRL_G1> trl_G1; //trial id g1
  array[NG2] int<lower=1, upper=NTRL_G2> trl_G2; //trial id g2
  array[NG1] real meas_G1;
  array[NG2] real meas_G2;
}
parameters {
  real mu_G1;
  real mu_G2;
  real<lower=0> sig_u_G1;
  real<lower=0> sig_u_G2;
  real<lower=0> sig_trl_G1;
  real<lower=0> sig_trl_G2;
  real<lower=0> sig_e_G1;
  real<lower=0> sig_e_G2;
  vector[JG1] u_raw_G1;
  vector[JG2] u_raw_G2;
  vector[NTRL_G1] trl_raw_G1;
  vector[NTRL_G2] trl_raw_G2;
}
transformed parameters {
  vector[JG1] u_G1;
  vector[JG2] u_G2;
  vector[NTRL_G1] trl_terms_G1;
  vector[NTRL_G2] trl_terms_G2;
  u_G1 = mu_G1 + sig_u_G1*u_raw_G1;
  u_G2 = mu_G2 + sig_u_G2*u_raw_G2;
  trl_terms_G1 = sig_trl_G1*trl_raw_G1;
  trl_terms_G2 = sig_trl_G2*trl_raw_G2;
}
model {
  u_raw_G1 ~ normal(0,1);
  trl_raw_G1 ~ normal(0,1);
  target += reduce_sum(partial_log_lik, seq_G1, grainsize_G1, meas_G1, u_G1[id_G1] + trl_terms_G1[trl_G1], sig_e_G1);
  u_raw_G2 ~ normal(0,1);
  trl_raw_G2 ~ normal(0,1);
  target += reduce_sum(partial_log_lik, seq_G2, grainsize_G2, meas_G2, u_G2[id_G2] + trl_terms_G2[trl_G2], sig_e_G2);
  
  mu_G1 ~ normal(0,100);
  mu_G2 ~ normal(0,100);
  sig_u_G1 ~ cauchy(0,40);
  sig_u_G2 ~ cauchy(0,40);
  sig_e_G1 ~ cauchy(0,40);
  sig_e_G2 ~ cauchy(0,40);
  sig_trl_G1 ~ cauchy(0,10);
  sig_trl_G2 ~ cauchy(0,10);
}
