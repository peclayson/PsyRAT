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
  array[NE1] real meas_E1;
  array[NE2] real meas_E2;
}
parameters {
  real mu_E1;
  real mu_E2;
  real<lower=0> sig_u_E1;
  real<lower=0> sig_u_E2;
  real<lower=0> sig_trl_E1;
  real<lower=0> sig_trl_E2;
  real<lower=0> sig_e_E1;
  real<lower=0> sig_e_E2;
  vector[JE1] u_raw_E1;
  vector[JE2] u_raw_E2;
  vector[NTRL_E1] trl_raw_E1;
  vector[NTRL_E2] trl_raw_E2;
}
transformed parameters {
  vector[JE1] u_E1;
  vector[JE2] u_E2;
  vector[NTRL_E1] trl_terms_E1;
  vector[NTRL_E2] trl_terms_E2;
  u_E1 = mu_E1 + sig_u_E1*u_raw_E1;
  u_E2 = mu_E2 + sig_u_E2*u_raw_E2;
  trl_terms_E1 = sig_trl_E1*trl_raw_E1;
  trl_terms_E2 = sig_trl_E2*trl_raw_E2;
}
model {
  u_raw_E1 ~ normal(0,1);
  trl_raw_E1 ~ normal(0,1);
    meas_E1 ~ normal(u_E1[id_E1] + trl_terms_E1[trl_E1], sig_e_E1);
  u_raw_E2 ~ normal(0,1);
  trl_raw_E2 ~ normal(0,1);
    meas_E2 ~ normal(u_E2[id_E2] + trl_terms_E2[trl_E2], sig_e_E2);
  
  mu_E1 ~ normal(0,100);
  mu_E2 ~ normal(0,100);
  sig_u_E1 ~ cauchy(0,40);
  sig_u_E2 ~ cauchy(0,40);
  sig_e_E1 ~ cauchy(0,40);
  sig_e_E2 ~ cauchy(0,40);
  sig_trl_E1 ~ cauchy(0,10);
  sig_trl_E2 ~ cauchy(0,10);
}
