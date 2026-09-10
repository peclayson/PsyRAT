data {
  int<lower=0> NOBS; //number of obs
  int<lower=0> NSUB; //number of subj
  int<lower=1> NTRL; //number of trials
  array[NOBS] int<lower=0, upper=NSUB> id; //subject id
  array[NOBS] int<lower=1, upper=NTRL> trl; //trial id
  array[NOBS] real meas;
}
parameters {
  real mu;
  real<lower=0> sig_u;
  real<lower=0> sig_trl;
  real<lower=0> sig_e;
  vector[NSUB] u_raw;
  vector[NTRL] trl_raw;
}
transformed parameters {
  vector[NSUB] u;
  vector[NTRL] trl_terms;
  u = mu + sig_u*u_raw;
  trl_terms = sig_trl*trl_raw;
}
model {
  u_raw ~ normal(0,1);
  trl_raw ~ normal(0,1);
    meas ~ normal(u[id] + trl_terms[trl], sig_e);
  
  mu ~ normal(0,123.5);
  sig_u ~ cauchy(0,33);
  sig_e ~ cauchy(0,27);
  sig_trl ~ cauchy(0,10);
}
