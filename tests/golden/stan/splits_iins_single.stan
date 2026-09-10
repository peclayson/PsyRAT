data {
  int<lower=0> NOBS; //number of obs
  int<lower=0> NSUB; //number of subj
  int<lower=1> NTRL; //number of splits
  array[NOBS] int<lower=0, upper=NSUB> id; //subject id
  array[NOBS] int<lower=1, upper=NTRL> trl; //split id
  array[NOBS] real meas;
  vector<lower=1>[NOBS] weight; //items per split (n_i), integer >= 1
}
parameters {
  real mu;
  real<lower=0> sig_u;
  real<lower=0> sig_trl;
  real<lower=0> sig_splitxid;
  real<lower=0> sig_err;
  vector[NSUB] u_raw;
  vector[NTRL] trl_raw;
  vector[NOBS] splitxid_raw;
}
transformed parameters {
  vector[NSUB] u;
  vector[NTRL] trl_terms;
  vector[NOBS] splitxid_terms;
  vector[NOBS] sig_obs;
  u = mu + sig_u*u_raw;
  trl_terms = sig_trl*trl_raw;
  splitxid_terms = sig_splitxid*splitxid_raw;
  sig_obs = sig_err ./ sqrt(weight);
}
model {
  u_raw ~ normal(0,1);
  trl_raw ~ normal(0,1);
  splitxid_raw ~ normal(0,1);
    meas ~ normal(u[id] + trl_terms[trl] + splitxid_terms, sig_obs);
  
  mu ~ normal(0,100);
  sig_u ~ cauchy(0,40);
  sig_trl ~ cauchy(0,10);
  sig_err ~ cauchy(0,40);
  sig_splitxid ~ cauchy(0,5);
}
