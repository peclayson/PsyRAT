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
  int<lower=0> NOBS; //number of obs
  int<lower=1> grainsize; //grainsize for reduce_sum
  array[NOBS] int<lower=1, upper=NOBS> seq; //index sequence
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
  vector[NOBS] y_hat = u[id] + trl_terms[trl];
  u_raw ~ normal(0,1);
  trl_raw ~ normal(0,1);
  target += reduce_sum(partial_log_lik, seq, grainsize, meas, y_hat, sig_e);
  
  mu ~ normal(0,100);
  sig_u ~ cauchy(0,40);
  sig_e ~ cauchy(0,40);
  sig_trl ~ cauchy(0,10);
}
