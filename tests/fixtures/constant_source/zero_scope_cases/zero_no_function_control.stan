data {
  int<lower=1> NSUB;
  int<lower=0,upper=1> estimate_rescor;
}
parameters {
  vector[NSUB] rho;
  vector[NSUB] wp_cov_ss;
  vector[1] rescor;
  matrix[2,2] err_varcov;
}
model {
  rho ~ normal(0,1);
  wp_cov_ss ~ normal(0,1);
  rescor ~ normal(0,1);
  to_vector(err_varcov) ~ normal(0,1);
}
