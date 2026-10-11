functions {
  real unused_zero_graph(int NSUB, int estimate_rescor, real rescor_mu,
                         real sd_rescor, vector z_rescor) {
    vector[NSUB] rho;
    matrix[NSUB,2] er_var_ss = rep_matrix(0,NSUB,2);
    vector[NSUB] wp_cov_ss;
    vector[1] rescor;
    matrix[2,2] err_varcov;
    for (s in 1:NSUB) {
      if (estimate_rescor == 1) {
        rho[s] = tanh(rescor_mu + sd_rescor * z_rescor[s]);
      } else {
        rho[s] = 0;
      }
    }
    rescor[1] = mean(rho);
    err_varcov = rep_matrix(0,2,2);
    for (s in 1:NSUB) {
      wp_cov_ss[s] = rho[s] * exp(er_var_ss[s,1]) * exp(er_var_ss[s,2]);
      err_varcov[1,1] += exp(2 * er_var_ss[s,1]) / NSUB;
      err_varcov[2,2] += exp(2 * er_var_ss[s,2]) / NSUB;
      err_varcov[1,2] += wp_cov_ss[s] / NSUB;
      err_varcov[2,1] += wp_cov_ss[s] / NSUB;
    }
    return rescor[1] + sum(wp_cov_ss) + sum(to_vector(err_varcov));
  }
}
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
