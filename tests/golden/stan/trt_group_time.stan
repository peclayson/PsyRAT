data { 
  int<lower=1> NOBS; //total number of observations 
  int<lower=1> NSUB; //total number of subjects
  int<lower=1> NOCC; //total number of occasions
  int<lower=1> NTRL; //total number of trials
  int<lower=1> NTID; // total unique for trial*id
  int<lower=1> NOID; // total unique for occ*id
  int<lower=1> NTO; // total unique for trial*occ
  array[NOBS] int<lower=1, upper=NSUB> id; //id variable
  array[NOBS] int<lower=1, upper=NOCC> occ; //occ variable
  array[NOBS] int<lower=1, upper=NTRL> trl; //trl variable
  array[NOBS] int<lower=1, upper=NTID> trlxid; //trlxid variable
  array[NOBS] int<lower=1, upper=NOID> occxid; //occxid variable
  array[NOBS] int<lower=1, upper=NTO> trlxocc; //trlxocc variable
  vector[NOBS] meas; //measurements
}
parameters {
  real pop_int; //population intercept
  real<lower=0> sig_id; //id-level std dev
  real<lower=0> sig_occ; //occ-level std dev
  real<lower=0> sig_trl; //trl-level std dev
  real<lower=0> sig_err; //residual std dev
  real<lower=0> sig_trlxid; //trlxid std dev
  real<lower=0> sig_occxid; //occxid std dev
  real<lower=0> sig_trlxocc; //trlxocc std dev
  vector[NSUB] id_raw; //id means
  vector[NOCC] occ_raw; //occ means
  vector[NTRL] trl_raw; //trl means
  vector[NTID] trlxid_raw; //trlxid means
  vector[NOID] occxid_raw; //occxid means
  vector[NTO] trlxocc_raw; //trlxocc means
}
transformed parameters {
  vector[NSUB] id_terms; //id-level terms
  vector[NOCC] occ_terms; //occ-level terms
  vector[NTRL] trl_terms; //trl-level terms
  vector[NTID] trlxid_terms; //trlxid terms
  vector[NOID] occxid_terms; //trlxocc terms
  vector[NTO] trlxocc_terms; //trlxocc terms
  id_terms = sig_id * id_raw;
  occ_terms = sig_occ * occ_raw;
  trl_terms = sig_trl * trl_raw;
  trlxid_terms = sig_trlxid * trlxid_raw;
  occxid_terms = sig_occxid * occxid_raw;
  trlxocc_terms = sig_trlxocc * trlxocc_raw;
}
model {
  vector[NOBS] y_hat;
    y_hat = pop_int + id_terms[id] + 
    occ_terms[occ] + trl_terms[trl] +
    trlxid_terms[trlxid] + occxid_terms[occxid] + 
    trlxocc_terms[trlxocc];
  meas ~ normal(y_hat,sig_err);
  id_raw ~ normal(0,1);
  occ_raw ~ normal(0,1);
  trl_raw ~ normal(0,1);
  trlxid_raw ~ normal(0,1);
  occxid_raw ~ normal(0,1);
  trlxocc_raw ~ normal(0,1);
  sig_id ~ cauchy(0,10);
  sig_occ ~ cauchy(0,5);
  sig_trl ~ cauchy(0,10);
  sig_err ~ cauchy(0,20);
  sig_trlxid ~ cauchy(0,5);
  sig_occxid ~ cauchy(0,5);
  sig_trlxocc ~ cauchy(0,1);
}
