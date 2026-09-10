data { 
  int<lower=1> NOBS; //total number of observations 
  int<lower=1> NSUB; //total number of subjects
  int<lower=1> NOCC; //total number of occasions
  int<lower=1> NTRL; //total number of splits
  int<lower=1> NTID; // total unique for split*id
  int<lower=1> NOID; // total unique for occ*id
  int<lower=1> NTO; // total unique for split*occ
  array[NOBS] int<lower=1, upper=NSUB> id; //id variable
  array[NOBS] int<lower=1, upper=NOCC> occ; //occ variable
  array[NOBS] int<lower=1, upper=NTRL> trl; //split variable
  array[NOBS] int<lower=1, upper=NTID> trlxid; //split*id variable
  array[NOBS] int<lower=1, upper=NOID> occxid; //occ*id variable
  array[NOBS] int<lower=1, upper=NTO> trlxocc; //split*occ variable
  vector[NOBS] meas; //measurements
  vector<lower=1>[NOBS] weight; //items per split (n_i), integer >= 1
}
parameters {
  real pop_int; //population intercept
  real<lower=0> sig_id; //id-level std dev
  real<lower=0> sig_occ; //occ-level std dev
  real<lower=0> sig_trl; //split-level std dev
  real<lower=0> sig_err; //per-item residual std dev
  real<lower=0> sig_trlxid; //split*id std dev (p:s)
  real<lower=0> sig_occxid; //occ*id std dev (p:o)
  real<lower=0> sig_trlxocc; //split*occ std dev (o:s)
  real<lower=0> sig_posxid; //occ*split*id std dev (p:o:s)
  vector[NSUB] id_raw; //id means
  vector[NOCC] occ_raw; //occ means
  vector[NTRL] trl_raw; //split means
  vector[NTID] trlxid_raw; //split*id means
  vector[NOID] occxid_raw; //occ*id means
  vector[NTO] trlxocc_raw; //split*occ means
  vector[NOBS] posxid_raw; //occ*split*id means (p:o:s, per-row)
}
transformed parameters {
  vector[NSUB] id_terms; //id-level terms
  vector[NOCC] occ_terms; //occ-level terms
  vector[NTRL] trl_terms; //split-level terms
  vector[NTID] trlxid_terms; //split*id terms
  vector[NOID] occxid_terms; //occ*id terms
  vector[NTO] trlxocc_terms; //split*occ terms
  vector[NOBS] posxid_terms; //occ*split*id terms (p:o:s)
  vector[NOBS] sig_obs; //n_i-weighted per-obs residual SD
  id_terms = sig_id * id_raw;
  occ_terms = sig_occ * occ_raw;
  trl_terms = sig_trl * trl_raw;
  trlxid_terms = sig_trlxid * trlxid_raw;
  occxid_terms = sig_occxid * occxid_raw;
  trlxocc_terms = sig_trlxocc * trlxocc_raw;
  posxid_terms = sig_posxid * posxid_raw;
  sig_obs = sig_err ./ sqrt(weight);
}
model {
  vector[NOBS] y_hat;
    y_hat = pop_int + id_terms[id] + 
    occ_terms[occ] + trl_terms[trl] +
    trlxid_terms[trlxid] + occxid_terms[occxid] + 
    trlxocc_terms[trlxocc] + posxid_terms;
  meas ~ normal(y_hat,sig_obs);
  id_raw ~ normal(0,1);
  occ_raw ~ normal(0,1);
  trl_raw ~ normal(0,1);
  trlxid_raw ~ normal(0,1);
  occxid_raw ~ normal(0,1);
  trlxocc_raw ~ normal(0,1);
  posxid_raw ~ normal(0,1);
  sig_id ~ cauchy(0,10);
  sig_occ ~ cauchy(0,5);
  sig_trl ~ cauchy(0,10);
  sig_err ~ cauchy(0,20);
  sig_trlxid ~ cauchy(0,5);
  sig_occxid ~ cauchy(0,5);
  sig_trlxocc ~ cauchy(0,1);
  sig_posxid ~ cauchy(0,5);
}
