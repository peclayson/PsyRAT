parameters {
  corr_matrix[2] R;
  real z;
}
model {
  R ~ lkj_corr(2);
  z ~ normal(0,1);
}
generated quantities {
  matrix[2,2] alias;
  if (z > 0) alias = block(R,1,1,2,2);
}
