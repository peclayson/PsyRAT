parameters { real z; }
model { z ~ normal(0,1); }
generated quantities {
  corr_matrix[2] R;
  matrix[2,2] alias;
  alias = block(R,1,1,2,2);
  R = diag_matrix(rep_vector(1,2));
}
