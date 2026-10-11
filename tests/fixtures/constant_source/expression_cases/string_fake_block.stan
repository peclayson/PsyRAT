parameters { real z; }
model { z ~ normal(0,1); }
generated quantities {
  matrix[2,2] pretend = rep_matrix(z,2,2);
  print("} generated quantities { corr_matrix[2] pretend;");
}
