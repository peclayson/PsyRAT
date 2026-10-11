parameters { corr_matrix[2] R; }
model { R ~ lkj_corr(2); }
generated quantities {
  matrix[2,2] alias;
  print("one literal backslash \");
  alias = block(R,1,1,2,2);
}
