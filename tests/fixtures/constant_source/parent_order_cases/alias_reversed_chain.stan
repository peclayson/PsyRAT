parameters { corr_matrix[2] R; }
model { R ~ lkj_corr(2); }
generated quantities {
  matrix[2,2] alias;
  matrix[2,2] second;
  second = block(alias,1,1,2,2);
  alias = block(R,1,1,2,2);
}
