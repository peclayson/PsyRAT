parameters { corr_matrix[2] R; }
model { R ~ lkj_corr(2); }
generated quantities {
  matrix[2,2] alias;
  print("ordinary text /* not a comment */ // } generated quantities {");
  alias = block(R,1,1,2,2);
}
