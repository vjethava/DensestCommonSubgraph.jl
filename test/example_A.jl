using DensestCommonSubgraph

for i=1:10
  d_true = 0.0; 
  G = get_dcs_test_instance(30);
  # (s_true, d_true) = dcs_parallel(G, 2);
  @time p_dcs_lp = dcs_lp(G);
  @time p_charikar_lp = dcs_charikar_lp(G);
  println("opt: $d_true DCS: $(p_dcs_lp.optval) CHARIKAR: $(p_charikar_lp.optval)")
end
