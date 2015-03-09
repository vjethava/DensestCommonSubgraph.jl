using DensestCommonSubgraph

for i=1:10
  d_true = 0.0; 
  G = get_dcs_test_instance(15);
  (s_true, d_true) = dcs_parallel(G, 2);
  p_dcs_lp, t_dcs, X_dcs, Y_dcs = dcs_lp(G);
  p_charikar_lp, t_c, X_c, Y_c = charikar_lp(dcs_intersection_graph(G));
  p_avg_lp, t_avg, X_avg, Y_avg = charikar_lp(dcs_avg_graph(G));
  println("opt: $d_true DCS: $(p_dcs_lp.optval) CHARIKAR: $(p_charikar_lp.optval) AVG: $(p_avg_lp.optval)")
  println("opt: $d_true DCS: $(p_dcs_lp.optval) CHARIKAR: $(p_charikar_lp.optval)")
end
