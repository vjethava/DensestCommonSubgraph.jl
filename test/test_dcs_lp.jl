
using DensestCommonSubgraph
using FactCheck
using Gadfly
using MAT
using DataFrames
using ProgressMeter
FactCheck.onlystats(true)
FactCheck.setstyle(:default)
using HDF5, JLD

function test_instances2(m::Integer, n::Integer, p::Real, k::Integer)
  G = Graph[]
  for i=1:m
    cG = planted_clique(n, p, k)
    push!(G , cG)
  end
  return G
end

function test_dcs_lp(nt::Integer=10 )
  m = 1
  Z = zeros(nt, m)
  D = zeros(nt, 4)  # Density of extracted subgraph
  V = zeros(nt, 4)  # LP optimum
  S = cell(nt, 4)   # subgraphs extracted
  I = cell(nt, 1)   # graph instances
  X = DataFrame(GID=Int64[] , RelDensity=Float64[], RelOptVal=Float64[], Jaccard=Float64[], Method=String[])
  p_bar = Progress(nt, 1, "Computing initial pass...", 50)
  for i=1:nt
    n = 20
    p = 0.5
    k = rand(1:ifloor(n/2))
    G = test_instances2(m, n, p, k);
    I[i] = G
    (s_true, d_val ) = dcs_parallel(G, 2);
    d_true = common_subgraph_density(G, s_true)
    for jj=1:m
      Z[i, jj] = subgraph_density(G[jj], s_true)
      println("i: $i j: $jj $(Z[i, jj])")
    end

    # @fact d_true => roughly(d_val/2.0)
    # println("OPT: density: $(lpad(d_true, 20, ' ')) subgraph: $s_true")

    D[i, 1] = d_true
    V[i, 1] = d_val
    S[i, 1] = s_true
    (p_dcs_lp, X_dcs, Y_dcs, t_dcs) = dcs_lp(G);
    s_dcs = extract_primal_solution(Y_dcs)
    d_dcs = common_subgraph_density(G, s_dcs)
    # println("DCS: density: $(lpad(d_dcs, 20, ' ')) subgraph: $s_dcs")

    D[i, 2] = d_dcs
    V[i, 2] = t_dcs
    S[i, 2] = s_dcs
    if t_dcs > d_val
      println("i: $i integrality gap: dcs: $t_dcs opt: $d_val  ratio: $(t_dcs/d_val)")
    end
    X = vcat(X, DataFrame(GID=i, RelDensity = d_dcs/d_true , RelOptVal=t_dcs/d_val, Jaccard=length(intersect(s_dcs , s_true))/length(union(s_dcs, s_true)), Method="DCS"))

    # @fact d_dcs => less_than_or_equal(d_true)

    (p_charikar_lp, Xc, Yc, t_c) = charikar_lp(dcs_intersection_graph(G))
    s_c = extract_primal_solution(Yc)
    d_c = common_subgraph_density(G, s_c)

    D[i, 3] = d_c
    V[i, 3] = t_c
    S[i, 3] = s_c

    X = vcat(X, DataFrame(GID=i, RelDensity = d_c/d_true , RelOptVal=t_c/d_val, Jaccard=length(intersect(s_c , s_true))/length(union(s_c, s_true)), Method="MIN"))

    # println("CHARIKAR: density: $(lpad(d_c, 20, ' ')) subgraph: $s_c")

    @fact d_c => less_than_or_equal(d_true)

    (p_avg_lp, Xavg, Yavg, t_avg) = charikar_lp(dcs_avg_graph(G))
    s_avg = extract_primal_solution(Yavg)
    d_avg = common_subgraph_density(G, s_avg)
    # println("AVG: density: $(lpad(d_avg, 20, ' ')) subgraph: $s_avg")

    D[i, 4] = d_avg
    V[i, 4] = t_avg
    S[i, 4] = s_avg

    X = vcat(X, DataFrame(GID=i, RelDensity = d_avg/d_true , RelOptVal=t_avg/d_val, Jaccard=length(intersect(s_avg , s_true))/length(union(s_avg, s_true)), Method="AVG"))

    # @fact d_avg => less_than_or_equal(d_true)
    println("$i/$nt opt: $(rpad(d_true, 20, ' ')) DCS: $(rpad(d_dcs, 20, ' ')) CHARIKAR: $(rpad(d_c, 20, ' ')) AVG: $(p_avg_lp.optval)\n")
    # println("\n")
    next!(p_bar)
  end
  return (I, X, S, Z)
end

# function evaluate_primal(G::GraphVec, methods::String[]= LP_METHODS)
#   m = length(G)
#   @assert m >= 1
#   n = size(G, 1)
#   S = Array{Int64, 1}[]
#   X = DataFrame(GID=Int64[] , RelDensity=Float64[], RelOptVal=Float64[], Jaccard=Float64[], Method=String[])
#   (s_true, d_val) = dcs_parallel(G, 2);
#   push!(S, s_true)
#   d_true = common_subgraph_density(G, s_true)
#   for m=methods
#     (cS, d, t) = dcs_lp_extract(G; method = m)
#     push!(S, cS)
#     X = vcat(X, DataFrame(GID=gid, RelDensity = d/d_true , RelOptVal=t/d_val, Jaccard=length(intersect(S , s_true))/length(union(cS, s_true)), Method=m))
#   end
#   return
# end

nt = 100;
(I, X, S, Z) = test_dcs_lp(nt)
legend = "I: graph instances, S: subgraphs returned by OPT,DCS,LP_min,LP_avg, X: stats, Z: densities by S_opt"
using HDF5, JLD
@save "output_test_dcs_lp_2.jld" I X S Z legend
p1 = plot(X, xgroup="Method", x="RelOptVal", Geom.subplot_grid(Geom.histogram), Guide.xlabel(nothing), Guide.title("Ratio of LP optimum to density of S_{OPT}"), color="Method")
p2 = plot(X, xgroup="Method", x="RelDensity", Geom.subplot_grid(Geom.histogram), Guide.xlabel(nothing), Guide.title("Ratio of densities of subgraph obtained by LP rounding to S_{OPT}"), color="Method")
draw(PDF("p1a.pdf", 9inch, 6inch), p1)
draw(PDF("p2a.pdf", 9inch, 6inch), p2)
