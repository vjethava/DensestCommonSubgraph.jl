using DensestCommonSubgraph
using FactCheck

FactCheck.onlystats(true)
FactCheck.setstyle(:default)

# write your own tests here
facts("Testing DensestCommonSubgraph") do 
  context("BinaryToGray conversion") do
    include("test_bitstr.jl")
    n = 5
    test_binary_to_gray(n)
  end

  context("DensityFunctions") do
    n = 1000

    p = 0.4
    k = 200
    G = erdos_renyi(n, p)
    for i=1:10
      ns = rand((n/10):n)
      S = randperm(n)[1:ns]
      d = subgraph_density(G, S)
      r = subgraph_relative_density(G, S)
      @fact d => roughly((ns-1)/2 * p; atol=1) "subgraph_density: $d expected: $((ns-1)/2 * p)"
      @fact r => roughly(p; atol=0.1) "relative_density: $r expected: $p"
    end
    
    
    H = planted_clique(n, p,  k)
    A = SparseMatrixCSC{Bool,Int64}[]
    push!(A, H)
    push!(A, G)
    for i=1:10
      ns = rand(1:n)
      S = randperm(n)[1:ns]     
      dg = subgraph_density(G, S)
      dh = subgraph_density(H, S)
      da = common_subgraph_density(A, S)
      @fact da => min(dh, dg)
    end

  end
          
  context("Faster DCS LP") do
    n= 30
    p= 1/n
    for i=1:10
      G  = get_dcs_test_instance(n, p);
      @time p2 = fast_dcs_lp(G);
      @time p3 = faster_dcs_lp(G); 
      @time p1 = dcs_lp(G);
      println("LP: $(p1.optval), FAST: $(p2.optval) FASTER: $(p3.optval)") 
      @fact p1.optval => roughly(p2.optval, atol=1e-2) "Error in fast LP!"
    end
  end
         
end

  
  
  
  
  
  
  
  
