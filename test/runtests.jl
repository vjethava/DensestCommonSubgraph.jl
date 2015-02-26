using DensestCommonSubgraph
using FactCheck

# write your own tests here
facts("Testing DensestCommonSubgraph") do 
  context("BinaryToGray conversion") do
    include("test_bitstr.jl")
    n = 10
    test_binary_to_gray(n)
  end
end

