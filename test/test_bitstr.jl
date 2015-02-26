using DensestCommonSubgraph
using Base.Test
using FactCheck 

function test_binary_to_gray(n::Integer=5)
  prev_code = "" 
  for i=0:(2^n-1)
    binary_code = bin(i, n)
    gray_code = binary_to_gray(binary_code)
    orig_code = gray_to_binary(gray_code)
    @fact orig_code => binary_code "test_binary_to_gray(): incorrect conversion of $i"
    if i >= 1 
      xx = filter((j) -> gray_code[j] != prev_code[j], [1:n])
      @fact length(xx) => 1 "test_binary_to_gray(): incorrect gray coding $(i-1), $i" 
    end
    
    # println("i: ", i, "\tbinary: ", binary_code, "\tgray_code: ", gray_code, "\toriginal", orig_code)
    prev_code = gray_code
  end
end

