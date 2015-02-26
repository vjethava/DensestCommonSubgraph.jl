export erdos_renyi, planted_clique, get_dcs_test_instance

function erdos_renyi(n::Integer, p::Real; is_directed::Bool=false)
  A = sprandbool(n, n, p)
  if is_directed==false
    A = triu(A, 1)
    A = A + A'
  else
    A = triu(A, 1) + tril(A, -1)
  end
  return A
end

function planted_clique(n::Integer, p::Real, k::Integer; is_directed=false)
  A = erdos_renyi(n, p; is_directed = is_directed)
  for i=1:k
    for j=1:k
      if i != j
        A[i, j] = 1
      end
    end
  end
  return A
end

function get_dcs_test_instance(n::Integer = 15, p::Real = 0.1)
  if VERSION.minor == 3
    k = max(5, iceil(sqrt(n)))
  else
    k = max(5, ceil(Integer, sqrt(n)))
  end
  
  G = SparseMatrixCSC{Bool,Int64}[]
  push!(G,erdos_renyi(n, p))
  push!(G,planted_clique(n, p, k))

  return G
end

