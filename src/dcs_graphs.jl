export Graph, erdos_renyi, planted_clique, get_dcs_test_instance

typealias TGraph{T} SparseMatrixCSC{T,Int64}
typealias Graph TGraph{Float64}
typealias GraphVec Array{Graph, 1}


function erdos_renyi(n::Integer, p::Real; is_directed::Bool=false)
  A = spones(sprand(n, n, p))
  if is_directed==false
    A = triu(A, 1)
    A = A + A'
  else
    for i=1:n
      A[i, i]=0.0
    end
  end
  return A
end

function erdos_renyi_fast(n::Integer, p::Real; is_directed::Bool=false)
  A = spzeros(Bool, n, n)
  if is_directed
    for i=1:n
      @simd for j=(i+1):n
        @inbounds A[i, j] = rand() < p
        @inbounds A[j, i] = rand() < p
      end
    end
  else
    for i=1:n
      B = sprandbool(1, n - i, p)
      # bi = find(B)
      @simd for j= 1 : (n-i)
        @inbounds A[i, j+i] = B[j]
        @inbounds A[j+i, i] = B[j]
      end

    end
  end
  return A
end


function planted_clique(n::Integer, p::Real, k::Integer; is_directed=false)
  A = erdos_renyi(n, p; is_directed = is_directed)
  S = randperm(n)[1:k;]
  for i=1:k
    for j=1:k
      if i != j
        A[S[i], S[j]] = 1
      end
    end
  end
  return A
end

function get_dcs_test_instance(n::Integer = 15, p::Real = 0.1, k::Integer = 0)
  if k==0 && VERSION.minor == 3
    k = max(5, iceil(sqrt(n)))
  elseif k==0
    k = max(5, ceil(Integer, sqrt(n)))
  end

  G = Graph[]
  push!(G,erdos_renyi(n, p))
  push!(G,planted_clique(n, p, k))

  return G
end

