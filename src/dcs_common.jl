export subgraph_density, subgraph_relative_density, common_subgraph_density, dcs_intersection_graph, dcs_union_graph, dcs_avg_graph
export in_degree, out_degree, degree

function getBytes(x)
  total = 0;
  fieldNames = typeof(x).names;
  if fieldNames == ()
    return sizeof(x);
  else
    for fieldName in fieldNames
      try
        total += getBytes(getfield(x,fieldName));
        println("getBytes()$fieldName, $total")
      catch e
        continue
      end
    end
    return total;
  end
end


function in_degree{T}(G::SparseMatrixCSC{T, Int64}, S::Array{Int64, 1}=Int64[])
  n = size(G, 1)
  if length(S) == 0
    S = 1:n
  end
  ns = length(S)
  result = zeros(Float64, ns , 1)
  H = G[S, S]
  (rr, cc) = findn(H)
  ne = length(rr)
  for i=1:ne
    crow = rr[i]
    ccol = cc[i]
    result[ccol] = result[ccol] + H[crow, ccol]
  end
  return result
end

function out_degree{T}(G::SparseMatrixCSC{T, Int64}, S::Array{Int64, 1}=Int64[])
  n = size(G, 1)
  if length(S) == 0
    S = 1:n
  end
  ns = length(S)
  result = zeros(Float64, ns , 1)
  H = G[S, S]
  (rr, cc) = findn(H)
  ne = length(rr)
  for i=1:ne
    crow = rr[i]
    ccol = cc[i]
    result[crow] = result[crow] + H[crow, ccol]
  end
  return result
end

function degree{T}(G::SparseMatrixCSC{T, Int64}, S::Array{Int64, 1}=Int64[]; is_directed=false, direction="both" )
  n = size(G, 1)
  if length(S) == 0
    S = 1:n
  end
  ns = length(S)
  result = zeros(Float64, ns , 1)
  H = G[S, S]
  (rr, cc) = findn(H)
  ne = length(rr)
  @inbounds for i=1:ne
    crow = rr[i]
    ccol = cc[i]
    cval = H[crow, ccol]
    if is_directed==false
      result[crow] = result[crow] + cval
    elseif direction == "both"
      result[crow] = result[crow] + cval
      result[ccol] = result[ccol] + cval
    elseif direction == "in"
      result[ccol] = result[ccol] + cval
    elseif direction == "out"
      result[crow] = result[crow] + cval
    else
      error("Invalid direction option!")
    end
  end
  return result
end

function subgraph_density{T}(G::SparseMatrixCSC{T, Int64}, S::Array{Int64, 1}=Int64[]; is_directed=false)
  n = size(G, 1)
  ns = length(S)
  result = 0.0::Float64
  if ns == 0 # return the graph density
    result= sum(G)/n
  else
    @assert minimum(S) >= 1
    @assert maximum(S) <= n
    result = sum(G[S, S])/ length(S)
  end

  if is_directed == false
    result /= 2.0
  end
  return result
end


function subgraph_relative_density{T}(G::SparseMatrixCSC{T, Int64}, S::Array{Int64, 1}=Int64[]; is_directed=false)
  relative=0.0
  deg = subgraph_density(G, S; is_directed=is_directed)
  n = size(G, 1)
  ns = length(S)
  if ns == 0
    relative = deg/(n-1)
  else
    relative = deg/(ns - 1)
  end

  if is_directed == false
    relative = relative * 2.0
  end

  return relative
end


function common_subgraph_density{T}(G::Array{SparseMatrixCSC{T, Int64}, 1}, S::Array{Int64, 1}; is_directed=false)
  min_density = 0.0
  m = length(G)
  densities = Float64[]
  for j=1:m
    curr_density = subgraph_density(G[j], S; is_directed=is_directed)
    push!(densities, curr_density)
  end
  (min_density, min_idx) = findmin(densities)
  return min_density
end


function dcs_intersection_graph{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
  m = length(H)
  n = size(H[1], 1)
  G = spzeros(T, n, n)

  ii, jj = findn(H[1])
  for k=1:length(ii)
    cr = ii[k]
    cc = jj[k]
    G[cr, cc] = H[1][cr, cc]
    for j=2:m
      G[cr, cc] = min(G[cr, cc], H[j][cr, cc])
    end
  end
  return G
end

function dcs_union_graph{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
  m = length(H)
  n = size(H[1], 1)
  G = spzeros(T, n, n)

  ii, jj = findn(H[1])
  for k=1:length(ii)
    cr = ii[k]
    cc = jj[k]
    G[cr, cc] = H[1][cr, cc]
    for j=2:m
      G[cr, cc] = max(G[cr, cc], H[j][cr, cc])
    end
  end
  return G
end

function dcs_avg_graph{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
  m = length(H)
  n = size(H[1], 1)
  G = spzeros(Float64, n, n)

  for j=1:m
    ii, jj= findn(H[j])
    for k=1:length(ii)
      cr = ii[k]
      cc = jj[k]
      G[cr, cc] += H[j][cr, cc]/m
    end
  end
  return G
end
