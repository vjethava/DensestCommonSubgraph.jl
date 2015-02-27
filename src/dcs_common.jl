export subgraph_density, subgraph_relative_density, common_subgraph_density

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
  return min_density, min_idx
end
