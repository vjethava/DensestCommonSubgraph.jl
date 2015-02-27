export dcs_lp, charikar_lp, dcs_charikar_lp, fast_dcs_lp, faster_dcs_lp

using Convex
using SCS

set_default_solver(SCSSolver(verbose=0))

function dcs_convert{T}(G::Array{SparseMatrixCSC{T, Int64}, 1})
  m = length(G)
  H = SparseMatrixCSC{Float64, Int64}[]
  for j=1:m
    push!(H, similar(G[j], Float64))
    ii, jj = findn(G[j])
    for k=1:length(ii)
       H[j][ii[k], jj[k] ]  = 1.0 
    end
  end
  return H
end


function dcs_lp{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
  G = dcs_convert(H)
  m = length(G)
  n = size(G[1], 1)
  X = Convex.Variable[]
  for j=1:m
    push!(X, Variable(n, n,Positive()))
  end 
  Y = Variable(n, Positive())
  t = Variable()
  z = sum(Y)
  constraints = Convex.Constraint[]
  push!(constraints, z <= 1.0)
    
    
  for j=1:m
    push!(constraints, trace(G[j] * X[j]) >= t)
    ii, jj = findn(G[j])
    for k=1:length(ii)
      crow = ii[k]
      ccol = jj[k]
      push!(constraints, X[j][crow, ccol] <= Y[crow])
      push!(constraints, X[j][crow, ccol] <= Y[ccol])
    end
  end

  problem  = maximize(t, constraints);
  solve!(problem);
  # println("Optval: ", problem.optval)
  return problem 
end


function fast_dcs_lp{T}(G::Array{SparseMatrixCSC{T, Int64}, 1})
  m = length(G)
  n = size(G[1], 1)
  nz = Int64[]
  @inbounds begin 
    iv = {}
    jv = {} 
    gv = {} 
    for j=1:m
      ii, jj= findn(G[j])
      push!(iv, ii)
      push!(jv, jj)
      push!(nz, length(ii))
      ss = similar(ii, Float64)
      for k=1:length(ii)
        ss[k] = G[j][ii[k], jj[k] ]
      end
      push!(gv, ss)    
    end
    ni = deepcopy(nz)
    insert!(ni, 1, 1)
    ni = cumsum(ni)

    nv = sum(nz)
    X = Convex.Variable(nv, 1 , Positive())


    Y = Variable(n, Positive())
    t = Variable()
    z = sum(Y)
    constraints = Convex.Constraint[]
    push!(constraints, z <= 1.0)
  
  
    for j=1:m
      cb = ni[j]
      ce = ni[j+1]-1
      # println("gv: $(length(gv[j])) ni: [$cb:$ce]")
      # @assert length(gv[j]) == length([cb:ce;])
      # @assert cb >= 1
      # @assert ce <= length(X)
      push!(constraints, sum(X[cb:ce] .* gv[j] ) >= t)
      # @assert length(gv[j]) == nz[j] 
      for k=1:nz[j]
        crow = iv[j][k]
        ccol = jv[j][k]
        ci = ni[j] - 1 + k
        # @assert (ci >= cb) && (ci <= ce)
        push!(constraints, X[ci] <= Y[crow])
        push!(constraints, X[ci] <= Y[ccol])
        # println("j: $j X[$ci] <= Y[$crow], Y[$ccol]")
      end
    end
  end 
  problem  = maximize(t, constraints);
  solve!(problem);
  # println("Optval: ", problem.optval)
  return problem 
end


function faster_dcs_lp{T}(G::Array{SparseMatrixCSC{T, Int64}, 1})
  m = length(G)
  n = size(G[1], 1)
  nz = Int64[]
  @inbounds begin 
    for j=1:m
      cval = nnz(G[j])
      push!(nz, cval)
    end
  
    ni = deepcopy(nz)
    insert!(ni, 1, 1)
    ni = cumsum(ni)

    X = Convex.Variable(sum(nz), 1 , Positive())
    Y = Variable(n, Positive())
    t = Variable()
    z = sum(Y)
    constraints = Convex.Constraint[]
    push!(constraints, z <= 1.0)
  
  
    for j=1:m
      cb = ni[j]
      ce = (ni[j+1]-1)
      gi = find(G[j])
      gval = map((x) -> convert(Float64, x), G[j][gi])
      push!(constraints, sum(X[cb:ce] .* gval) >= t)
      for k=1:nz[j]
        ccol = iceil(gi[k]/n)
        crow = (gi[k] % n)
        if crow == 0
          crow = n
        end
        ci = ni[j] - 1 + k

        push!(constraints, X[ci] <= Y[crow])
        push!(constraints, X[ci] <= Y[ccol])
        # println("j: $j X[$ci] <= Y[$crow], Y[$ccol]")
      end
    end
  end
  problem  = maximize(t, constraints);
  solve!(problem);
  return problem 
end





function charikar_lp{T}(G::SparseMatrixCSC{T, Int64})
  n = size(G, 1)
  X = Variable(n, n, Positive())
  Y = Variable(n, Positive())
  constraints = Convex.Constraint[]

  push!(constraints, sum(Y) <= 1)
  
  ii, jj = findn(G)
  for k=1:length(ii)
    crow = ii[k]
    ccol = jj[k]
    push!(constraints, X[crow, ccol] <= Y[crow])
    push!(constraints, X[crow, ccol] <= Y[ccol])
  end

  problem = maximize(trace(G * X), constraints); 
  solve!(problem);
  return problem 
end

function dcs_charikar_lp{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
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
  return charikar_lp(G)
  
end
