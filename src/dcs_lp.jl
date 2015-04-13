export dcs_lp, charikar_lp, dcs_lp_extract, extract_primal_solution, charikar_JuMP_LP
export LP_METHODS

using JuMP
using Gurobi
using Convex
LP_METHODS=["DCS"; "MIN"; "AVG"]




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


function old_dcs_lp{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
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
  return problem, X, Y, t
end


function dcs_lp{T}(G::Array{SparseMatrixCSC{T, Int64}, 1})

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
    w = sum(Y)
    constraints = Convex.Constraint[]
    push!(constraints, w <= 1.0)


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

  Z = SparseMatrixCSC{Float64, Int64}[]
  typeof(Z)
  for j=1:m
    push!(Z, spzeros(Float64, n, n))
  end
  for j=1:m
    for k=1:nz[j]
      crow = iv[j][k]
      ccol = jv[j][k]
      ci = ni[j] - 1 + k
      # @assert crow >= 1 && crow <= n
      # @assert ccol >= 1 && ccol <= n
      # @assert G[j][crow, ccol] != 0
      # println(X.value[ci])
      Z[j][crow, ccol] = X.value[ci]
    end
  end

  # println("Optval: ", problem.optval)
  return problem, Z, Y.value, t.value
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
  return problem, X, Y, t
end

# function DCSP(G::Graph[])
#   m = length(G)
#   n = size(G[1], 1)
#   R = Model(solver=GurobiSolver(OutputFlag=1))
#   @defVar(R, Y[1:n] >= 0)
#   for i=1:m
#     [gi, gj] = findn(G[i])
#     ne = length(gi)
#     @defVar(R, X[i=i, j=1:ne] >= 0)

# end

function charikar_JuMP_LP(G::Graph)
  n = size(G, 1)
  tic()
  ii, jj = findn(G)
  ne = length(ii)
  println("graph with $n nodes and $ne edges. Setting up Charikar LP ... ")
  R = Model(solver=GurobiSolver(OutputFlag=0))
  @defVar(R, X[1:ne] >= 0)
  @defVar(R, Y[1:n] >= 0)
  @addConstraint(R, gamma, sum{Y[i], i=1:n} <= 1)

  @addConstraint(R, alpha[k=1:ne], X[k] <= Y[ ii[k] ] )
  @addConstraint(R, beta[k=1:ne], X[k] <= Y[ jj[k] ] )
#   for k=1:ne
#     crow = ii[k]
#     ccol = jj[k]
#     @addConstraint(R, X[crow, ccol] <= Y[crow])
#     @addConstraint(R, X[crow, ccol] <= Y[ccol])
#   end
  @setObjective(R, Max, sum{X[k], k=1:ne })
  t_elapsed = toc();
  println("Completed setup, starting to solve LP with Gurobi ")
  status = solve(R);
  if status == :Optimal
    println("solved LP successfully, returning primal variables")
  end

  t = R.objVal
  X1 = spzeros(n, n)
  Y1 = spzeros(n, 1)
  for k=1:ne
    X1[ii[k], jj[k]]= getValue(X[k])
  end
  for j=1:n
    Y1[j] = getValue(Y[j])
  end
  return X1, Y1, t, status
end


function charikar_lp{T}(G::SparseMatrixCSC{T, Int64})

  set_default_solver(GurobiSolver())
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
  t= problem.optval
  return problem, X.value, Y.value, t
end


function extract_primal_solution(Y)
  V = [j > 0.0 for j in Y]
  S = find(V)
  return S
end

# Uses LP primal solution to find dense common subgraph
#
#
function dcs_lp_extract(G::GraphVec; method::String = "DCS" )
  m = length(G)
  @assert m >= 1
  n = size(G{1}, 1)

  if method == "DCS"
    p, X, Y, t = dcs_lp(G)
  elseif method == "MIN"
    H = dcs_intersection_graph(G)
    p, X, Y, t =  charikar_lp(H)
  elseif method == "AVG"
    H = dcs_avg_graph(G)
    p, X, Y, t = charikar_lp(H)
  else
    error("Undefined method: $method") #  ('dcs', 'charikar', 'weighted')"
  end
  S = extract_primal_solution(Y)
  d = common_subgraph_density(S)
  return (S, d, t)
end


# function dcs_charikar_lp{T}(H::Array{SparseMatrixCSC{T, Int64}, 1})
#   m = length(H)
#   n = size(H[1], 1)
#   G = spzeros(T, n, n)

#   ii, jj = findn(H[1])
#   for k=1:length(ii)
#     cr = ii[k]
#     cc = jj[k]
#     G[cr, cc] = H[1][cr, cc]
#     for j=2:m
#       G[cr, cc] = min(G[cr, cc], H[j][cr, cc])
#     end
#   end
#   return charikar_lp(G)

  # end
