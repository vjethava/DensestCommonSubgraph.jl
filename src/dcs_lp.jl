Pkg.add("Convex")
Pkg.add("Gurobi")
Pkg.add("Mosek")

include("dcs.jl")

using Convex
using Mosek 

set_default_solver(MosekSolver())

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

H = test_dcs(10);
G = dcs_convert(H);

function dcs_lp{T}(G::Array{SparseMatrixCSC{T, Int64}, 1})
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
        ## trace_val = 0.0
        ## for k=1:n
        ##    for l=1:n
        ##       trace_val = trace_val + X[j][k, l] * G[j][k, l]
        ##    end
        ## end
        push!(constraints, trace(G[j] * X[j]) >= t)
        ii, jj = findn(G[j])
        for k=1:length(ii)
            crow = ii[k]
            ccol = jj[k]
            push!(constraints, X[j][crow, ccol] <= Y[crow])
            push!(constraints, X[j][crow, ccol] <= Y[ccol])
        end
    end

    problem  = maximize(t, constraints)

end



