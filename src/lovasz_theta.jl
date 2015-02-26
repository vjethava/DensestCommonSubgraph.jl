# This files implements Lovasz Theta computation described in the following paper:
#
# An SDP Primal-Dual Algorithm for Approximating the Lovász-Theta Function
# T.-H. Hubert Chan, Kevin L. Chang and Rajiv Raman
# Algorithmica (2014) 69:605-618
# DOI: 10.1007/s00453-013-9756-5
#
using Base.Test
Pkg.add("Graphs")
Pkg.add("SVM")


using SVM

function test_expm()
  a = 2e10
  b = 4e8/6
  c = 200/3
  d = 3
  e = 1e-8
  A = [0 e 0; -(a+b) -d a; c 0 -c]
  E = expm(A)
  E_true = [  0.446849 1.54044e-9 0.462811; -5.74307e6 -0.015283 -4.52654e6;    0.447723 1.5427e-9 0.463481]
  Z = E_true ./ E # test relative error
  @test_approx_eq_eps ones(3,3) Z 1e-3
end
test_expm()


function ORACLE(W, β, E)
  n = size(W, 1)
  Ø = zeros(n, n)
  if isequal(W, Ø)
    z = β
    Y = Ø
    return ("dual", z, Y)
  else
    J = ones(n, n)
    t = trace(W)
    X = W ./ t
    σ = sum (J .* W) - β
    if σ <= 0.0
      z = β
      Y = Ø
      return ("dual", z, Y)
    end
    X_E = X .* E
    f = normfro(X_E)
    if f >= σ/(β*n)
      z = β
      Y = σ/(f^2) * X_E
      return ("dual", z, Y)
    end
    X_tilde = X - X_E + σ/(β*n)*J
    X_hat  = X_tilde ./ trace(X_tilde)
    return ("primal", null, X_hat)
  end
end

function ϑ(E, ρ, δ, β)
  n = size(E, 1)
  J = ones(n, n)
  I = eye(n)
  S_M = zeros(n, n)
  S_Y = zeros(n, n)
  W = I
  ϵ = δ/(2*ρ)
  ϵ_1 = -log(1 - ϵ)
  T = iceil ( 8 * ρ^2 * log(n) /  δ^2 )
  for t=1:T
    (solution, val, mat) = ORACLE(W, β, E)
    if solution == "primal"
      return (solution, val, mat)
    end
    Y = sum( mat .* E)
    M = ( β*I + Y - J + ρ*I ) / (2*ρ)
    S_M += M
    S_Y += Y
    print(W)
    W = expm(-ϵ_1 * S_M)
  end
  z = β + δ
  Y_bar = 1\T * S_Y
  return ("dual", z, Y_bar)
end


using Graphs

# function test_ϑ()
  n = 10
  p = 0.5
  G = Graphs.erdos_renyi_graph(n, p)
  E = int(Graphs.adjacency_matrix(G))
  β = n
  ρ = β*n
  δ = 0.5
  W = eye(n, n)
  (sol, val, mat) = ORACLE(W, β, E)
  print(sol)
(solution, val, mat) = ϑ(E, ρ, δ, β)
# end


test_ϑ()
