using Symbolics, LazySets
using HybridSystems
using MosekTools
using Ipopt
using Plots
include("../pwa/product.jl")
include("../pwa/gmp.jl")
include("../pwa/qcqp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))
end
x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]

h, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(h), " modes")
println(HybridSystems.ntransitions(h), " transitions")

c(x,u) = sum(x.^2) + sum(u.^2)
policy = GMPPolicy(h, c; optimizer=Mosek.Optimizer)
p, E, m = action(policy, (q0, x0), (qT, xT))
P = decode(E, p, nmodes(h)+1, nmodes(h)+2)
P = P[2:end-1]

qpolicy = QCQPPolicy(h, c; T=10, optimizer=Ipopt.Optimizer)
uq, (xq, qq), m = action(qpolicy, (P[1], x0), (P[end], xT), P)
scatter(xq[:,1],xq[:,2])
savefig("test/qcqp.pdf")
