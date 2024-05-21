using Symbolics, LazySets
using HybridSystems
using Gurobi
using Ipopt
using Plots
include("../pwa/product.jl")
include("../pwa/miqp.jl")
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
policy = MIQPPolicy(h, c, h=0.1, T=50, optimizer=Gurobi.Optimizer)
up, (xp, qp), m = action(policy, (q0, x0), (qT, xT))
P = argmax.(eachrow(qp))
P = [p for (i,p) in enumerate(P) if i == 1 || p != P[i-1]]
@show P

# P = [1, 5, 2, 8]

qpolicy = QCQPPolicy(h, c; T=10, optimizer=Ipopt.Optimizer)
uq, (xq, qq), m = action(qpolicy, (P[1], x0), (P[end], xT), P)
scatter(xq[:,1],xq[:,2])
savefig("test/qcqp.pdf")
