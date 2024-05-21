using Symbolics, LazySets
using HybridSystems
using MosekTools
using Ipopt
using Plots
include("pwa/product.jl")
include("pwa/gmp.jl")
include("pwa/qcqp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    lane1 = HalfSpace(x[2] >= -3.0, x) & HalfSpace(x[2] <= -1.0, x) & HalfSpace(x[1] <= -5.0, x)
    lane2 = HalfSpace(x[2] >= -1.0, x) & HalfSpace(x[2] <=  1.0, x)
    l = G(lane1 | lane2)
end
x0 = [-10.0, -2.0, 10.0, 0.0]
xT = [-0.0, -0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

c(x,u) = sum(x.^2) + sum(u.^2)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
C, p, E, m = action(policy, (q0, x0), (qT, xT))
P = decode(E, log.(p.+1e-6), nmodes(s)+1, nmodes(s)+2)
P = P[2:end-1]

qpolicy = QCQPPolicy(s, c; T=20, optimizer=Ipopt.Optimizer)
uq, (xq, qq), m = action(qpolicy, (P[1], x0), (P[end], xT), P)
scatter(xq[:,1],xq[:,2])
savefig("img/gmp-merge.pdf")
write("img/gmp-merge.txt","$(objective_value(m))")
write("img/gmpp-merge.txt","$(p)")
