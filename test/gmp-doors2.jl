include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots
using .PWA

c(x,u) = x'*x + u'*u + 1
A = [0 1; 0 0]
B = [0; 1]
l = let x = Symbolics.variables(:x, 1:2)
    k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x)
    d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x)
    room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x)
    lf = G(room) & U(!d1,k1)# & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)
end
x0 = [-0.4, 0.0]
xT = [-0.0, 0.0]

s, q0, qT = PPWA(A, B, l)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")
savefig("img/gmp-doors2.pdf")
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
