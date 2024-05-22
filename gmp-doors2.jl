using Symbolics, LazySets
using HybridSystems
using MosekTools
using Ipopt
using Plots
include("pwa/product.jl")
include("pwa/gmp.jl")
include("pwa/qcqp.jl")

A = [0 1; 0 0]
B = [0; 1]
l = let x = Symbolics.variables(:x, 1:2)
    k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x)
    d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x)
    room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x)
    lf = G(room) & U(!d1,k1)# & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)
end
ltlf, dict = translate(LTLTranslator(), l)
x0 = [-0.4, 0.0]
xT = [-0.0, 0.0]

s, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

c(x,u) = x'*x + u'*u + 1
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
C, p, E, m = action(policy, (q0, x0), (qT, xT))
P = decode(E, log.(p.+1e-6), nmodes(s)+1, nmodes(s)+2)
P = P[2:end-1]

qpolicy = QCQPPolicy(s, c; T=20, optimizer=Ipopt.Optimizer)
uq, (xq, qq), m = action(qpolicy, (P[1], x0), (P[end], xT), P)
scatter(xq[:,1],xq[:,2])
savefig("img/gmp-doors2.pdf")
write("img/gmp-doors2.txt","$(objective_value(m))")
write("img/gmpp-doors2.txt","$(p)")
