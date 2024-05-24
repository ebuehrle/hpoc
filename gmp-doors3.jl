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
    k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.4 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.2 <= x[2], x) & HalfSpace(x[2] <= 0.3, x)
    k2 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.1, x)
    d2 = HalfSpace(-0.2 <= x[1], x) & HalfSpace(x[1] <= -0.1, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
    w1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
    w2 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    w3 = HalfSpace(-0.1 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.1 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    lf = G(room) & U(!d1,k1) & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)
end
x0 = [-0.41, 0.41, 0.0, 0.0]
xT = [-0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

c(x,u) = x'*x + u'*u + 0.1
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
C, p, E, m = _action(policy, (q0, x0), (qT, xT))
P = decode(E, log.(clamp.(p, 1e-6, 1-1e-6)), nmodes(s)+1, nmodes(s)+2)
P = P[2:end-1]
println(P)
[println(HybridSystems.mode(s,p).X.constraints) for p in P]

qpolicy = QCQPPolicy(s, c; T=20, optimizer=Ipopt.Optimizer)
uq, (xq, qq), m2 = action(qpolicy, (P[1], x0), (P[end], xT), P)
scatter(xq[:,1],xq[:,2],label="$(round(objective_value(m), digits=2)) / $(round(objective_value(m2), digits=2))")
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.0, 0.0, 0.2, 0.2, 0.0], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.3, 0.3, 0.5, 0.5, 0.3], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.1, -0.0, -0.0, -0.1, -0.1], [0.1, 0.1, 0.5, 0.5, 0.1], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.6, -0.5, -0.5, -0.6, -0.6], [0.4, 0.4, 0.5, 0.5, 0.4], color=:green, linestyle=:dash, label=false)
plot!([-0.6, -0.5, -0.5, -0.6, -0.6], [0.0, 0.0, 0.1, 0.1, 0.0], color=:green, linestyle=:dash, label=false)
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.2, 0.2, 0.3, 0.3, 0.2], color=:red, linestyle=:dash, label=false)
plot!([-0.2, -0.1, -0.1, -0.2, -0.2], [0.0, 0.0, 0.2, 0.2, 0.0], color=:red, linestyle=:dash, label=false)
plot!([-0.1, -0.0, -0.0, -0.1, -0.1], [0.0, 0.0, 0.1, 0.1, 0.0], color=:blue, linestyle=:dash, label=false)
plot!([-0.6, -0.0, -0.0, -0.6, -0.6], [0.0, 0.0, 0.5, 0.5, 0.0], color=:black, label=false)
savefig("img/gmp-doors3.pdf")
write("img/gmp-doors3.txt","$(objective_value(m))")
write("img/gmpp-doors3.txt","$(p)")
