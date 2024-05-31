include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots
using .PWA

c(x,u) = x'*x + u'*u
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    y = HalfSpace(-0.9 <= x[1], x) & HalfSpace(x[1] <= -0.8, x) & HalfSpace(-0.6 <= x[2], x) & HalfSpace(x[2] <= -0.4, x)
    b = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.4, x) & HalfSpace(-0.6 <= x[2], x) & HalfSpace(x[2] <= -0.4, x)
    r = HalfSpace(-0.5 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(-0.9 <= x[2], x) & HalfSpace(x[2] <= -0.8, x)
    g = HalfSpace(-0.4 <= x[1], x) & HalfSpace(x[1] <= -0.3, x) & HalfSpace(-0.2 <= x[2], x) & HalfSpace(x[2] <= -0.1, x)
    F(r) & F(g) & G(!b)
end
x0 = [-0.9, -0.9, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")
plot!([-0.9, -0.8, -0.8, -0.9, -0.9], [-0.6, -0.6, -0.4, -0.4, -0.6], color=:yellow, fill=true, fillalpha=0.2, label=false)
plot!([-0.6, -0.4, -0.4, -0.6, -0.6], [-0.6, -0.6, -0.4, -0.4, -0.6], color=:blue, fill=true, fillalpha=0.2, label=false)
plot!([-0.5, -0.2, -0.2, -0.5, -0.5], [-0.9, -0.9, -0.8, -0.8, -0.9], color=:red, fill=true, fillalpha=0.2, label=false)
plot!([-0.4, -0.3, -0.3, -0.4, -0.4], [-0.2, -0.2, -0.1, -0.1, -0.2], color=:green, fill=true, fillalpha=0.2, label=false)
plot!([-1.1,  0.1,  0.1, -1.1, -1.1], [-1.1, -1.1,  0.1,  0.1, -1.1], color=:black, label=false)
savefig("img/stlcg-1.pdf")