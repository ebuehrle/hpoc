include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots
using .PWA
using Interpolations

c(x,u) = x'*x + u'*u
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    y = HalfSpace(-0.9 <= x[1], x) & HalfSpace(x[1] <= -0.8, x) & HalfSpace(-0.6 <= x[2], x) & HalfSpace(x[2] <= -0.4, x)
    b = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.4, x) & HalfSpace(-0.6 <= x[2], x) & HalfSpace(x[2] <= -0.4, x)
    r = HalfSpace(-0.5 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(-0.9 <= x[2], x) & HalfSpace(x[2] <= -0.8, x)
    g = HalfSpace(-0.5 <= x[1], x) & HalfSpace(x[1] <= -0.25, x) & HalfSpace(-0.25 <= x[2], x) & HalfSpace(x[2] <= -0.0, x)
    F(y) & G(!g) & G(!b)
end
x0 = [-0.9, -0.9, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, tq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

x1 = linear_interpolation(tq, xq[:,1])
x2 = linear_interpolation(tq, xq[:,2])
tt = tq[1]:0.02:tq[end]
scatter(x1.(tt),x2.(tt),label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")
plot!([-0.9, -0.8, -0.8, -0.9, -0.9], [-0.6, -0.6, -0.4, -0.4, -0.6], linestyle=:dash, color=:yellow, fill=true, fillalpha=0.2, label=false)
plot!([-0.6, -0.4, -0.4, -0.6, -0.6], [-0.6, -0.6, -0.4, -0.4, -0.6], linestyle=:dash, color=:blue, fill=true, fillalpha=0.2, label=false)
plot!([-0.5, -0.2, -0.2, -0.5, -0.5], [-0.9, -0.9, -0.8, -0.8, -0.9], linestyle=:dash, color=:red, fill=true, fillalpha=0.2, label=false)
plot!([-0.5, -0.25, -0.25, -0.5, -0.5], [-0.25, -0.25, -0.0, -0.0, -0.25], linestyle=:dash, color=:green, fill=true, fillalpha=0.2, label=false)
plot!([-1.1,  0.1,  0.1, -1.1, -1.1], [-1.1, -1.1,  0.1,  0.1, -1.1], color=:black, label=false)
savefig("img/stlcg-2.pdf")
