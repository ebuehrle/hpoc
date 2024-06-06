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
    k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.4 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.2 <= x[2], x) & HalfSpace(x[2] <= 0.3, x)
    k2 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.1, x)
    d2 = HalfSpace(-0.2 <= x[1], x) & HalfSpace(x[1] <= -0.1, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
    w1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
    w2 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    w3 = HalfSpace(-0.1 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.1 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    G(room) & U(!d1,k1) & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)
end
x0 = [-0.41, 0.41, 0.0, 0.0]
xT = [-0.0, 0.0, 0.0, 0.0]
x02 = [-0.41, 0.19, 0.0, 0.0]
xT2 = [-0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)
uq2, xq2, qq2, mq2, m2 = extract(policy, (q0, x02), (qT, xT2); T=20, optimizer=Ipopt.Optimizer)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")
scatter!(xq2[:,1],xq2[:,2],label="J = $(round(objective_value(mq2), digits=2)) ($(round(objective_value(m2), digits=2)))")
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.0, 0.0, 0.2, 0.2, 0.0], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.3, 0.3, 0.5, 0.5, 0.3], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.1, -0.0, -0.0, -0.1, -0.1], [0.1, 0.1, 0.5, 0.5, 0.1], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.6, -0.5, -0.5, -0.6, -0.6], [0.4, 0.4, 0.5, 0.5, 0.4], color=:green, linestyle=:dash, label=false)
plot!([-0.6, -0.5, -0.5, -0.6, -0.6], [0.0, 0.0, 0.1, 0.1, 0.0], color=:green, linestyle=:dash, label=false)
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.2, 0.2, 0.3, 0.3, 0.2], color=:red, linestyle=:dash, label=false)
plot!([-0.2, -0.1, -0.1, -0.2, -0.2], [0.0, 0.0, 0.2, 0.2, 0.0], color=:red, linestyle=:dash, label=false)
plot!([-0.1, -0.0, -0.0, -0.1, -0.1], [0.0, 0.0, 0.1, 0.1, 0.0], color=:blue, linestyle=:dash, label=false)
plot!([-0.6, -0.0, -0.0, -0.6, -0.6], [0.0, 0.0, 0.5, 0.5, 0.0], color=:black, label=false)
savefig("img/gmp-doors.pdf")
