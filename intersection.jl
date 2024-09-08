include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots; ENV["GKSwstype"] = "100"
using .PWA
using Interpolations

c(x,u) = x'*x + u'*u
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    r = HalfSpace(-1.0 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(-1.0 <= x[2], x) & HalfSpace(x[2] <= -0.0, x)
    d = HalfSpace(-0.7 <= x[1], x) & HalfSpace(x[1] <= -0.3, x) & HalfSpace(-0.7 <= x[2], x) & HalfSpace(x[2] <= -0.3, x)
    e = HalfSpace(-0.3 <= x[1], x) & HalfSpace(-0.3 <= x[2], x)
    G(r) & G(!d) & F(e)
end
x0 = [-1.0, -0.9, 0.0, 0.0]
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
tt = tq[1]:0.04:tq[end]
scatter(x1.(tt),x2.(tt),label=false,markersize=2)
plot!([-0.7, -0.3, -0.3, -0.7, -0.7], [-0.7, -0.7, -0.3, -0.3, -0.7], linestyle=:dash, color=:orange, fill=true, fillalpha=0.2, label=false)
plot!([-1.0, -0.7, -0.7], [-0.7, -0.7, -1.0], linestyle=:dash, color=:black, label=false)
plot!([-0.3, -0.3, -0.0], [-1.0, -0.7, -0.7], linestyle=:dash, color=:black, label=false)
plot!([-1.0, -0.7, -0.7], [-0.3, -0.3, -0.0], linestyle=:dash, color=:black, label=false)
plot!([-0.3, -0.3, -0.0], [-0.0, -0.3, -0.3], linestyle=:dash, color=:black, label=false)
plot!([-1.0, -0.0, -0.0, -1.0, -1.0], [-1.0, -1.0, -0.0, -0.0, -1.0], linestyle=:dash, color=:black, label=false)
plot!(ratio=:equal)
savefig("img/intersection.pdf")
