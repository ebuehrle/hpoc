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
    G(!(HalfSpace(x[1] >= -0.3, x) & HalfSpace(x[1] <= 0.1, x) 
    & HalfSpace(x[2] >= -0.1, x) & HalfSpace(x[2] <= 0.3, x)))
end
x0 = [-1.0, -0.8, 0.5, 0.5]
xT = [ 1.0,  1.0, 0.5, 0.5]

s, q0, qT = PPWA(A, B, l)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

scatter(xq[:,1], zeros(size(xq[:,1])),label=false)
scatter!(zeros(size(xq[:,2])), xq[:,2],label=false)
plot!([-1.0, -0.3, -0.3], [-0.1, -0.1, -1.0], color=:black,label=false)
plot!([0.1, 0.1, 1.0], [-1.0, -0.1, -0.1], color=:black,label=false)
plot!([-1.0, -0.3, -0.3], [0.3, 0.3, 1.0], color=:black,label=false)
plot!([0.1, 0.1, 1.0], [1.0, 0.3, 0.3], color=:black,label=false)
plot!([-1.0, -0.3, -0.3], [0.1, 0.1, -0.1], color=:black, linestyle=:dash,label=false)
plot!([0.1, 0.1, 1.0], [0.3, 0.1, 0.1], color=:black, linestyle=:dash,label=false)
plot!([-0.1, -0.1, 0.1], [-1.0, -0.1, -0.1], color=:black, linestyle=:dash,label=false)
plot!([-0.1, -0.1, -0.3], [1.0, 0.3, 0.3], color=:black, linestyle=:dash,label=false)
plot!(background_color_inside=:lightgrey)
savefig("img/gmp-intersection.pdf")
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
