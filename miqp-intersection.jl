using Symbolics, LazySets
using Gurobi
using Plots
include("pwa/product.jl")
include("pwa/miqp.jl")
include("pwa/simulate.jl")

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
policy = MIQPPolicy(s, c, h=0.1, T=50, optimizer=Gurobi.Optimizer)
up, (xp, qp), m = action(policy, (q0, x0), (qT, xT))

scatter(xp[:,1], zeros(size(xp[:,1])),label=false)
scatter!(zeros(size(xp[:,2])), xp[:,2],label=false)
plot!([-1.0, -0.3, -0.3], [-0.1, -0.1, -1.0], color=:black,label=false)
plot!([0.1, 0.1, 1.0], [-1.0, -0.1, -0.1], color=:black,label=false)
plot!([-1.0, -0.3, -0.3], [0.3, 0.3, 1.0], color=:black,label=false)
plot!([0.1, 0.1, 1.0], [1.0, 0.3, 0.3], color=:black,label=false)
plot!([-1.0, -0.3, -0.3], [0.1, 0.1, -0.1], color=:black, linestyle=:dash,label=false)
plot!([0.1, 0.1, 1.0], [0.3, 0.1, 0.1], color=:black, linestyle=:dash,label=false)
plot!([-0.1, -0.1, 0.1], [-1.0, -0.1, -0.1], color=:black, linestyle=:dash,label=false)
plot!([-0.1, -0.1, -0.3], [1.0, 0.3, 0.3], color=:black, linestyle=:dash,label=false)
plot!(background_color_inside=:lightgrey)
savefig("img/miqp-intersection.pdf")
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
