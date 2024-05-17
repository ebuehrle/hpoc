using Symbolics, LazySets
using Gurobi
using Plots
include("pwa/product.jl")
include("pwa/miqp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))
end

h, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(h), " modes")
println(HybridSystems.ntransitions(h), " transitions")

c(x,u) = sum(x.^2) + sum(u.^2)
policy = MIQPPolicy(h, c, h=0.1, T=50, optimizer=Gurobi.Optimizer)

x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]
(u, (x, q), m) = action(policy, (q0, x0), (qT, xT))
scatter(x[:,1], x[:,2], label=objective_value(m))
savefig("img/miqp.pdf")
