using Symbolics, LazySets
using Gurobi
using Plots
include("pwa/product.jl")
include("pwa/miqp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
      lane1 = HalfSpace(x[2] >= -3.0, x) & HalfSpace(x[2] <= -1.0, x) & HalfSpace(x[1] <= -5.0, x)
      lane2 = HalfSpace(x[2] >= -1.0, x) & HalfSpace(x[2] <=  1.0, x)
      G(lane1 | lane2)
end

h, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(h), " modes")
println(HybridSystems.ntransitions(h), " transitions")

c(x,u) = sum(x.^2) + sum(u.^2)
policy = MIQPPolicy(h, c, h=0.1, T=50, optimizer=Gurobi.Optimizer)

x0 = [-10.0, -2.0, 10.0, 0.0]
xT = [-0.0, -0.0, 0.0, 0.0]
(u, (x, q), m) = action(policy, (q0, x0), (qT, xT))
scatter(x[:,1], x[:,2], label=objective_value(m))
savefig("img/miqp-merge.pdf")