include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using Gurobi
using Plots
using .PWA

c(x,u) = x'*x + u'*u
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
      lane1 = HalfSpace(x[2] >= -3.0, x) & HalfSpace(x[2] <= -1.0, x) & HalfSpace(x[1] <= -5.0, x)
      lane2 = HalfSpace(x[2] >= -1.0, x) & HalfSpace(x[2] <=  1.0, x)
      G(lane1 | lane2)
end
x0 = [-10.0, -2.0, 10.0, 0.0]
xT = [-0.0, -0.0, 0.0, 0.0]

h, q0, qT = PPWA(A, B, l)
policy = MIQPPolicy(h, c, h=0.3, T=30, optimizer=Gurobi.Optimizer)
up, (xp, qp), m = action(policy, (q0, x0), (qT, xT))

scatter(xp[:,1], xp[:,2], label=objective_value(m))
savefig("img/miqp-merge.pdf")
