using Symbolics, LazySets
using Gurobi
using Plots
include("pwa/product.jl")
include("pwa/miqp.jl")
include("pwa/simulate.jl")

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

up, (xp, qp), m = action(policy, (q0, x0), (qT, xT))

x, q, u = simulate(
    EulerSimulator(policy.h, 50, Gurobi.Optimizer),
    h,
    ((q0,x0),) -> let (u, _, _) = action(policy, (argmax(q0), x0), (qT, xT)); u[1,:] end,
    (qp[1,:], x0)
)

scatter(xp[:,1], xp[:,2], label=objective_value(m))
scatter!(x[:,1], x[:,2], label="closed-loop")
savefig("img/miqp.pdf")
