using Symbolics, LazySets
using Gurobi
using Plots
include("pwa/product.jl")
include("pwa/miqp.jl")
include("pwa/simulate.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
      k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.4 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
      d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.2 <= x[2], x) & HalfSpace(x[2] <= 0.3, x)
      k2 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.1, x)
      d2 = HalfSpace(-0.2 <= x[1], x) & HalfSpace(x[1] <= -0.1, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
      w1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
      w2 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
      w3 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
      room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
      G(room) & U(!d1,k1)# & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)
end

h, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(h), " modes")
println(HybridSystems.ntransitions(h), " transitions")

c(x,u) = sum(x.^2) + sum(u.^2)
policy = MIQPPolicy(h, c, h=0.3, T=30, optimizer=Gurobi.Optimizer)

x0 = [-0.4, 0.4, 0.0, 0.0]
xT = [-0.0, 0.0, 0.0, 0.0]

up, (xp, qp), m = action(policy, (q0, x0), (qT, xT))

x, q, u = simulate(
    EulerSimulator(0.1, 50, Gurobi.Optimizer),
    h,
    ((q0,x0),) -> let (u, _, _) = action(policy, (argmax(q0), x0), (qT, xT)); u[1,:] end,
    (qp[1,:], x0)
)

scatter(xp[:,1], xp[:,2], label=objective_value(m))
scatter!(x[:,1], x[:,2], label="closed-loop")
savefig("img/miqp-doors.pdf")
