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

h, q0, qT = PPWA(A, B, l)
policy = MIQPPolicy(h, c, h=0.3, T=30, optimizer=Gurobi.Optimizer)
up, (xp, qp), m = action(policy, (q0, x0), (qT, xT))

scatter(xp[:,1], xp[:,2], label=objective_value(m))
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.0, 0.0, 0.2, 0.2, 0.0], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.3, 0.3, 0.5, 0.5, 0.3], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.1, -0.0, -0.0, -0.1, -0.1], [0.1, 0.1, 0.5, 0.5, 0.1], color=:black, fill=true, fillalpha=0.2, label=false)
plot!([-0.6, -0.5, -0.5, -0.6, -0.6], [0.4, 0.4, 0.5, 0.5, 0.4], color=:green, linestyle=:dash, label=false)
plot!([-0.6, -0.5, -0.5, -0.6, -0.6], [0.0, 0.0, 0.1, 0.1, 0.0], color=:green, linestyle=:dash, label=false)
plot!([-0.3, -0.2, -0.2, -0.3, -0.3], [0.2, 0.2, 0.3, 0.3, 0.2], color=:red, linestyle=:dash, label=false)
plot!([-0.2, -0.1, -0.1, -0.2, -0.2], [0.0, 0.0, 0.2, 0.2, 0.0], color=:red, linestyle=:dash, label=false)
plot!([-0.1, -0.0, -0.0, -0.1, -0.1], [0.0, 0.0, 0.1, 0.1, 0.0], color=:blue, linestyle=:dash, label=false)
plot!([-0.6, -0.0, -0.0, -0.6, -0.6], [0.0, 0.0, 0.5, 0.5, 0.0], color=:black, label=false)
savefig("img/miqp-doors.pdf")
println(HybridSystems.nmodes(h), " modes")
println(HybridSystems.ntransitions(h), " transitions")
