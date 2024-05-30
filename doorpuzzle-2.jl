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
    k1 = And_([HalfSpace(0.21-1 <= x[1], x), HalfSpace(x[1] <= 0.29-1, x), HalfSpace( 0.40 <= x[2], x), HalfSpace(x[2] <=  0.50, x)])
    k2 = And_([HalfSpace(0.11-1 <= x[1], x), HalfSpace(x[1] <= 0.19-1, x), HalfSpace( 0.40 <= x[2], x), HalfSpace(x[2] <=  0.50, x)])
    k3 = And_([HalfSpace(0.00-1 <= x[1], x), HalfSpace(x[1] <= 0.10-1, x), HalfSpace( 0.21 <= x[2], x), HalfSpace(x[2] <=  0.29, x)])
    k4 = And_([HalfSpace(0.00-1 <= x[1], x), HalfSpace(x[1] <= 0.10-1, x), HalfSpace(-0.21 <= x[2], x), HalfSpace(x[2] <= -0.29, x)])
    k5 = And_([HalfSpace(0.11-1 <= x[1], x), HalfSpace(x[1] <= 0.19-1, x), HalfSpace(-0.40 <= x[2], x), HalfSpace(x[2] <= -0.50, x)])
    k6 = And_([HalfSpace(0.21-1 <= x[1], x), HalfSpace(x[1] <= 0.29-1, x), HalfSpace(-0.40 <= x[2], x), HalfSpace(x[2] <= -0.50, x)])

    d1 = And_([HalfSpace(0.88-1 <= x[1], x), HalfSpace(x[1] <= 0.92-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d2 = And_([HalfSpace(0.78-1 <= x[1], x), HalfSpace(x[1] <= 0.82-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d3 = And_([HalfSpace(0.68-1 <= x[1], x), HalfSpace(x[1] <= 0.72-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d4 = And_([HalfSpace(0.58-1 <= x[1], x), HalfSpace(x[1] <= 0.62-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d5 = And_([HalfSpace(0.48-1 <= x[1], x), HalfSpace(x[1] <= 0.52-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d6 = And_([HalfSpace(0.38-1 <= x[1], x), HalfSpace(x[1] <= 0.42-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])

    room = And_([HalfSpace(-1.00 <= x[1], x), HalfSpace(x[1] <= 0.00, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])

    G(room) & U(!d1,k1) & U(!d2,k2) & U(!d3,k3) & U(!d4,k4) & U(!d5,k5) & U(!d6,k6)
end
x0 = [-0.8, 0.0, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")

plot!([-0.79, -0.71, -0.71, -0.79, -0.79], [ 0.40,  0.50,   0.50,  0.40,  0.40], color=:green, linestyle=:dash, label=false)
plot!([-0.89, -0.81, -0.81, -0.89, -0.89], [ 0.40,  0.50,   0.50,  0.40,  0.40], color=:green, linestyle=:dash, label=false)
plot!([-1.00, -0.90, -0.90, -1.00, -1.00], [ 0.21,  0.29,   0.29,  0.21,  0.21], color=:green, linestyle=:dash, label=false)
plot!([-1.00, -0.90, -0.90, -1.00, -1.00], [-0.21, -0.29,  -0.29, -0.21, -0.21], color=:green, linestyle=:dash, label=false)
plot!([-0.89, -0.81, -0.81, -0.89, -0.89], [-0.40, -0.50,  -0.50, -0.40, -0.40], color=:green, linestyle=:dash, label=false)
plot!([-0.79, -0.71, -0.71, -0.79, -0.79], [-0.40, -0.50,  -0.50, -0.40, -0.40], color=:green, linestyle=:dash, label=false)

plot!([-0.12, -0.08, -0.08, -0.12, -0.12], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, linestyle=:dash, label=false)
plot!([-0.22, -0.18, -0.18, -0.22, -0.22], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, linestyle=:dash, label=false)
plot!([-0.32, -0.28, -0.28, -0.32, -0.32], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, linestyle=:dash, label=false)
plot!([-0.42, -0.38, -0.38, -0.42, -0.42], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, linestyle=:dash, label=false)
plot!([-0.52, -0.48, -0.48, -0.52, -0.52], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, linestyle=:dash, label=false)
plot!([-0.62, -0.58, -0.58, -0.62, -0.62], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, linestyle=:dash, label=false)

plot!([0.0, 1.0, 1.0, 0.0, 0.0], [-1.0, -1.0, 0.0, 0.0, -1.0], color=:black, label=false)

savefig("img/doorpuzzle-2.pdf")
