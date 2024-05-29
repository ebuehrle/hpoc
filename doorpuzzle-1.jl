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
    k1 = HalfSpace(0.55 <= x[1], x) & HalfSpace(x[1] <= 0.70, x) & HalfSpace(-0.95 <= x[2], x) & HalfSpace(x[2] <= -0.80, x)
    k2 = HalfSpace(0.30 <= x[1], x) & HalfSpace(x[1] <= 0.45, x) & HalfSpace(-0.70 <= x[2], x) & HalfSpace(x[2] <= -0.55, x)
    k3 = HalfSpace(0.05 <= x[1], x) & HalfSpace(x[1] <= 0.20, x) & HalfSpace(-0.95 <= x[2], x) & HalfSpace(x[2] <= -0.80, x)
    k4 = HalfSpace(0.80 <= x[1], x) & HalfSpace(x[1] <= 0.95, x) & HalfSpace(-0.95 <= x[2], x) & HalfSpace(x[2] <= -0.80, x)
    k5 = HalfSpace(0.30 <= x[1], x) & HalfSpace(x[1] <= 0.45, x) & HalfSpace(-0.45 <= x[2], x) & HalfSpace(x[2] <= -0.30, x)

    d1 = HalfSpace(0.23 <= x[1], x) & HalfSpace(x[1] <= 0.27, x) & HalfSpace(-1.00 <= x[2], x) & HalfSpace(x[2] <= -0.75, x)
    d2 = HalfSpace(0.73 <= x[1], x) & HalfSpace(x[1] <= 0.77, x) & HalfSpace(-0.50 <= x[2], x) & HalfSpace(x[2] <= -0.25, x)
    d3 = HalfSpace(0.75 <= x[1], x) & HalfSpace(x[1] <= 1.00, x) & HalfSpace(-0.77 <= x[2], x) & HalfSpace(x[2] <= -0.73, x)
    d4 = HalfSpace(0.00 <= x[1], x) & HalfSpace(x[1] <= 0.25, x) & HalfSpace(-0.52 <= x[2], x) & HalfSpace(x[2] <= -0.48, x)
    d5 = HalfSpace(0.75 <= x[1], x) & HalfSpace(x[1] <= 1.00, x) & HalfSpace(-0.27 <= x[2], x) & HalfSpace(x[2] <= -0.23, x)

    w1 = HalfSpace(0.00 <= x[1], x) & HalfSpace(x[1] <= 0.75, x) & HalfSpace(-0.27 <= x[2], x) & HalfSpace(x[2] <= -0.23, x)
    w2 = HalfSpace(0.48 <= x[1], x) & HalfSpace(x[1] <= 0.52, x) & HalfSpace(-0.50 <= x[2], x) & HalfSpace(x[2] <= -0.25, x)
    w3 = HalfSpace(0.25 <= x[1], x) & HalfSpace(x[1] <= 0.50, x) & HalfSpace(-0.52 <= x[2], x) & HalfSpace(x[2] <= -0.48, x)
    w4 = HalfSpace(0.23 <= x[1], x) & HalfSpace(x[1] <= 0.27, x) & HalfSpace(-0.75 <= x[2], x) & HalfSpace(x[2] <= -0.50, x)
    w5 = HalfSpace(0.73 <= x[1], x) & HalfSpace(x[1] <= 0.77, x) & HalfSpace(-1.00 <= x[2], x) & HalfSpace(x[2] <= -0.50, x)

    room = HalfSpace(0.00 <= x[1], x) & HalfSpace(x[1] <= 1.00, x) & HalfSpace(-1.00 <= x[2], x) & HalfSpace(x[2] <= -0.00, x)

    #G(room) & G(!w1) & G(!w2) & G(!w3) & U(!d5,k5) & U(!d4,k4)
    G(room) & G(!w1) & G(!w2) & G(!w3) & G(!w4) & U(!d5,k5) & U(!d4,k4)
    #U(!d1,k1) & U(!d2,k2) & U(!d3,k3) & U(!d4,k4) & U(!d5,k5)
end
x0 = [0.55, -0.75, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=true)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")
plot!([0.55, 0.7, 0.7, 0.55, 0.55], [-0.95, -0.95, -0.8, -0.8, -0.95],  color=:green, linestyle=:dash, label=false)
plot!([0.3, 0.45, 0.45, 0.3, 0.3], [-0.7, -0.7, -0.55, -0.55, -0.7],    color=:green, linestyle=:dash, label=false)
plot!([0.05, 0.2, 0.2, 0.05, 0.05], [-0.95, -0.95, -0.8, -0.8, -0.95],  color=:green, linestyle=:dash, label=false)
plot!([0.8, 0.95, 0.95, 0.8, 0.8], [-0.95, -0.95, -0.8, -0.8, -0.95],   color=:green, linestyle=:dash, label=false)
plot!([0.3, 0.45, 0.45, 0.3, 0.3], [-0.45, -0.45, -0.3, -0.3, -0.45],   color=:green, linestyle=:dash, label=false)
# annotate!(0.7, -0.95, text("(1)", :green, :right))
# annotate!(0.45, -0.7, text("(2)", :green, :right))
# annotate!(0.2, -0.95, text("(3)", :green, :right))
# annotate!(0.95, -0.95, text("(4)", :green, :right))
# annotate!(0.45, -0.45, text("(5)", :green, :right))

plot!([0.23, 0.27, 0.27, 0.23, 0.23], [-1.0, -1.0, -0.75, -0.75, -1.0],     color=:red, linestyle=:dash, label=false)
plot!([0.73, 0.77, 0.77, 0.73, 0.73], [-0.5, -0.5, -0.25, -0.25, -0.5],     color=:red, linestyle=:dash, label=false)
plot!([0.75, 1.0, 1.0, 0.75, 0.75], [-0.77, -0.77, -0.73, -0.73, -0.77],    color=:red, linestyle=:dash, label=false)
plot!([0.0, 0.25, 0.25, 0.0, 0.0], [-0.52, -0.52, -0.48, -0.48, -0.52],     color=:red, linestyle=:dash, label=false)
plot!([0.75, 1.0, 1.0, 0.75, 0.75], [-0.27, -0.27, -0.23, -0.23, -0.27],    color=:red, linestyle=:dash, label=false)

plot!([0.0, 0.75, 0.75, 0.0, 0.0], [-0.27, -0.27, -0.23, -0.23, -0.27],     color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.48, 0.52, 0.52, 0.48, 0.48], [-0.5, -0.5, -0.25, -0.25, -0.5],     color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.25, 0.5, 0.5, 0.25, 0.25], [-0.52, -0.52, -0.48, -0.48, -0.52],    color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.23, 0.27, 0.27, 0.23, 0.23], [-0.75, -0.75, -0.5, -0.5, -0.75],    color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.73, 0.77, 0.77, 0.73, 0.73], [-1.0, -1.0, -0.5, -0.5, -1.0],       color=:black, fill=true, fillalpha=0.3, label=false)

plot!([0.0, 1.0, 1.0, 0.0, 0.0], [-1.0, -1.0, 0.0, 0.0, -1.0], color=:black, label=false)

savefig("img/doorpuzzle-1.pdf")
