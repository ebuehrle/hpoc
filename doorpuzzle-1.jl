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
    k1 = And_([HalfSpace(0.52 <= x[1], x), HalfSpace(x[1] <= 0.73, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.77, x)])
    k2 = And_([HalfSpace(0.27 <= x[1], x), HalfSpace(x[1] <= 0.48, x), HalfSpace(-0.72 <= x[2], x), HalfSpace(x[2] <= -0.52, x)])
    k3 = And_([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.23, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.77, x)])
    k4 = And_([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.77, x)])
    k5 = And_([HalfSpace(0.27 <= x[1], x), HalfSpace(x[1] <= 0.48, x), HalfSpace(-0.48 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])

    d1 = And_([HalfSpace(0.23 <= x[1], x), HalfSpace(x[1] <= 0.27, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.77, x)])
    d2 = And_([HalfSpace(0.73 <= x[1], x), HalfSpace(x[1] <= 0.77, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    d3 = And_([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.77 <= x[2], x), HalfSpace(x[2] <= -0.73, x)])
    d4 = And_([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.23, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.48, x)])
    d5 = And_([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.27 <= x[2], x), HalfSpace(x[2] <= -0.23, x)])

    w1 = And_([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.77, x), HalfSpace(-0.27 <= x[2], x), HalfSpace(x[2] <= -0.23, x)])
    w2 = And_([HalfSpace(0.48 <= x[1], x), HalfSpace(x[1] <= 0.52, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    w3 = And_([HalfSpace(0.23 <= x[1], x), HalfSpace(x[1] <= 0.48, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.48, x)])
    w4 = And_([HalfSpace(0.23 <= x[1], x), HalfSpace(x[1] <= 0.27, x), HalfSpace(-0.77 <= x[2], x), HalfSpace(x[2] <= -0.52, x)])
    w5 = And_([HalfSpace(0.73 <= x[1], x), HalfSpace(x[1] <= 0.77, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.48, x)])

    room = And_([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.00, x)])

    G(room) & G(!w1) & G(!w2) & G(!w3) & G(!w4) & G(!w5) & U(!d5,k5) & U(!d4,k4)# & U(!d3,k3) & U(!d2,k2) & U(!d1,k1)
end
x0 = [0.55, -0.75, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")
plot!([0.52, 0.73, 0.73, 0.52, 0.52], [-1.00, -1.00, -0.77, -0.77, -1.00], color=:green, linestyle=:dash, label=false)
plot!([0.27, 0.48, 0.48, 0.27, 0.27], [-0.77, -0.77, -0.52, -0.52, -0.77], color=:green, linestyle=:dash, label=false)
plot!([0.00, 0.23, 0.23, 0.00, 0.00], [-1.00, -1.00, -0.77, -0.77, -1.00], color=:green, linestyle=:dash, label=false)
plot!([0.77, 1.00, 1.00, 0.77, 0.77], [-1.00, -1.00, -0.77, -0.77, -1.00], color=:green, linestyle=:dash, label=false)
plot!([0.27, 0.48, 0.48, 0.27, 0.27], [-0.48, -0.48, -0.27, -0.27, -0.48], color=:green, linestyle=:dash, label=false)
# annotate!(0.7, -0.95, text("(1)", :green, :right))
# annotate!(0.45, -0.7, text("(2)", :green, :right))
# annotate!(0.2, -0.95, text("(3)", :green, :right))
# annotate!(0.95, -0.95, text("(4)", :green, :right))
# annotate!(0.45, -0.45, text("(5)", :green, :right))

plot!([0.23, 0.27, 0.27, 0.23, 0.23], [-1.00, -1.00, -0.77, -0.77, -1.00], color=:red, linestyle=:dash, label=false)
plot!([0.73, 0.77, 0.77, 0.73, 0.73], [-0.52, -0.52, -0.27, -0.27, -0.52], color=:red, linestyle=:dash, label=false)
plot!([0.77, 1.00, 1.00, 0.77, 0.77], [-0.77, -0.77, -0.73, -0.73, -0.77], color=:red, linestyle=:dash, label=false)
plot!([0.00, 0.23, 0.23, 0.00, 0.00], [-0.52, -0.52, -0.48, -0.48, -0.52], color=:red, linestyle=:dash, label=false)
plot!([0.77, 1.00, 1.00, 0.77, 0.77], [-0.27, -0.27, -0.23, -0.23, -0.27], color=:red, linestyle=:dash, label=false)

plot!([0.00, 0.77, 0.77, 0.00, 0.00], [-0.27, -0.27, -0.23, -0.23, -0.27], color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.48, 0.52, 0.52, 0.48, 0.48], [-0.52, -0.52, -0.27, -0.27, -0.52], color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.23, 0.48, 0.48, 0.23, 0.23], [-0.52, -0.52, -0.48, -0.48, -0.52], color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.23, 0.27, 0.27, 0.23, 0.23], [-0.77, -0.77, -0.52, -0.52, -0.77], color=:black, fill=true, fillalpha=0.3, label=false)
plot!([0.73, 0.77, 0.77, 0.73, 0.73], [-1.00, -1.00, -0.52, -0.52, -1.00], color=:black, fill=true, fillalpha=0.3, label=false)

plot!([0.0, 1.0, 1.0, 0.0, 0.0], [-1.0, -1.0, 0.0, 0.0, -1.0], color=:black, label=false)

savefig("img/doorpuzzle-1.pdf")
