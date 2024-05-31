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
l = "G(room) & (!d5 U k5) & (!d4 U k4) & (!d3 U k3) & (!d2 U k2) & (!d1 U k1)"

V = let x = Symbolics.variables(:x, 1:4)
    k1 = HPolyhedron([HalfSpace(0.55 <= x[1], x), HalfSpace(x[1] <= 0.70, x), HalfSpace(-0.95 <= x[2], x), HalfSpace(x[2] <= -0.80, x)])
    k2 = HPolyhedron([HalfSpace(0.30 <= x[1], x), HalfSpace(x[1] <= 0.45, x), HalfSpace(-0.70 <= x[2], x), HalfSpace(x[2] <= -0.55, x)])
    k3 = HPolyhedron([HalfSpace(0.05 <= x[1], x), HalfSpace(x[1] <= 0.20, x), HalfSpace(-0.95 <= x[2], x), HalfSpace(x[2] <= -0.80, x)])
    k4 = HPolyhedron([HalfSpace(0.80 <= x[1], x), HalfSpace(x[1] <= 0.95, x), HalfSpace(-0.95 <= x[2], x), HalfSpace(x[2] <= -0.80, x)])
    k5 = HPolyhedron([HalfSpace(0.30 <= x[1], x), HalfSpace(x[1] <= 0.45, x), HalfSpace(-0.45 <= x[2], x), HalfSpace(x[2] <= -0.30, x)])

    d1 = HPolyhedron([HalfSpace(0.23 <= x[1], x), HalfSpace(x[1] <= 0.27, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.77, x)])
    d2 = HPolyhedron([HalfSpace(0.73 <= x[1], x), HalfSpace(x[1] <= 0.77, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    d3 = HPolyhedron([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.77 <= x[2], x), HalfSpace(x[2] <= -0.73, x)])
    d4 = HPolyhedron([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.23, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.48, x)])
    d5 = HPolyhedron([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.27 <= x[2], x), HalfSpace(x[2] <= -0.23, x)])

    w1 = HPolyhedron([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.77, x), HalfSpace(-0.27 <= x[2], x), HalfSpace(x[2] <= -0.23, x)])
    w2 = HPolyhedron([HalfSpace(0.48 <= x[1], x), HalfSpace(x[1] <= 0.52, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    w3 = HPolyhedron([HalfSpace(0.23 <= x[1], x), HalfSpace(x[1] <= 0.48, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.48, x)])
    w4 = HPolyhedron([HalfSpace(0.23 <= x[1], x), HalfSpace(x[1] <= 0.27, x), HalfSpace(-0.77 <= x[2], x), HalfSpace(x[2] <= -0.52, x)])
    w5 = HPolyhedron([HalfSpace(0.73 <= x[1], x), HalfSpace(x[1] <= 0.77, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.48, x)])

    r1 = HPolyhedron([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.23 <= x[2], x), HalfSpace(x[2] <= -0.00, x)])
    r2 = HPolyhedron([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.48, x), HalfSpace(-0.48 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    r3 = HPolyhedron([HalfSpace(0.52 <= x[1], x), HalfSpace(x[1] <= 0.73, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    r4 = HPolyhedron([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.52 <= x[2], x), HalfSpace(x[2] <= -0.27, x)])
    r5 = HPolyhedron([HalfSpace(0.00 <= x[1], x), HalfSpace(x[1] <= 0.23, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.52, x)])
    r6 = HPolyhedron([HalfSpace(0.27 <= x[1], x), HalfSpace(x[1] <= 0.73, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.52, x)])
    r7 = HPolyhedron([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-0.73 <= x[2], x), HalfSpace(x[2] <= -0.52, x)])
    r8 = HPolyhedron([HalfSpace(0.77 <= x[1], x), HalfSpace(x[1] <= 1.00, x), HalfSpace(-1.00 <= x[2], x), HalfSpace(x[2] <= -0.77, x)])

    [k1, k2, k3, k4, k5, d1, d2, d3, d4, d5, r1, r2, r3, r4, r5, r6, r7, r8]
end
O = [
    ["k1", "room"], ["k2", "room"], ["k3", "room"], ["k4", "room"], ["k5", "room"],
    ["d1", "room"], ["d2", "room"], ["d3", "room"], ["d4", "room"], ["d5", "room"],
    ["room"], ["room"], ["room"], ["room"], ["room"], ["room"], ["room"], ["room"],
]

x0 = [0.6, -0.6, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false, V=V, O=O)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

policy = MIQPPolicy(s, c, h=0.3, T=30, optimizer=Gurobi.Optimizer)
uq, (xq, qq), mq = action(policy, (q0, x0), (qT, xT))

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2))")
plot!([0.55, 0.70, 0.70, 0.55, 0.55], [-0.95, -0.95, -0.80, -0.80, -0.95], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.30, 0.45, 0.45, 0.30, 0.30], [-0.70, -0.70, -0.55, -0.55, -0.70], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.05, 0.20, 0.20, 0.05, 0.05], [-0.95, -0.95, -0.80, -0.80, -0.95], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.80, 0.95, 0.95, 0.80, 0.80], [-0.95, -0.95, -0.80, -0.80, -0.95], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.30, 0.45, 0.45, 0.30, 0.30], [-0.45, -0.45, -0.30, -0.30, -0.45], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
# annotate!(0.7, -0.95, text("(1)", :green, :right))
# annotate!(0.45, -0.7, text("(2)", :green, :right))
# annotate!(0.2, -0.95, text("(3)", :green, :right))
# annotate!(0.95, -0.95, text("(4)", :green, :right))
# annotate!(0.45, -0.45, text("(5)", :green, :right))

plot!([0.23, 0.27, 0.27, 0.23, 0.23], [-1.00, -1.00, -0.77, -0.77, -1.00], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.73, 0.77, 0.77, 0.73, 0.73], [-0.52, -0.52, -0.27, -0.27, -0.52], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.77, 1.00, 1.00, 0.77, 0.77], [-0.77, -0.77, -0.73, -0.73, -0.77], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.00, 0.23, 0.23, 0.00, 0.00], [-0.52, -0.52, -0.48, -0.48, -0.52], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([0.77, 1.00, 1.00, 0.77, 0.77], [-0.27, -0.27, -0.23, -0.23, -0.27], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)

plot!([0.00, 0.77, 0.77, 0.00, 0.00], [-0.27, -0.27, -0.23, -0.23, -0.27], color=:black, fill=true, fillalpha=0.4, label=false)
plot!([0.48, 0.52, 0.52, 0.48, 0.48], [-0.52, -0.52, -0.27, -0.27, -0.52], color=:black, fill=true, fillalpha=0.4, label=false)
plot!([0.23, 0.48, 0.48, 0.23, 0.23], [-0.52, -0.52, -0.48, -0.48, -0.52], color=:black, fill=true, fillalpha=0.4, label=false)
plot!([0.23, 0.27, 0.27, 0.23, 0.23], [-0.77, -0.77, -0.52, -0.52, -0.77], color=:black, fill=true, fillalpha=0.4, label=false)
plot!([0.73, 0.77, 0.77, 0.73, 0.73], [-1.00, -1.00, -0.52, -0.52, -1.00], color=:black, fill=true, fillalpha=0.4, label=false)

plot!([0.0, 1.0, 1.0, 0.0, 0.0], [-1.0, -1.0, 0.0, 0.0, -1.0], color=:black, label=false)

savefig("img/doorpuzzle-1b-miqp.pdf")
