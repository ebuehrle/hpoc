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
l = "G(room) & (!d6 U k6) & (!d5 U k5) & (!d4 U k4) & (!d3 U k3) & (!d2 U k2) & (!d1 U k1)"

V = let x = Symbolics.variables(:x, 1:4)
    k1 = HPolyhedron([HalfSpace(0.21-1 <= x[1], x), HalfSpace(x[1] <= 0.29-1, x), HalfSpace( 0.40 <= x[2], x), HalfSpace(x[2] <=  0.50, x)])
    k2 = HPolyhedron([HalfSpace(0.11-1 <= x[1], x), HalfSpace(x[1] <= 0.19-1, x), HalfSpace( 0.40 <= x[2], x), HalfSpace(x[2] <=  0.50, x)])
    k3 = HPolyhedron([HalfSpace(0.00-1 <= x[1], x), HalfSpace(x[1] <= 0.10-1, x), HalfSpace( 0.21 <= x[2], x), HalfSpace(x[2] <=  0.29, x)])
    k4 = HPolyhedron([HalfSpace(0.00-1 <= x[1], x), HalfSpace(x[1] <= 0.10-1, x), HalfSpace(-0.29 <= x[2], x), HalfSpace(x[2] <= -0.21, x)])
    k5 = HPolyhedron([HalfSpace(0.11-1 <= x[1], x), HalfSpace(x[1] <= 0.19-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= -0.40, x)])
    k6 = HPolyhedron([HalfSpace(0.21-1 <= x[1], x), HalfSpace(x[1] <= 0.29-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= -0.40, x)])

    d1 = HPolyhedron([HalfSpace(0.88-1 <= x[1], x), HalfSpace(x[1] <= 0.92-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d2 = HPolyhedron([HalfSpace(0.78-1 <= x[1], x), HalfSpace(x[1] <= 0.82-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d3 = HPolyhedron([HalfSpace(0.68-1 <= x[1], x), HalfSpace(x[1] <= 0.72-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d4 = HPolyhedron([HalfSpace(0.58-1 <= x[1], x), HalfSpace(x[1] <= 0.62-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d5 = HPolyhedron([HalfSpace(0.48-1 <= x[1], x), HalfSpace(x[1] <= 0.52-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    d6 = HPolyhedron([HalfSpace(0.38-1 <= x[1], x), HalfSpace(x[1] <= 0.42-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])

    r1 = HPolyhedron([HalfSpace(0.00-1 <= x[1], x), HalfSpace(x[1] <= 0.38-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    r2 = HPolyhedron([HalfSpace(0.42-1 <= x[1], x), HalfSpace(x[1] <= 0.48-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    r3 = HPolyhedron([HalfSpace(0.52-1 <= x[1], x), HalfSpace(x[1] <= 0.58-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    r4 = HPolyhedron([HalfSpace(0.62-1 <= x[1], x), HalfSpace(x[1] <= 0.68-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    r5 = HPolyhedron([HalfSpace(0.72-1 <= x[1], x), HalfSpace(x[1] <= 0.78-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    r6 = HPolyhedron([HalfSpace(0.82-1 <= x[1], x), HalfSpace(x[1] <= 0.88-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])
    r7 = HPolyhedron([HalfSpace(0.92-1 <= x[1], x), HalfSpace(x[1] <= 1.00-1, x), HalfSpace(-0.50 <= x[2], x), HalfSpace(x[2] <= 0.50, x)])

    [k1, k2, k3, k4, k5, k6, d1, d2, d3, d4, d5, d6, r1, r2, r3, r4, r5, r6, r7]
end
O = [
    ["k1", "room"], ["k2", "room"], ["k3", "room"], ["k4", "room"], ["k5", "room"], ["k6", "room"],
    ["d1", "room"], ["d2", "room"], ["d3", "room"], ["d4", "room"], ["d5", "room"], ["d6", "room"],
    ["room"], ["room"], ["room"], ["room"], ["room"], ["room"], ["room"],
]
x0 = [-0.8, 0.0, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false, V=V, O=O)

println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")

plot!([-0.79, -0.71, -0.71, -0.79, -0.79], [ 0.40,   0.40,  0.50,  0.50,  0.40], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.89, -0.81, -0.81, -0.89, -0.89], [ 0.40,   0.40,  0.50,  0.50,  0.40], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-1.00, -0.90, -0.90, -1.00, -1.00], [ 0.21,   0.21,  0.29,  0.29,  0.21], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-1.00, -0.90, -0.90, -1.00, -1.00], [-0.29,  -0.29, -0.21, -0.21, -0.29], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.89, -0.81, -0.81, -0.89, -0.89], [-0.50,  -0.50, -0.40, -0.40, -0.50], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.79, -0.71, -0.71, -0.79, -0.79], [-0.50,  -0.50, -0.40, -0.40, -0.50], color=:green, fill=true, fillalpha=0.2, linestyle=:dash, label=false)

plot!([-0.12, -0.08, -0.08, -0.12, -0.12], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.22, -0.18, -0.18, -0.22, -0.22], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.32, -0.28, -0.28, -0.32, -0.32], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.42, -0.38, -0.38, -0.42, -0.42], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.52, -0.48, -0.48, -0.52, -0.52], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)
plot!([-0.62, -0.58, -0.58, -0.62, -0.62], [-0.50, -0.50, 0.50, 0.50, -0.50], color=:red, fill=true, fillalpha=0.2, linestyle=:dash, label=false)

plot!([-1.0, 0.0, 0.0, -1.0, -1.0], [-0.5, -0.5, 0.5, 0.5, -0.5], color=:black, label=false)

savefig("img/doorpuzzle-2.pdf")
