include("PWA.jl/PWA.jl")
using .PWA
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots; ENV["GKSwstype"] = "100"
using Interpolations

hpoly(h, x, s) = HPolyhedron([LazySets.HalfSpace(c.a' * x[s] <= c.b, x) for c in h.constraints])

include("interaction.jl/maps.jl")
lanelet_map = planelets("interaction-dataset/", scenario="DR_DEU_Roundabout_OF")
lanelets = [translate(l.region, [-1000., -1000.]) for l in lanelet_map]
segments = [
    HPolyhedron(polyhedron(convexhull([-40, 18], [-40, 22], [-20, 8], [-18, 11]))),
    HPolyhedron(polyhedron(convexhull([-20, 8], [-18, 11], [-16, 1], [-10, 2]))),
    HPolyhedron(polyhedron(convexhull([-16, 1], [-10, 2], [-9, -10], [-6, -5]))),
    HPolyhedron(polyhedron(convexhull([-9, -10], [-6, -5], [5, -5], [10, -10]))),
    HPolyhedron(polyhedron(convexhull([5, -5], [10, -10], [16, -6], [18, -2]))),
    HPolyhedron(polyhedron(convexhull([16, -6], [18, -2], [40, -9], [40, -6]))),
]

c(x,u) = x[3:4]'*x[3:4] + 100*u'*u + 0.1
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
V = let x = Symbolics.variables(:x, 1:4)
    p1 = [hpoly(h, x, 1:2) for h in segments]
    p1
end
O = [
    ["v1l1"], ["v1l2"], ["v1l3"], ["v1l4"], ["v1l5"], ["v1l6"],
]
l = "G(v1l1 | v1l2 | v1l3 | v1l4 | v1l5 | v1l6)"
x0 = [975, 1012, 0.0, 0.0] + [-1000, -1000, 0, 0]
xT = [1025, 995, 0.0, 0.0] + [-1000, -1000, 0, 0]

s, q0, qT = PPWA(A, B, l, V=V, O=O, merge_modes=false)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, tq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

x1 = linear_interpolation(tq, xq[:,1])
x2 = linear_interpolation(tq, xq[:,2])
tt = tq[1]:5.0:tq[end]

plot_all(lanelets, alpha=0.2, label=false)
plot_all!(segments, linecolor=:red, fillcolor=nothing, fillalpha=0, label=false)
scatter!(x1.(tt),x2.(tt),label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))",legend=:topright)
savefig("img/roundabout.pdf")
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
