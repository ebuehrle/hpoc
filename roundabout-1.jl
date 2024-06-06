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
s1 = [
    HPolyhedron(polyhedron(convexhull([-40, 18], [-40, 22], [-20, 9], [-18, 11]))),
    HPolyhedron(polyhedron(convexhull([-20, 9], [-18, 11], [-15, 1], [-10, 2]))),
    HPolyhedron(polyhedron(convexhull([-15, 1], [-10, 2], [-9, -10], [-6, -7]))),
    HPolyhedron(polyhedron(convexhull([-9, -10], [-6, -7], [6, -7], [10, -10]))),
    HPolyhedron(polyhedron(convexhull([6, -7], [10, -10], [16, -5], [18, -2]))),
    HPolyhedron(polyhedron(convexhull([16, -5], [18, -2], [40, -9], [40, -6]))),
]
o1 = ["v1l1", "v1l2", "v1l3", "v1l4", "v1l5", "v1l6"]
s2 = [
    HPolyhedron(polyhedron(convexhull([-2, 11], [-2, 17], [-12, 11], [-9, 9]))),
    HPolyhedron(polyhedron(convexhull([-12, 11], [-9, 9], [-15, 1], [-12, 2]))),
    HPolyhedron(polyhedron(convexhull([-15, 1], [-12, 2], [-8, -9], [-5.5, -7]))),
    HPolyhedron(polyhedron(convexhull([-8, -9], [-5.5, -7], [-2.5, -18], [-1.5, -18]))),
    HPolyhedron(polyhedron(convexhull([-2.5, -18], [-1.5, -18], [5, -40], [8, -40]))),
    HPolyhedron(polyhedron(convexhull([5, -40], [8, -40], [12, -55], [15, -55]))),
]
o2 = ["v2l1", "v2l2", "v2l3", "v2l4", "v2l5", "v2l6"]

c(x,u) = x[3:4]'*x[3:4] + 2*x[7:8]'*x[7:8] + 100*u'*u + 0.1
A1 = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B1 = [0 0; 0 0; 1 0; 0 1]
A = [A1 zeros(4,4); zeros(4,4) A1]
B = [B1 zeros(4,2); zeros(4,2) B1]

V = let x = Symbolics.variables(:x, 1:8)
    p1 = [hpoly(h, x, 1:2) for h in s1]
    p2 = [hpoly(h, x, 5:6) for h in s2]
    [intersection(h1,h2) for h1 in p1 for h2 in p2]
end
O = [[h1, h2] for h1 in o1 for h2 in o2]
l = "G(v1l1 | v1l2 | v1l3 | v1l4 | v1l5 | v1l6) & G(v2l1 | v2l2 | v2l3 | v2l4 | v2l5 | v2l6) & G!(v1l2 & v2l1) & G!(v1l2 & v2l2)"

x0 = [-40, 20, 0.0, 0.0, -2, 14, 0, 0]
xT = [25, -5, 0, 0, 6, -40, 0, 0]

s, q0, qT = PPWA(A, B, l, V=V, O=O, merge_modes=false)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")

policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, tq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

x1 = linear_interpolation(tq, xq[:,1])
x2 = linear_interpolation(tq, xq[:,2])
x5 = linear_interpolation(tq, xq[:,5])
x6 = linear_interpolation(tq, xq[:,6])
tt = tq[1]:5.0:tq[end]

plot_all(lanelets, alpha=0.2, label=false, color=:grey)
#plot_all!(s1, linecolor=:red, fillcolor=nothing, fillalpha=0, label=false)
#plot_all!(s2, linecolor=:red, fillcolor=nothing, fillalpha=0, label=false)
scatter!(x5.(tt),x6.(tt),label=false,color=:orange)
scatter!(x1.(tt),x2.(tt),label=false,color=:blue)
savefig("img/roundabout-1.pdf")
