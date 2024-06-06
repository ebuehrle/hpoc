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
lanelet_map = planelets("interaction-dataset/", scenario="DR_DEU_Merging_MT")
lanelets = [translate(l.region, [-900., -1007.]) for l in lanelet_map]
s1 = [
    HPolyhedron(polyhedron(convexhull([100, -2], [100, 1.5], [75, 1.5], [75, 3]))),
    HPolyhedron(polyhedron(convexhull([75, 1.5], [75, 3], [40, -2], [40, -3]))),
    HPolyhedron(polyhedron(convexhull([40, -2], [40, -3], [0, -1], [0, 1.5]))),
]
o1 = ["v1l1", "v1l2", "v1l3"]
s2 = [
    HPolyhedron(polyhedron(convexhull([100, -4], [100, -2], [75, -1], [75, 0]))),
    HPolyhedron(polyhedron(convexhull([75, -1], [75, 0], [40, -2], [40, -3]))),
    HPolyhedron(polyhedron(convexhull([40, -2], [40, -3], [0, -1], [0, 1.5]))),
]
o2 = ["v2l1", "v2l2", "v2l3"]

c(x,u) = x[3:4]'*x[3:4] + x[7:8]'*x[7:8] + 100*u'*u + 0.1
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
l = "G(v1l1 | v1l2 | v1l3) & G(v2l1 | v2l2 | v2l3) & G!(v1l2 & v2l2)"

x0 = [100, 0, 0.0, 0.0, 100, -3, 0, 0]
xT = [0, 0, -1, 0, 0, 0, -1, 0]

s, q0, qT = PPWA(A, B, l, V=V, O=O, merge_modes=false)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")

policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, tq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

x1 = linear_interpolation(tq, xq[:,1])
x2 = linear_interpolation(tq, xq[:,2])
x5 = linear_interpolation(tq, xq[:,5])
x6 = linear_interpolation(tq, xq[:,6])
tt = tq[1]:2:tq[end]

plot_all(lanelets, alpha=0.2, label=false, ratio=:equal, color=:grey)
#plot_all!(s1, linecolor=:red, fillcolor=nothing, fillalpha=0, label=false)
#plot_all!(s2, linecolor=:red, fillcolor=nothing, fillalpha=0, label=false)
scatter!(x5.(tt),x6.(tt),label=false,color=:blue)
scatter!(x1.(tt),x2.(tt),label=false,color=:orange)
xlims!(25,90)
savefig("img/merge-1.pdf")
