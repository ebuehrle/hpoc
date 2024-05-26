include("PWA.jl/PWA.jl")
using .PWA
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots

include("interaction.jl/maps.jl")
smap = planelets("interaction-dataset/", scenario="DR_DEU_Roundabout_OF")
lanelets = [translate(l.region, [-1000., -1000.]) for l in smap]
segments = [
    polyhedron(convexhull([-40, 18], [-40, 22], [-20, 8], [-18, 11])),
    polyhedron(convexhull([-20, 8], [-18, 11], [-16, 1], [-10, 2])),
    polyhedron(convexhull([-16, 1], [-10, 2], [-9, -10], [-6, -5])),
    polyhedron(convexhull([-9, -10], [-6, -5], [5, -5], [10, -10])),
    polyhedron(convexhull([5, -5], [10, -10], [16, -6], [18, -2])),
    polyhedron(convexhull([16, -6], [18, -2], [40, -9], [40, -6])),
]
poly(h, x, s) = reduce(&, [LazySets.HalfSpace(c.a' * x[s] <= c.b, x) for c in h.constraints])

c(x,u) = x[3:4]'*x[3:4] + 100*u'*u + 0.1
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    polys = [poly(HPolyhedron(s), x, 1:2) for s in segments]
    road = reduce(|, polys)
    G(road)
end
x0 = [975, 1012, 0.0, 0.0] + [-1000, -1000, 0, 0]
xT = [1025, 995, 0.0, 0.0] + [-1000, -1000, 0, 0]

s, q0, qT = PPWA(A, B, l)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

plot_all(lanelets, alpha=0.2, label=false)
plot_all!(segments, linecolor=:red, fillcolor=false, label=false)
scatter!(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))",legend=:topright)
savefig("img/roundabout.pdf")
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
