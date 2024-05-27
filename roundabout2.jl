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
    
    polyhedron(convexhull([0, 11], [0, 17], [-6, 10], [-7, 16])),
    polyhedron(convexhull([-6, 10], [-7, 16], [-14, 9], [-8, 6])),
    polyhedron(convexhull([-14, 9], [-8, 6], [-16, 1], [-10, 2])),
    polyhedron(convexhull([-9, -10], [-6, -5], [-4, -18], [0, -18])),
    polyhedron(convexhull([-4, -18], [0, -18], [4, -40], [9, -40])),
]
poly(h, x, s) = reduce(&, [LazySets.HalfSpace(c.a' * x[s] <= c.b, x) for c in h.constraints])

c(x,u) = x[3:4]'*x[3:4] + 100*u[1:2]'*u[1:2] + 0.1 + x[7:8]'*x[7:8] + 100*u[3:4]'*u[3:4] + 0.1
A1 = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B1 = [0 0; 0 0; 1 0; 0 1]
A = [A1 zeros(4,4); zeros(4,4) A1]
B = [B1 zeros(4,2); zeros(4,2) B1]
l = let x = Symbolics.variables(:x, 1:8)
    poly1 = [poly(HPolyhedron(s), x, 1:2) for s in segments]
    road1 = reduce(|, poly1)
    poly2 = [poly(HPolyhedron(s), x, 7:8) for s in segments]
    road2 = reduce(|, poly2)
    G(road1) & G(road2)
end
x01 = [975, 1012, 0.0, 0.0] + [-1000, -1000, 0, 0]
x02 = [1000, 1015, 0, 0] + [-1000, -1000, 0, 0]
xT1 = [1025, 995, 0.0, 0.0] + [-1000, -1000, 0, 0]
xT2 = [1006, 960, 0, 0] + [-1000, -1000, 0, 0]
x0 = [x01; x02]
xT = [xT1; xT2]

s, q0, qT = PPWA(A, B, l)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

plot_all(lanelets, alpha=0.2, label=false)
plot_all!(segments, linecolor=:red, fillcolor=false, label=false)
scatter!(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))",legend=:topright)
scatter!(xq[:,5],xq[:,6],label=false)
savefig("img/roundabout2.pdf")
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
