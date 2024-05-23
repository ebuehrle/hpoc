include("interaction.jl/maps.jl")
include("pwa/product.jl")
include("pwa/gmp.jl")
include("pwa/qcqp.jl")
using MosekTools
using Ipopt

poly(l::Polyhedra.DefaultPolyhedron, x, s=nothing) = let
    s = (isnothing(s) ? (1:length(x)) : s)
    l = remove_redundant_constraints(HPolyhedron(l))
    reduce(&, [LazySets.HalfSpace(c.a'*x[s] <= c.b, x) for c in l.constraints]) 
end

lanelets = let p = planelets(scenario="DR_DEU_Roundabout_OF")
    Dict(l.id => translate(l.region, [-1000., -1000.]) for l in p)
end

x0 = [960, 1020, 0, 0] - [1000, 1000, 0, 0]
xT = [975, 1012, 0, 0] - [1000, 1000, 0, 0]

plot_all(collect(values(lanelets)), alpha=0.2)
plot!(lanelets["30025"])
scatter!(Tuple(x0[1:2]))
scatter!(Tuple(xT[1:2]))

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f = let x = Symbolics.variables(:x, 1:4)
    l1 = poly(lanelets["30025"], x, 1:2)
    G(l1)
end
c(x,u) = x[3:4]'*x[3:4] + u'*u

s, q0, qT = PPWA(A, B, f)
policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, (xq, qq), m = action(policy, (q0, x0), (qT, xT))
scatter!(xq[:,1],xq[:,2])
savefig("img/roundabout.pdf")
