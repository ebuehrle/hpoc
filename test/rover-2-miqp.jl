include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using Gurobi
using Plots
using .PWA
using Spot

c(x,u) = x[3:4]'*x[3:4] + x[7:8]'*x[7:8] + u'*u + 0.1
A1 = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B1 = [0 0; 0 0; 1 0; 0 1]
A = [A1 zeros(4,4); zeros(4,4) A1]
B = [B1 zeros(4,2); zeros(4,2) B1]
#l = "G(ra) & G((r1o1 | r1o2 | r1o3) -> F(r1t1)) & G((r2o1 | r2o2 | r2o3) -> F(r2t1)) & F(r1o1 | r2o1) & F(r1o2 | r2o2) & F(r1o3 | r2o3)"# & F(r1o4 | r2o4) & F(r1t1 | r2t1) & F(r1t2 | r2t2)
l = "G(ra) & G((r1o1 | r1o2 | r1o3) -> F(r1t1)) & G((r2o1 | r2o2 | r2o3) -> F(r2t1)) & F(r1o1 | r2o1) & F(r1o2 | r2o2) & F(r1o3 | r2o3) & (!(r1o1 | r2o1 | r1o3 | r2o3) U (r1o2 | r2o2))"

V = let x = Symbolics.variables(:x, 1:8)
    r1o1 = HPolyhedron([HalfSpace(-0.4 <= x[1], x), HalfSpace(x[1] <= -0.3, x), HalfSpace(-0.4 <= x[2], x), HalfSpace(x[2] <= -0.3, x)])
    r1o2 = HPolyhedron([HalfSpace( 0.3 <= x[1], x), HalfSpace(x[1] <=  0.4, x), HalfSpace(-0.4 <= x[2], x), HalfSpace(x[2] <= -0.3, x)])
    r1o3 = HPolyhedron([HalfSpace( 0.3 <= x[1], x), HalfSpace(x[1] <=  0.4, x), HalfSpace( 0.3 <= x[2], x), HalfSpace(x[2] <=  0.4, x)])
    # r1o4 = HPolyhedron([HalfSpace(-0.4 <= x[1], x), HalfSpace(x[1] <= -0.3, x), HalfSpace( 0.3 <= x[2], x), HalfSpace(x[2] <=  0.4, x)])
    r2o1 = HPolyhedron([HalfSpace(-0.4 <= x[5], x), HalfSpace(x[5] <= -0.3, x), HalfSpace(-0.4 <= x[6], x), HalfSpace(x[6] <= -0.3, x)])
    r2o2 = HPolyhedron([HalfSpace( 0.3 <= x[5], x), HalfSpace(x[5] <=  0.4, x), HalfSpace(-0.4 <= x[6], x), HalfSpace(x[6] <= -0.3, x)])
    r2o3 = HPolyhedron([HalfSpace( 0.3 <= x[5], x), HalfSpace(x[5] <=  0.4, x), HalfSpace( 0.3 <= x[6], x), HalfSpace(x[6] <=  0.4, x)])
    # r2o4 = HPolyhedron([HalfSpace(-0.4 <= x[5], x), HalfSpace(x[5] <= -0.3, x), HalfSpace( 0.3 <= x[6], x), HalfSpace(x[6] <=  0.4, x)])

    r1t1 = HPolyhedron([HalfSpace(-0.4 <= x[1], x), HalfSpace(x[1] <= -0.3, x), HalfSpace(-0.05 <= x[2], x), HalfSpace(x[2] <= 0.05, x)])
    #r1t2 = HPolyhedron([HalfSpace( 0.3 <= x[1], x), HalfSpace(x[1] <=  0.4, x), HalfSpace(-0.05 <= x[2], x), HalfSpace(x[2] <= 0.05, x)])
    r2t1 = HPolyhedron([HalfSpace(-0.4 <= x[5], x), HalfSpace(x[5] <= -0.3, x), HalfSpace(-0.05 <= x[6], x), HalfSpace(x[6] <= 0.05, x)])
    #r2t2 = HPolyhedron([HalfSpace( 0.3 <= x[5], x), HalfSpace(x[5] <=  0.4, x), HalfSpace(-0.05 <= x[6], x), HalfSpace(x[6] <= 0.05, x)])

    # r1c = HPolyhedron([HalfSpace(-0.2 <= x[1], x), HalfSpace(x[1] <= 0.2, x), HalfSpace(-0.2 <= x[2], x), HalfSpace(x[2] <= 0.2, x)])
    # r2c = HPolyhedron([HalfSpace(-0.2 <= x[5], x), HalfSpace(x[5] <= 0.2, x), HalfSpace(-0.2 <= x[6], x), HalfSpace(x[6] <= 0.2, x)])
    
    ra = HPolyhedron([
        HalfSpace(-1.0 <= x[1], x), HalfSpace(x[1] <= 1.0, x), HalfSpace(-1.0 <= x[2], x), HalfSpace(x[2] <= 1.0, x),
        HalfSpace(-1.0 <= x[5], x), HalfSpace(x[5] <= 1.0, x), HalfSpace(-1.0 <= x[6], x), HalfSpace(x[6] <= 1.0, x),
    ])

    [
        r1o1, r1o2, r1o3,
        r2o1, r2o2, r2o3,
        r1t1, #r1t2,
        r2t1, #r2t2,
        ra,
    ]

end
O = [
    ["r1o1", "ra"], ["r1o2", "ra"], ["r1o3", "ra"],
    ["r2o1", "ra"], ["r2o2", "ra"], ["r2o3", "ra"],
    ["r1t1", "ra"], #["r1t2", "ra"],
    ["r2t1", "ra"], #["r2t2", "ra"],
    ["ra"],
]
x0 = [-0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, LTLTranslator(deterministic=true, state_based_acceptance=true, buchi=true), V=V, O=O, merge_modes=false)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

policy = MIQPPolicy(s, c, h=0.3, T=30, optimizer=Gurobi.Optimizer)
uq, (xq, qq), mq = action(policy, (q0, x0), (qT, xT))

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2))")
scatter!(xq[:,5],xq[:,6],label=false)

plot!([-0.4, -0.3, -0.3, -0.4, -0.4], [-0.4, -0.4, -0.3, -0.3, -0.4], color=:green, linestyle=:dash, fill=true, fillalpha=0.2, label=false)
plot!([ 0.3,  0.4,  0.4,  0.3,  0.3], [-0.4, -0.4, -0.3, -0.3, -0.4], color=:green, linestyle=:dash, fill=true, fillalpha=0.2, label=false)
plot!([ 0.3,  0.4,  0.4,  0.3,  0.3], [ 0.3,  0.3,  0.4,  0.4,  0.3], color=:green, linestyle=:dash, fill=true, fillalpha=0.2, label=false)
#plot!([-0.4, -0.3, -0.3, -0.4, -0.4], [ 0.3,  0.3,  0.4,  0.4,  0.3], color=:green, linestyle=:dash, fill=true, fillalpha=0.2, label=false)

plot!([-0.4, -0.3, -0.3, -0.4, -0.4], [-0.05, -0.05, 0.05, 0.05, -0.05], color=:yellow, linestyle=:dash, fill=true, fillalpha=0.2, label=false)
#plot!([ 0.3,  0.4,  0.4,  0.3,  0.3], [-0.05, -0.05, 0.05, 0.05, -0.05], color=:yellow, linestyle=:dash, fill=true, fillalpha=0.2, label=false)

#plot!([-0.2, 0.2, 0.2, -0.2, -0.2], [-0.2, -0.2, 0.2, 0.2, -0.2], color=:blue, linestyle=:dash, fill=true, fillalpha=0.2, label=false)

savefig("img/rover-2-miqp.pdf")
