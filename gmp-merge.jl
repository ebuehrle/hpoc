using Symbolics, LazySets
using Clarabel
include("pwa/product.jl")
include("pwa/gmp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    lane1 = HalfSpace(x[2] >= -3.0, x) & HalfSpace(x[2] <= -1.0, x) & HalfSpace(x[1] <= -5.0, x)
    lane2 = HalfSpace(x[2] >= -1.0, x) & HalfSpace(x[2] <=  1.0, x)
    l = G(lane1 | lane2)
end

s, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

c(x,u) = x'*x + u'*u
policy = GMPPolicy(s, c, optimizer=Clarabel.Optimizer)

x0 = [-10.0, -2.0, 10.0, 0.0]
xT = [-0.0, -0.0, 0.0, 0.0]
m, p = action(policy, (q0,x0), (qT,xT))
write("img/gmp-merge.txt","$(objective_value(m))")
write("img/gmpp-merge.txt","$(p)")
