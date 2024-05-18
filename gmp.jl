using Symbolics, LazySets
using Clarabel
include("pwa/product.jl")
include("pwa/gmp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))
end

s, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

c(x,u) = x'*x + u'*u
policy = GMPPolicy(s, q0, qT, c, optimizer=Clarabel.Optimizer)
println(HybridSystems.nmodes(policy.s.fh), " modes")
println(HybridSystems.ntransitions(policy.s.fh), " transitions")
println(policy.s.fq0)
println(policy.s.fqT)

x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]
m, p = action(policy, (q0,x0), (qT,xT))
write("img/gmp.txt","$(objective_value(m))")
write("img/gmpp.txt","$(p)")
