using Symbolics, LazySets
using Clarabel
include("pwa/product.jl")
include("pwa/gmp.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]

l = let x = Symbolics.variables(:x, 1:4)
    k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.4 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.2 <= x[2], x) & HalfSpace(x[2] <= 0.3, x)
    k2 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.1, x)
    d2 = HalfSpace(-0.2 <= x[1], x) & HalfSpace(x[1] <= -0.1, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
    w1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
    w2 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    w3 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
    lf = G(room) & U(!d1,k1)# & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)
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

x0 = [-0.4, 0.4, 0.0, 0.0]
xT = [-0.0, 0.0, 0.0, 0.0]
m, p = action(policy, (q0,x0), (qT,xT))
write("img/gmp-doors.txt","$(objective_value(m))")
write("img/gmpp-doors.txt","$(p)")
