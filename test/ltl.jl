using MomentOpt, DynamicPolynomials, Clarabel
using Symbolics, LazySets, SemialgebraicSets
include("../pwa/pwa.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u
x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]

Symbolics.@variables x[1:4]
l = G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))

@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [x0; 0.0; 0.0])
μT = DiracMeasure([x;u], [xT; 0.0; 0.0])

V, E, K, q0, qT = ltla(LTLTranslator(), l)
hs, q0, qT = pwa(A, B, l, x0, xT, LTLTranslator())
