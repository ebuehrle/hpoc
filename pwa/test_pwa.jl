include("pwa.jl")

@variables x[1:4]
f = G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))

V, E, K, q0, qT = pwa(LTLTranslator(), f)

using SemialgebraicSets, DynamicPolynomials

@polyvar x[1:4]
K = [BasicSemialgebraicSet(FullSpace(), -A*x+b) for (A,b) in K]
