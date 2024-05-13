using MomentOpt, DynamicPolynomials, Clarabel
using Symbolics, LazySets, SemialgebraicSets
include("pwa/pwa.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u

Symbolics.@variables x[1:4]
l = G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))

@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [-1.0, -1.0, 0.0, 0.5, 0.0, 0.0])
μT = DiracMeasure([x;u], [-0.0, -0.0, 0.0, 0.0, 0.0, 0.0])

V, E, K, q0, qT = pwa(LTLTranslator(), l)
K = [BasicSemialgebraicSet(FullSpace(), -A*x+b) for (A,b) in K]
Ks = [BasicSemialgebraicSet(FullSpace(), Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Float64}[]) for q in q0]
Kt = [BasicSemialgebraicSet(FullSpace(), Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Float64}[]) for q in qT]
K = [K;Ks;Kt]
Es = [(length(V)+1,q) for q in q0]
Et = [(q,length(V)+2) for q in qT]
E = [E;Es;Et]

eout(i) = [k == i for (k,_) in E]
einc(i) = [k == i for (_,k) in E]

m = GMPModel(Clarabel.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 2)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f(x,u)

@variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i])
@objective m Min sum(Mom.(c(x,u),μ[:,2]))
@constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
@constraint m sum(μ[eout(length(V)+1),1]) == μ0
@constraint m sum(μ[einc(length(V)+2),3]) == μT
@constraint m [i=1:length(V)] sum(μ[eout(i),1]) == sum(μ[einc(i),3])

optimize!(m)
write("img/gmp.txt","$(objective_value(m))")
p = integrate.(1,μ[:,1])
write("img/gmpp.txt","$(p)")
