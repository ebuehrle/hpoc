using MomentOpt, DynamicPolynomials, Clarabel
using Symbolics, LazySets, SemialgebraicSets
include("pwa/pwa.jl")

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

hs, q0, qT, K0, KT = fpwa(A, B, l, x0, xT, LTLTranslator())
set(X::Union{HPolytope,HPolyhedron}, x::Vector) = BasicSemialgebraicSet(FullSpace(), 
            - stack([h.a for h in X.constraints])'*x 
            + stack([h.b for h in X.constraints]))
Kf = [(set(HybridSystems.mode(hs,source(hs,t)).X, x),
        set(HybridSystems.resetmap(hs,t).X, x), 
        set(HybridSystems.mode(hs,target(hs,t)).X, x)) for t in HybridSystems.transitions(hs)]
Kt = [(set(HybridSystems.mode(hs,q).X, x),
        set(K, x),
        set(K, x)) for (q,K) in zip(qT,KT)]
Ks = [(set(K, x),
        set(K, x),
        set(HybridSystems.mode(hs,q).X, x)) for (q,K) in zip(q0,K0)]
K = [Kf; Kt; Ks]

Ef = stack([source(hs, t), target(hs, t)] for t in HybridSystems.transitions(hs))'
Et = stack([q, nmodes(hs)+2] for q in qT)'
Es = stack([nmodes(hs)+1, q] for q in q0)'
E = [Ef; Et; Es]

eout(i) = E[:,1] .== i
einc(i) = E[:,2] .== i

m = GMPModel(Clarabel.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 2)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f(x,u)

@variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i][j])
@objective m Min sum(Mom.(c(x,u),μ[:,2]))
@constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
@constraint m [i=1:nmodes(hs)] sum(μ[eout(i),1]) == sum(μ[einc(i),3])
@constraint m sum(μ[eout(nmodes(hs)+1),1]) == μ0
@constraint m sum(μ[einc(nmodes(hs)+2),3]) == μT

optimize!(m)
write("img/gmp.txt","$(objective_value(m))")
p = integrate.(1,μ[:,1])
write("img/gmpp.txt","$(p)")
