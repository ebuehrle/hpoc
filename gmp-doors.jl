using MomentOpt, DynamicPolynomials, Clarabel
using Symbolics, LazySets, SemialgebraicSets
include("pwa/pwa.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u
x0 = [-0.4, 0.4, 0.0, 0.0]
xT = [-0.0, 0.0, 0.0, 0.0]

Symbolics.@variables x[1:4]
k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.4 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.2 <= x[2], x) & HalfSpace(x[2] <= 0.3, x)
k2 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.1, x)
d2 = HalfSpace(-0.2 <= x[1], x) & HalfSpace(x[1] <= -0.1, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
w1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
w2 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
w3 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
lf = G(room) & U(!d1,k1)# & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)

@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [x0; 0.0; 0.0])
μT = DiracMeasure([x;u], [xT; 0.0; 0.0])

hs, q0, qT, K0, KT = fpwa(A, B, lf, x0, xT, LTLTranslator())
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
@objective m Min sum(Mom.(c(x,u),μ[:,2])) + 0.01*sum(Mom.(1, μ))
@constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
@constraint m [i=1:nmodes(hs)] sum(μ[eout(i),1]) == sum(μ[einc(i),3])
@constraint m sum(μ[eout(nmodes(hs)+1),1]) == μ0
@constraint m sum(μ[einc(nmodes(hs)+2),3]) == μT
@constraint m Mom.(1, μ[:,1]) .<= 1
@constraint m Mom.(1, μ[:,3]) .<= 1

optimize!(m)
write("img/gmp-doors.txt","$(objective_value(m))")
p = integrate.(1,μ[:,1])
write("img/gmpp-doors.txt","$(p)")
