using MomentOpt, DynamicPolynomials
using HybridSystems
using SemialgebraicSets
include("fpwa.jl")

set(X::Union{HPolytope,HPolyhedron}, x::Vector) = BasicSemialgebraicSet(FullSpace(), 
    - stack([h.a for h in X.constraints])'*x 
    + stack([h.b for h in X.constraints]))

eout(E,i) = E[:,1] .== i
einc(E,i) = E[:,2] .== i

struct GMPPolicy
    ps::HybridSystem
    pq0::Vector{Int}
    pqT::Vector{Int}
    s::FPWA
    c::Function
    optimizer
end

GMPPolicy(s, q0, qT, c; optimizer) = let f = FPWA(s, q0, qT)
    GMPPolicy(s, q0, qT, f, c, optimizer)
end

function action(p::GMPPolicy, (q0,x0), (qT,xT))

    A = HybridSystems.mode(p.ps,1).A
    B = HybridSystems.mode(p.ps,1).B
    nx, nu = size(B)

    @polyvar x[1:nx] u[1:nu]
    μ0 = DiracMeasure([x;u], [x0; zeros(nu)])
    μT = DiracMeasure([x;u], [xT; zeros(nu)])

    f(x,u) = A*x + B*u

    Q0 = [q for q in q0 if x0 ∈ HybridSystems.mode(p.ps,q).X]
    q0 = [i for (i,(j,k)) in enumerate(p.s.V) if (j in Q0) || (k in Q0)]
    K0 = [HybridSystems.mode(p.ps, first(Set((j,k)) ∩ Set(Q0))).X for (j,k) in p.s.V if (j in Q0) || (k in Q0)]
    println(q0)

    QT = [q for q in qT if xT ∈ HybridSystems.mode(p.ps,q).X]
    qT = [i for (i,(j,k)) in enumerate(p.s.V) if (j in QT) || (k in QT)]
    KT = [HybridSystems.mode(p.ps, first(Set((j,k)) ∩ Set(QT))).X for (j,k) in p.s.V if (j in QT) || (k in QT)]
    println(qT)

    Kf = [(set(HybridSystems.mode(p.s.fh,source(p.s.fh,t)).X, x),
           set(HybridSystems.resetmap(p.s.fh,t).X, x), 
           set(HybridSystems.mode(p.s.fh,target(p.s.fh,t)).X, x)) for t in HybridSystems.transitions(p.s.fh)]
    Kt = [(set(HybridSystems.mode(p.s.fh,q).X, x),
           set(k, x),
           set(k, x)) for (q,k) in zip(qT,KT)]
    Ks = [(set(k, x),
           set(k, x),
           set(HybridSystems.mode(p.s.fh,q).X, x)) for (q,k) in zip(q0,K0)]
    K = [Kf; Kt; Ks]

    Ef = stack([source(p.s.fh, t), target(p.s.fh, t)] for t in HybridSystems.transitions(p.s.fh))'
    Et = stack([q, nmodes(p.s.fh)+2] for q in qT)'
    Es = stack([nmodes(p.s.fh)+1, q] for q in q0)'
    E = [Ef; Et; Es]

    m = GMPModel(p.optimizer)
    set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
    set_approximation_degree(m, 2)
    b = monomials(x, 0:approximation_degree(m))
    dbdt = differentiate(b, x) * f(x,u)

    @variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i][j])
    @objective m Min sum(Mom.(p.c(x,u),μ[:,2])) + 0.01*sum(Mom.(1, μ))
    @constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
    @constraint m [i=1:nmodes(p.s.fh)] sum(μ[eout(E,i),1]) == sum(μ[einc(E,i),3])
    @constraint m sum(μ[eout(E,nmodes(p.s.fh)+1),1]) == μ0
    @constraint m sum(μ[einc(E,nmodes(p.s.fh)+2),3]) == μT
    @constraint m Mom.(1, μ[:,1]) .<= 1
    @constraint m Mom.(1, μ[:,3]) .<= 1

    optimize!(m)

    p = integrate.(1,μ[:,1])

    return m, p

end
