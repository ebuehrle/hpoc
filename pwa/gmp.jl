using MomentOpt, DynamicPolynomials
using HybridSystems
using SemialgebraicSets
include("fpwa.jl")

set(A::AbstractMatrix, b::Vector, x::Vector) = BasicSemialgebraicSet(FullSpace(), -A*x + b)
set((A, b)::Tuple{Matrix, Vector}, x::Vector) = set(A, b, x)
set(X::Union{HPolytope,HPolyhedron}, x::Vector) = let 
    A = stack([h.a for h in X.constraints])'
    b = stack([h.b for h in X.constraints])
    set(A, b, x)
end

eout(E,i) = E[:,1] .== i
einc(E,i) = E[:,2] .== i

struct GMPPolicy
    s::HybridSystem
    c::Function
    K::Vector{Tuple{Tuple{Matrix, Array}, Tuple{Matrix, Array}}}
    E::Matrix
    optimizer
end

GMPPolicy(s, c; optimizer) = let 
    E = stack([(source(s,t),target(s,t)) for t in HybridSystems.transitions(s)])'
    Ks= [(stack([h.a for h in HybridSystems.mode(s,i).X.constraints])',
          stack([h.b for h in HybridSystems.mode(s,i).X.constraints])) for (i,_) in eachrow(E)]
    Kt= [(stack([h.a for h in HybridSystems.mode(s,i).X.constraints])',
          stack([h.b for h in HybridSystems.mode(s,i).X.constraints])) for (_,i) in eachrow(E)]
    K = collect(zip(Ks, Kt))
    GMPPolicy(s, c, K, E, optimizer)
end

function action(p::GMPPolicy, (q0,x0), (qT,xT))

    if !(q0 isa Array) q0 = [q0] end
    if !(qT isa Array) qT = [qT] end

    q0 = [q for q in q0 if x0 ∈ HybridSystems.mode(p.s,q).X]
    qT = [q for q in qT if xT ∈ HybridSystems.mode(p.s,q).X]

    A = HybridSystems.mode(p.s,1).A
    B = HybridSystems.mode(p.s,1).B
    nx, nu = size(B)

    @polyvar x[1:nx] u[1:nu]
    μ0 = DiracMeasure([x;u], [x0; zeros(nu)])
    μT = DiracMeasure([x;u], [xT; zeros(nu)])

    f(x,u) = A*x + B*u

    Kf = [(set(k1, x), set(k2, x)) for (k1, k2) in p.K]
    Kt = [(set(HybridSystems.mode(p.s, q).X, x),
           set(HybridSystems.mode(p.s, q).X, x)) for q in qT]
    Ks = [(set(HybridSystems.mode(p.s, q).X, x),
           set(HybridSystems.mode(p.s, q).X, x)) for q in q0]
    K = [Kf; Kt; Ks]

    Et = stack([q, nmodes(p.s)+2] for q in qT)'
    Es = stack([nmodes(p.s)+1, q] for q in q0)'
    E = [p.E; Et; Es]

    m = GMPModel(p.optimizer)
    set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
    set_approximation_degree(m, 2)
    b = monomials(x, 0:approximation_degree(m))
    dbdt = differentiate(b, x) * f(x,u)

    supps = [(k1,k1,k2,k2) for (k1,k2) in K]
    modes = collect(Set(p.E))
    @variable m μ[i=1:length(K),j=1:4] Meas([x;u], support=supps[i][j])
    @objective m Min sum(Mom.(p.c(x,u), μ[:,2] + μ[:,3])) + 0.01*sum(Mom.(1, μ))
    @constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2] + μ[i,3]) .== Mom.(b, μ[i,4]) - Mom.(b, μ[i,1])
    @constraint m [i=modes] sum(μ[eout(E,i),1]) == sum(μ[einc(E,i),4])
    @constraint m sum(μ[eout(E,nmodes(p.s)+1),1]) == μ0
    @constraint m sum(μ[einc(E,nmodes(p.s)+2),4]) == μT
    @constraint m Mom.(1, μ[:,1]) .<= 1
    @constraint m Mom.(1, μ[:,4]) .<= 1

    optimize!(m)

    p = integrate.(1,μ[:,1])

    return m, p

end
