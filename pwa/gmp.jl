using MomentOpt, DynamicPolynomials
using HybridSystems
using SemialgebraicSets
using HiGHS
using Ipopt
include("graph.jl")
include("qcqp.jl")

set(A::AbstractMatrix, b::Vector, x::Vector) = BasicSemialgebraicSet(FullSpace(), -A*x + b)
set((A, b)::Tuple{Matrix, Vector}, x::Vector) = set(A, b, x)
set(X::Union{HPolytope,HPolyhedron}, x::Vector) = let 
    A = stack([h.a for h in X.constraints])'
    b = stack([h.b for h in X.constraints])
    set(A, b, x)
end

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
    Ks = [(set(HybridSystems.mode(p.s, q).X, x),
           set(HybridSystems.mode(p.s, q).X, x)) for q in q0]
    Kt = [(set(HybridSystems.mode(p.s, q).X, x),
           set(HybridSystems.mode(p.s, q).X, x)) for q in qT]
    K = [Kf; Ks; Kt]

    Es = stack([nmodes(p.s)+1, q] for q in q0)'
    Et = stack([q, nmodes(p.s)+2] for q in qT)'
    E = [p.E; Es; Et]

    m = GMPModel(p.optimizer)
    set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
    set_approximation_degree(m, 2)
    b = monomials(x, 0:approximation_degree(m))
    dbdt = differentiate(b, x) * f(x,u)

    modes = collect(Set(p.E))
    @variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i][1])
    @objective m Min sum(Mom.(p.c(x,u), μ[:,2])) + 0.01*sum(Mom.(1, μ))
    cn = @constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
    @constraint m [i=modes] sum(μ[eout(E,i),1]) == sum(μ[einc(E,i),3])
    @constraint m sum(μ[eout(E,nmodes(p.s)+1),3]) == μ0
    @constraint m sum(μ[einc(E,nmodes(p.s)+2),3]) == μT
    @constraint m Mom.(1, μ[:,1]) .<= 1
    @constraint m Mom.(1, μ[:,3]) .<= 1

    optimize!(m)

    return μ, E, K, m

end

function _decode(E, p, s, t)

    v = s
    P = [s]
    for _ in 1:length(p)

        e = argmax(p .* eout(E, v))
        v = E[e, 2]
        push!(P, v)

        if v == t break end

    end

    @assert P[end] == t

    return P

end

function decode(E, C, s, t; optimizer=HiGHS.Optimizer)

    modes = setdiff(Set(E), Set([s, t]))

    m = Model(optimizer)
    @variable m 0 .<= w[1:length(C)] .<= 1
    @objective m Min -w'*C
    @constraint m [k=modes] w'*eout(E,k) == w'*einc(E,k)
    @constraint m w'*eout(E,s) == w'*einc(E,s) + 1
    @constraint m w'*eout(E,t) == w'*einc(E,t) - 1

    optimize!(m)

    return _decode(E, value.(w), s, t)

end

function extract(s, c, μ, E::AbstractMatrix, x0, xT; T=20, optimizer=Ipopt.Optimizer)

    pp = integrate.(1,μ[:,1])
    P = decode(E, log.(clamp.(pp, 1e-6, 1-1e-6)), nmodes(s)+1, nmodes(s)+2)
    P = P[2:end-1]
    qpolicy = QCQPPolicy(s, c; T=T, optimizer=optimizer)
    return action(qpolicy, (P[1], x0), (P[end], xT), P)

end

function extract(p::GMPPolicy, (q0, x0)::Tuple{Vector{Int}, Vector}, (qT, xT)::Tuple{Vector{Int}, Vector}; T=20, optimizer=Ipopt.Optimizer)

    μ, E, K, m = action(p, (q0, x0), (qT, xT))
    x, u, q, o = extract(p.s, p.c, μ, E, x0, xT, T=T, optimizer=optimizer)
    return x, u, q, o, m

end
