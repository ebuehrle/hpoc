using SumOfSquares, DynamicPolynomials
using MultivariateMoments
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
    K::Vector{Tuple{Matrix, Array}}
    E::Matrix
    optimizer
end

GMPPolicy(s, c; optimizer) = let 
    E = stack([(source(s,t),target(s,t)) for t in HybridSystems.transitions(s)])'
    K = [(stack([h.a for h in HybridSystems.mode(s,i).X.constraints])',
          stack([h.b for h in HybridSystems.mode(s,i).X.constraints])) for i in 1:nmodes(s)]
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

    Kf = [set(HybridSystems.mode(p.s, q).X, x) for q in 1:nmodes(p.s)]
    Ks = [set(HybridSystems.mode(p.s, q).X, x) for q in q0]
    Kt = [set(HybridSystems.mode(p.s, q).X, x) for q in qT]
    K = [Kf; Ks; Kt]

    Es = stack([nmodes(p.s)+1, q] for q in q0)'
    Et = stack([q, nmodes(p.s)+2] for q in qT)'
    E = [p.E; Es; Et]

    f(x,u) = A*x + B*u

    println("formulating SOS")
    m = SOSModel(p.optimizer)
    @variable m V[i=1:length(p.K)] Poly(monomials(x, 0:2))
    @variable m V0 Poly(monomials(x, 0:2))
    @variable m VT Poly(monomials(x, 0:2))
    c2 = @constraint m [i=1:length(p.K)] differentiate(V[i],x)'*f(x,u) >= -p.c(x,u) domain=K[i]
    c1 = @constraint m [i=1:size(p.E,1)] V[E[i,1]] <= V[E[i,2]] domain=@set K[E[i,1]] && K[E[i,2]]
    cT = @constraint m [i=1:length(qT)] V[qT[i]] <= VT
    c0 = @constraint m [i=1:length(q0)] V0 <= V[q0[i]]
    @constraint m VT(xT) <= 0
    @objective m Max V0(x0)

    println("formulating SDP")
    optimize!(m)

    μf = collect(zip(moments.(c1), [moments(c2[e[1]]) for e in eachrow(p.E)]))
    μ0 = collect(zip(moments.(c0), moments.(c0)))
    μT = collect(zip(moments.(cT), moments.(cT)))
    μ = [μf; μ0; μT]
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

    p = [m[1].a[1] for m in μ]
    H = [m[2].a[1] for m in μ]
    P = decode(E, log.(clamp.(p, 1e-6, 1-1e-6)), nmodes(s)+1, nmodes(s)+2)
    P = P[2:end-1]
    qpolicy = QCQPPolicy(s, c; T=T, optimizer=optimizer)
    return action(qpolicy, (P[1], x0), (P[end], xT), P, H)

end

function extract(p::GMPPolicy, (q0, x0)::Tuple{Vector{Int}, Vector}, (qT, xT)::Tuple{Vector{Int}, Vector}; T=20, optimizer=Ipopt.Optimizer)

    μ, E, K, m = action(p, (q0, x0), (qT, xT))
    u, x, q, t, o = extract(p.s, p.c, μ, E, x0, xT, T=T, optimizer=optimizer)
    return u, x, q, t, o, m

end
