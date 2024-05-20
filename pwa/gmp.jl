using MomentOpt, DynamicPolynomials
using HybridSystems
using SemialgebraicSets

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
    @objective m Min sum(Mom.(p.c(x,u), μ[:,2]))# + 0.01*sum(Mom.(1, μ))
    cn = @constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
    @constraint m [i=modes] sum(μ[eout(E,i),1]) == sum(μ[einc(E,i),3])
    @constraint m sum(μ[eout(E,nmodes(p.s)+1),3]) == μ0
    @constraint m sum(μ[einc(E,nmodes(p.s)+2),3]) == μT
    #@constraint m Mom.(1, μ[:,1]) .<= 1
    #@constraint m Mom.(1, μ[:,3]) .<= 1

    optimize!(m)

    pp = integrate.(1,μ[:,1])

    q0 = Es[argmax(pp[eout(E,nmodes(p.s)+1)]),2]
    p0 = pp .* eout(E,q0)
    c0 = cn[argmax(p0)]
    v0 = first.(-dual.(c0))' * b
    dv0 = differentiate(v0, x)
    dv0 = [d(x0) for d in dv0]

    return dv0, m, pp

end
