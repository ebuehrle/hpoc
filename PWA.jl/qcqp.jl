using HybridSystems
using Graphs
using LinearAlgebra

struct QCQPPolicy
    s::HybridSystem
    c::Function
    T::Integer
    K::Vector{Tuple{Matrix, Array}}
    E::Matrix
    optimizer
end

QCQPPolicy(s::HybridSystem, c::Function; T=20, optimizer) = QCQPPolicy(
    s, c, T,
    [(stack([h.a for h in HybridSystems.mode(s,i).X.constraints])', 
      stack([h.b for h in HybridSystems.mode(s,i).X.constraints])) for i in 1:nmodes(s)],
    adjacency_matrix(s.automaton.G) + I,
    optimizer
)

function action(p::QCQPPolicy, (q0, x0), (qT, xT), P, H=nothing; S=1e3)

    @assert P[1] == q0
    @assert P[end] == qT

    if isnothing(H) H = ones(size(P)) end

    A = HybridSystems.mode(p.s, 1).A
    B = HybridSystems.mode(p.s, 1).B
    nx, nu = size(B)
    M = length(P)

    m = Model(p.optimizer)
    @variable m x[1:M, 1:p.T, 1:nx]
    @variable m u[1:M, 1:p.T, 1:nu]
    @variable m h[i=1:M] .>= 1e-4 start=H[i]/p.T
    @variable m 0 .<= s[1:M, 1:p.T] .<= 0.01
    
    @objective m Min sum(p.c(x[k,t,:], u[k,t,:]) * h[k] for t=1:p.T for k=1:M) + S*sum(s.^2)

    @constraint m [k=1:M] x[k,2:end,:]' .== x[k,1:end-1,:]' + h[k]*A*x[k,1:end-1,:]' + h[k]*B*u[k,1:end-1,:]'
    @constraint m [k=1:M,t=1:p.T] p.K[P[k]][1] / norm(p.K[P[k]][1]) * x[k,t,:] .<= p.K[P[k]][2] / norm(p.K[P[k]][1]) .+ s[k,t]
    @constraint m [k=2:M] x[k,1,:] .== x[k-1,end,:] + h[k-1]*A*x[k-1,end,:] + h[k-1]*B*u[k-1,end,:]
    @constraint m x[1,1,:] .== x0
    @constraint m x[end,end,:] .== xT

    optimize!(m)

    xr = reshape(permutedims(value.(x), [2, 1, 3]), (M*p.T,nx))
    ur = reshape(permutedims(value.(u), [2, 1, 3]), (M*p.T,nu))
    qr = repeat(P, inner=p.T)
    tr = [0; cumsum(repeat(value.(h), inner=p.T))[1:end-1]]

    return ur, xr, qr, tr, m

end
