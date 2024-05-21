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

QCQPPolicy(s::HybridSystem, c::Function; T=10, optimizer) = QCQPPolicy(
    s, c, T,
    [(stack([h.a for h in HybridSystems.mode(s,i).X.constraints])', 
      stack([h.b for h in HybridSystems.mode(s,i).X.constraints])) for i in 1:nmodes(s)],
    adjacency_matrix(s.automaton.G) + I,
    optimizer
)

function action(p::QCQPPolicy, (q0, x0), (qT, xT), P, H=nothing)

    @assert P[1] == q0
    @assert P[end] == qT

    if isnothing(H) H = ones(size(P)) end
    h = H / p.T

    A = HybridSystems.mode(p.s, 1).A
    B = HybridSystems.mode(p.s, 1).B
    nx, nu = size(B)
    M = length(P)

    m = Model(p.optimizer)
    @variable m x[1:M, 1:p.T, 1:nx]
    @variable m u[1:M, 1:p.T, 1:nu]
    
    @objective m Min sum(p.c(x[k,:,:], u[k,:,:]) * h[k] for k=1:M)

    @constraint m [k=1:M] x[k,2:end,:]' .== x[k,1:end-1,:]' + h[k]*A*x[k,1:end-1,:]' + h[k]*B*u[k,1:end-1,:]'
    @constraint m [k=1:M,t=1:p.T] p.K[P[k]][1] * x[k,t,:] .<= p.K[P[k]][2]
    @constraint m [k=2:M] x[k-1,end,:] .== x[k,1,:]
    @constraint m [k=2:M] u[k-1,end,:] .== u[k,1,:]
    @constraint m x[1,1,:] .== x0
    @constraint m x[end,end,:] .== xT

    optimize!(m)

    xr = reshape(value.(x), (M*p.T,nx))
    ur = reshape(value.(u), (M*p.T,nu))
    qr = repeat(P, p.T)

    return ur, (xr, qr), m

end
