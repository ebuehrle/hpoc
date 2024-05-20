using HybridSystems
using JuMP
using Graphs
using LinearAlgebra
using ProgressBars

function simulate(s::HybridSystem, Cx::Array, Cq::Array, Cu::Array, xo::Matrix, 
    qo::Matrix, uo::Matrix; dt::Float64, optimizer)

    K = [(stack([h.a for h in HybridSystems.mode(s, i).X.constraints])', 
          stack([h.b for h in HybridSystems.mode(s, i).X.constraints])) for i in 1:nmodes(s)]
    E = adjacency_matrix(s.automaton.G) + I

    A = HybridSystems.mode(s, 1).A
    B = HybridSystems.mode(s, 1).B
    nx, nu = size(B)

    @assert size(xo, 1) == size(qo, 1)
    @assert size(xo, 1) == size(uo, 1)
    T = size(xo, 1)

    if ndims(Cx) == 2 Cx = repeat(Cx, outer=(T,1,1)) end
    if ndims(Cq) == 2 Cq = repeat(Cq, outer=(T,1,1)) end
    if ndims(Cu) == 2 Cu = repeat(Cu, outer=(T,1,1)) end

    m = Model(optimizer)
    @variable m x[1:T,1:nx]
    @variable m u[1:T,1:nu]
    @variable m q[1:T,1:nmodes(s)] Bin

    @constraint m x[2:end,:]' .== x[1:end-1,:]' + dt*A*x[1:end-1,:]' + dt*B*u[1:end-1,:]'
    @constraint m [t=1:T,i=1:nmodes(s),j=1:size(K[i][1],1)] q[t,i] => {K[i][1][j,:]'*x[t,:] <= K[i][2][j]}
    @constraint m [t=1:T-1] q[t+1,:] .<= E'*q[t,:]
    @constraint m sum(q,dims=2) .== 1
    
    @constraint m [t=1:T] xo[t,:] .== Cx[t,:,:] * x[t,:]
    @constraint m [t=1:T] qo[t,:] .== Cq[t,:,:] * q[t,:]
    @constraint m [t=1:T] uo[t,:] .== Cu[t,:,:] * u[t,:]

    optimize!(m)

    return value.(x), value.(q), value.(u)

end

function euler(s::HybridSystem, x0::Vector, q0::Vector, u::AbstractMatrix; dt::Float64, optimizer)

    nx = size(x0, 1)
    nq = size(q0, 1)
    nu = size(u, 2)
    T = size(u, 1) + 1

    x = [x0'; zeros(T-1,nx)]
    q = [q0'; zeros(T-1,nq)]
    u = [u; zeros(1, nu)]

    Cx = stack([[I(nx)]; [zeros(nx, nx) for _ in 1:T-1]], dims=1)
    Cq = stack([[I(nq)]; [zeros(nq, nq) for _ in 1:T-1]], dims=1)
    Cu = stack([[I(nu) for _ in 1:T-1]; [zeros(nu, nu)]], dims=1)

    xs, qs, us = simulate(s, Cx, Cq, Cu, x, q, u, dt=dt, optimizer=optimizer)

    return xs, qs, us

end

euler(s::HybridSystem, x0::Vector, q0::Vector, u::Vector; dt::Float64, optimizer) = 
    euler(s, x0, q0, u', dt=dt, optimizer=optimizer)

struct EulerSimulator
    dt::Float64
    T::Int
    optimizer
end

function simulate(sim::EulerSimulator, s::HybridSystem, p::Function, (q0, x0)::Tuple{Vector, Vector})

    nx = length(x0)
    nq = length(q0)
    nu = length(p((q0, x0)))

    x = [x0'; zeros(sim.T-1,nx)]
    q = [q0'; zeros(sim.T-1,nq)]
    u = zeros(sim.T, nu)

    for t in ProgressBar(2:sim.T)
        u0 = p((q0, x0))
        x0, q0, _ = euler(s, x0, q0, u0, dt=sim.dt, optimizer=sim.optimizer)
        x0 = x0[end,:]
        q0 = q0[end,:]
        x[t,:] = x0
        q[t,:] = q0
        u[t,:] = u0
    end

    return x, q, u

end
