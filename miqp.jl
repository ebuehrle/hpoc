using JuMP, Gurobi
using Symbolics, LazySets, SemialgebraicSets
using Graphs, LinearAlgebra
using Plots
include("pwa/product.jl")

struct MIQPPolicy
    s::HybridSystem
    c::Function
    h::Float64
    T::Integer
    K::Vector{Tuple{Matrix, Array}}
    E::Matrix
end

MIQPPolicy(s, c; h = 0.1, T = 50) = MIQPPolicy(
    s,
    c,
    h,
    T,
    [(stack([h.a for h in HybridSystems.mode(s,i).X.constraints])', 
      stack([h.b for h in HybridSystems.mode(s,i).X.constraints])) for i in 1:nmodes(s)],
    adjacency_matrix(s.automaton.G) + I
)

function action(p::MIQPPolicy, (q0,x0), (qT,xT))

    if !(q0 isa Array) q0 = [q0] end
    if !(qT isa Array) qT = [qT] end

    A = HybridSystems.mode(p.s, 1).A
    B = HybridSystems.mode(p.s, 1).B
    nx, nu = size(B)

    m = Model(Gurobi.Optimizer)
    @variable m x[1:p.T,1:nx]
    @variable m u[1:p.T,1:nu]
    @variable m q[1:p.T,1:nmodes(p.s)] Bin
    @objective m Min p.c(x,u) * p.h
    @constraint m x[1,:] .== x0
    @constraint m x[end,:] .== xT
    @constraint m x[2:end,:]' .== x[1:end-1,:]' + p.h*A*x[1:end-1,:]' + p.h*B*u[1:end-1,:]'
    @constraint m [t=1:p.T,i=1:nmodes(p.s),j=1:size(p.K[i][1],1)] q[t,i] => {p.K[i][1][j,:]'*x[t,:] <= p.K[i][2][j]}
    @constraint m [t=1:p.T-1] q[t+1,:] .<= p.E'*q[t,:]
    @constraint m sum(q,dims=2) .== 1
    @constraint m q[1,:] .<= [q in q0 for q=1:nmodes(p.s)]
    @constraint m q[end,:] .<= [q in qT for q=1:nmodes(p.s)]

    optimize!(m)

    return value.(u), (value.(x), value.(q)), m

end

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))
end

h, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(h), " modes")
println(HybridSystems.ntransitions(h), " transitions")

c(x,u) = sum(x.^2) + sum(u.^2)
policy = MIQPPolicy(h, c, h=0.1, T=50)

x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]
(u, (x, q), m) = action(policy, (q0, x0), (qT, xT))
scatter(x[:,1], x[:,2], label=objective_value(m))
savefig("img/miqp.pdf")
