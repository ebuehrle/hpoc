using JuMP, Gurobi
using Symbolics, LazySets, SemialgebraicSets
using Graphs, LinearAlgebra
using Plots
include("pwa/pwa.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
h = 0.1
T = 50
x0 = [-10.0, -2.0, 10.0, 0.0]
xT = [-0.0, -0.0, 0.0, 0.0]

Symbolics.@variables x[1:4]
lane1 = HalfSpace(x[2] >= -3.0, x) & HalfSpace(x[2] <= -1.0, x) & HalfSpace(x[1] <= -5.0, x)
lane2 = HalfSpace(x[2] >= -1.0, x) & HalfSpace(x[2] <=  1.0, x)
lf = G(lane1 | lane2)

hs, q0, qT = pwa(A, B, lf, x0, xT, LTLTranslator())
K = [(stack([h.a for h in HybridSystems.mode(hs,i).X.constraints])', 
      stack([h.b for h in HybridSystems.mode(hs,i).X.constraints])) for i in 1:nmodes(hs)]
E = adjacency_matrix(hs.automaton.G) + I

m = Model(Gurobi.Optimizer)
@variable m x[1:T,1:4]
@variable m u[1:T,1:2]
@variable m q[1:T,1:nmodes(hs)] Bin
@objective m Min sum(x.^2)*h + sum(u.^2)*h
@constraint m x[1,:] .== x0
@constraint m x[end,:] .== xT
@constraint m x[2:end,:]' .== x[1:end-1,:]' + h*A*x[1:end-1,:]' + h*B*u[1:end-1,:]'
@constraint m [t=1:T,i=1:nmodes(hs),j=1:size(K[i][1],1)] q[t,i] => {K[i][1][j,:]'*x[t,:] <= K[i][2][j]}
@constraint m [t=1:T-1] q[t+1,:] .<= E'*q[t,:]
@constraint m sum(q,dims=2) .== 1
@constraint m q[1,:] .<= [q in q0 for q=1:nmodes(hs)]
@constraint m q[end,:] .<= [q in qT for q=1:nmodes(hs)]

optimize!(m)
scatter(value.(x[:,1]), value.(x[:,2]), label=objective_value(m))
savefig("img/miqp-merge.pdf")
