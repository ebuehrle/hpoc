using JuMP, Gurobi
using Symbolics, LazySets, SemialgebraicSets
using Graphs, LinearAlgebra
using Plots
include("pwa/product.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
h = 0.3
T = 30
x0 = [-0.4, 0.4, 0.0, 0.0]
xT = [-0.0, 0.0, 0.0, 0.0]

Symbolics.@variables x[1:4]
k1 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.4 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
d1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.2 <= x[2], x) & HalfSpace(x[2] <= 0.3, x)
k2 = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.5, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.1, x)
d2 = HalfSpace(-0.2 <= x[1], x) & HalfSpace(x[1] <= -0.1, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
w1 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.2, x)
w2 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
w3 = HalfSpace(-0.3 <= x[1], x) & HalfSpace(x[1] <= -0.2, x) & HalfSpace(0.3 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
room = HalfSpace(-0.6 <= x[1], x) & HalfSpace(x[1] <= -0.0, x) & HalfSpace(0.0 <= x[2], x) & HalfSpace(x[2] <= 0.5, x)
lf = G(room) & U(!d1,k1)# & U(!d2,k2) & G(!w1) & G(!w2) & G(!w3)

hs, q0, qT = ppwa(A, B, lf, x0, xT, LTLTranslator())
println(HybridSystems.nmodes(hs), " modes")
println(HybridSystems.ntransitions(hs), " transitions")
println("initial states ", q0)
println("terminal states ", qT)
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
savefig("img/miqp-doors.pdf")
