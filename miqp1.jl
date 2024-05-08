using JuMP, Gurobi
using Plots

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
h = 0.1
T = 100

m = Model(Gurobi.Optimizer)
@variable m x[1:T,1:4]
@variable m u[1:T,1:2]
@variable m q[1:T,1:2] Bin
@constraint m x[2:end,:]' .== x[1:end-1,:]' + h*A*x[1:end-1,:]' + h*B*u[1:end-1,:]'
@constraint m x[1,:] .== [-1.0, -1.0, 0.0, 0.5]
@constraint m x[end,:] .== [0.0, 0.0, 0.0, 0.0]
@constraint m q[:,1] .=> {x[:,2] .<= -0.7}
@constraint m q[:,2] .=> {-0.3 .<= x[:,1]}
@constraint m sum(q, dims=2) .== 1
@objective m Min sum(x.^2)*h + sum(u.^2)*h

optimize!(m)

scatter(value.(x[:,1]), value.(x[:,2]), label="J=$(objective_value(m))")
savefig("img/miqp1.pdf")
