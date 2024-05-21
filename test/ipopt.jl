using JuMP, Ipopt
using Plots

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
T = 10

m = Model(Ipopt.Optimizer)
@variable m x[1:5,1:T,1:4]
@variable m u[1:5,1:T,1:2]
@variable m h[1:5] .>= 0 start=1/T
@objective m Min sum(h[i]*(sum(x[i,:,:].^2) + sum(u[i,:,:].^2)) for i=1:5) + sum(h*T)
@constraint m [i=1:5] x[i,2:end,:]' .== x[i,1:end-1,:]' + h[i]*A*x[i,1:end-1,:]' + h[i]*B*u[i,1:end-1,:]'

@constraint m x[1,1,:] .== [-1, -1, 0, 0]
@constraint m x[1,end,:] .== x[2,1,:]
@constraint m x[2,end,:] .== x[3,1,:]
@constraint m x[3,end,:] .== x[4,1,:]
@constraint m x[4,end,:] .== x[5,1,:]
@constraint m x[5,end,:] .== [0, 0, 0, 0]

@constraint m [-1.0 -1.0] .<= x[1,:,1:2]; @constraint m x[1,:,1:2] .<= [-0.7 -0.7]
@constraint m [-1.0 -0.7] .<= x[2,:,1:2]; @constraint m x[2,:,1:2] .<= [-0.7 -0.3]
@constraint m [-1.0 -0.3] .<= x[3,:,1:2]; @constraint m x[3,:,1:2] .<= [-0.7 -0.0]
@constraint m [-0.7 -0.3] .<= x[4,:,1:2]; @constraint m x[4,:,1:2] .<= [-0.3 -0.0]
@constraint m [-0.3 -0.3] .<= x[5,:,1:2]; @constraint m x[5,:,1:2] .<= [-0.0 -0.0]

optimize!(m)

plot([-0.7, -0.3, -0.3, -0.7, -0.7], [-0.7, -0.7, -0.3, -0.3, -0.7], linestyle=:dash)
scatter!(value.(x[1,:,1]), value.(x[1,:,2]))
scatter!(value.(x[2,:,1]), value.(x[2,:,2]))
scatter!(value.(x[3,:,1]), value.(x[3,:,2]))
scatter!(value.(x[4,:,1]), value.(x[4,:,2]))
scatter!(value.(x[5,:,1]), value.(x[5,:,2]))
savefig("test/ipopt.pdf")
