using JuMP, HiGHS

w = [4, 1, 2, 2]
e = [
    1 2;
    2 3;
    4 3;
    1 4;
]
out = stack([e[:,1] .== i for i=1:size(e,1)])'
inc = stack([e[:,2] .== i for i=1:size(e,1)])'

m = Model(HiGHS.Optimizer)
@variable m x[1:size(e,1)]
@constraint m 0 .<= x
@constraint m x .<= 1
@objective m Min w'*x
@constraint m sum(x[out[1,:]]) == 1
@constraint m sum(x[inc[3,:]]) == 1
@constraint m [i in [2,4]] sum(x[out[i,:]]) == sum(x[inc[i,:]])

optimize!(m)
println(objective_value(m))
