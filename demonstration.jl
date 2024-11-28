using MomentOpt
using COSMO
using DynamicPolynomials
using DifferentialEquations
using Plots; default(fontfamily="Computer Modern", framestyle=:box);
using JLD2

@polyvar x[1:4] u[1:2]
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]

f = A*x + B*u
c = x'*x + u'*u + 0.1
µ0 = DiracMeasure(x, [-1, -0.5, 0, 0])
µT = DiracMeasure(x, [0.0, 0.0, 0, 0])

m = GMPModel(COSMO.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
mx = monomials(x,0:6)

@variable m µ1  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[1] <= -0.8)
@variable m µ10 Meas(x, support=@set x'*x <= 5 && x[1] <= -0.8)
@variable m µ1T Meas(x, support=@set x'*x <= 5 && x[1] <= -0.8)
c1 = @constraint m Mom.(differentiate(mx,x)*f,µ1) .== Mom.(mx,µ1T) - Mom.(mx,µ10)
@constraint m Mom.(mx,μ10) .== integrate.(mx,μ0)

@variable m µ2  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[2] <= -0.8)
@variable m µ20 Meas(x, support=@set x'*x <= 5 && x[2] <= -0.8)
@variable m µ2T Meas(x, support=@set x'*x <= 5 && x[2] <= -0.8)
c2 = @constraint m Mom.(differentiate(mx,x)*f,µ2) .== Mom.(mx,µ2T) - Mom.(mx,µ20)
@constraint m Mom.(mx,µ20) .== Mom.(mx,µ1T)

@variable m µ3  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[1] >= -0.2)
@variable m µ30 Meas(x, support=@set x'*x <= 5 && x[1] >= -0.2)
@variable m µ3T Meas(x, support=@set x'*x <= 5 && x[1] >= -0.2)
c3 = @constraint m Mom.(differentiate(mx,x)*f,µ3) .== Mom.(mx,µ3T) - Mom.(mx,µ30)
@constraint m Mom.(mx,µ30) .== Mom.(mx,µ2T)
@constraint m Mom.(mx,µ3T) .== integrate.(mx,µT)

@objective m Min Mom(c,µ1+µ2+µ3)
optimize!(m)

v1 = -first.(dual.(c1))'*mx
π1(s) = -0.5*B'*map(e->e(s),differentiate(v1,x))

v2 = -first.(dual.(c2))'*mx
π2(s) = -0.5*B'*map(e->e(s),differentiate(v2,x))

v3 = -first.(dual.(c3))'*mx
π3(s) = -0.5*B'*map(e->e(s),differentiate(v3,x))

g(s,p,t) = (A*s + B*π1(s))*all([-1.0,-0.7] .<= s[1:2] .<= [-0.6,-0.3]) +
    (A*s + B*π2(s))*all([-1.0,-1.0] .<= s[1:2] .<= [-0.4,-0.7]) +
    (A*s + B*π3(s))*all([-0.4,-1.0] .<= s[1:2] .<= [-0.0,-0.0])
X = solve(ODEProblem(g, [-1.0, -0.5, 0.0, 0.0], (0.0,20.0), dt=0.05), Euler())

contourf(range(-1.1,-0.7,100),range(-0.7,-0.3,100),(s...)->v1(s...,0,0))
contourf!(range(-1.1,-0.7,100),range(-0.3,0.1,100),(s...)->v1(s...,0,0))
contourf!(range(-1.1,-0.7,100),range(-1.1,-0.7,100),(s...)->v1(s...,0,0))
contourf!(range(-0.7,-0.3,100),range(-1.1,-0.7,100),(s...)->v2(s...,0,0))
contourf!(range(-0.3,0.1,100),range(-1.1,-0.7,100),(s...)->v3(s...,0,0))
contourf!(range(-0.3,0.1,100),range(-0.7,-0.3,100),(s...)->v3(s...,0,0))
contourf!(range(-0.3,0.1,100),range(-0.3,0.1,100),(s...)->v3(s...,0,0))
scatter!(stack(X.u)[1,:],stack(X.u)[2,:],color=1,label=false,grid=false,colorbar=false)
xlims!(-1.1,0.1)
ylims!(-1.1,0.1)
#savefig("out/demonstration.pdf")

mxu = monomials([x;u],0:6)
save("out/demonstration.jld2", 
    "monomials", mxu, 
    "moments", integrate.(mxu,µ1)+integrate.(mxu,µ2)+integrate.(mxu,µ3),
    "t", X.t,
    "u", X.u,
)
