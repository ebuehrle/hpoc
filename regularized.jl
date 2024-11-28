using MomentOpt
using COSMO
using DynamicPolynomials
using DifferentialEquations
using Plots; default(fontfamily="Computer Modern", framestyle=:box);
using JLD2
using MultivariateMoments

@polyvar x[1:4] u[1:2]
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]

f = A*x + B*u
c = x'*x + u'*u + 0.1
µ0 = DiracMeasure(x, [-1, -0.5, 0, 0])
µT = DiracMeasure(x, [0.0, 0.0, 0, 0])
μGT_moments = load("out/demonstration.jld2", "moments")
μGT_monomials = monomials([x;u],0:6)
μGT = moment_measure(MultivariateMoments.Measure(μGT_moments, μGT_monomials))

m = GMPModel(COSMO.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
mx = monomials(x,0:6)

@variable m µ1  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[1] <= -0.7)
@variable m µ10 Meas(x, support=@set x'*x <= 5 && x[1] <= -0.7)
@variable m µ1T Meas(x, support=@set x'*x <= 5 && x[1] <= -0.7)
c1 = @constraint m Mom.(differentiate(mx,x)*f,µ1) .== Mom.(mx,µ1T) - Mom.(mx,µ10)

@variable m µ14  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[1] <= -0.7)
@variable m µ140 Meas(x, support=@set x'*x <= 5 && x[1] <= -0.7)
@variable m µ14T Meas(x, support=@set x'*x <= 5 && x[1] <= -0.7)
c14 = @constraint m Mom.(differentiate(mx,x)*f,µ14) .== Mom.(mx,µ14T) - Mom.(mx,µ140)
@constraint m Mom.(mx,µ10) + Mom.(mx,µ140) .== integrate.(mx,µ0)

@variable m µ2  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[2] <= -0.7)
@variable m µ20 Meas(x, support=@set x'*x <= 5 && x[2] <= -0.7)
@variable m µ2T Meas(x, support=@set x'*x <= 5 && x[2] <= -0.7)
c2 = @constraint m Mom.(differentiate(mx,x)*f,µ2) .== Mom.(mx,µ2T) - Mom.(mx,µ20)
@constraint m Mom.(mx,µ20) .== Mom.(mx,µ1T)

@variable m µ4  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[2] >= -0.3)
@variable m µ40 Meas(x, support=@set x'*x <= 5 && x[2] >= -0.3)
@variable m µ4T Meas(x, support=@set x'*x <= 5 && x[2] >= -0.3)
c4 = @constraint m Mom.(differentiate(mx,x)*f,µ4) .== Mom.(mx,µ4T) - Mom.(mx,µ40)
@constraint m Mom.(mx,µ40) .== Mom.(mx,µ14T)

@variable m µ3  Meas([x;u], support=@set [x;u]'*[x;u] <= 5 && x[1] >= -0.3)
@variable m µ30 Meas(x, support=@set x'*x <= 5 && x[1] >= -0.3)
@variable m µ3T Meas(x, support=@set x'*x <= 5 && x[1] >= -0.3)
c3 = @constraint m Mom.(differentiate(mx,x)*f,µ3) .== Mom.(mx,µ3T) - Mom.(mx,µ30)
@constraint m Mom.(mx,µ30) .== Mom.(mx,µ2T) + Mom.(mx,µ4T)
@constraint m Mom.(mx,µ3T) .== integrate.(mx,µT)

@polyvar e
@variable m μe[1:length(mx)] Meas([e], support=@set e'*e <= 1)
ce = @constraint m Mom.(differentiate(mx,x)*f,µ1+µ2+µ3+μ14+μ4) + Mom.(e,μe) .== integrate.(differentiate(mx,x)*f,μGT)
@constraint m Mom.(1,μe) .== 1
@constraint m [i=1:length(μe)] Mom.(monomials([e],3:6),μe[i]) .== 0

@objective m Min Mom(c,µ1+µ2+µ3+μ14+μ4) + 10^4*Mom(e^2,sum(μe))
optimize!(m)

@show integrate(1,µ10)
@show integrate(1,µ140)

ve = -first.(dual.(ce))'*mx

v1 = -first.(dual.(c1))'*mx + ve
π1(s) = -0.5*B'*map(e->e(s),differentiate(v1,x))

v14 = -first.(dual.(c14))'*mx + ve
π14(s) = -0.5*B'*map(e->e(s),differentiate(v14,x))

v2 = -first.(dual.(c2))'*mx + ve
π2(s) = -0.5*B'*map(e->e(s),differentiate(v2,x))

v3 = -first.(dual.(c3))'*mx + ve
π3(s) = -0.5*B'*map(e->e(s),differentiate(v3,x))

v4 = -first.(dual.(c4))'*mx + ve
π4(s) = -0.5*B'*map(e->e(s),differentiate(v4,x))

g(s,p,t) = (A*s + B*π1(s))*all([-1.1,-0.7] .<= s[1:2] .<= [-0.6,-0.3]) +
    (A*s + B*π2(s))*all([-1.0,-1.0] .<= s[1:2] .<= [-0.4,-0.7]) +
    (A*s + B*π3(s))*all([-0.4,-1.0] .<= s[1:2] .<= [-0.0,-0.0]) +
    (A*s + B*π4(s))*all([-1.0,-0.3] .<= s[1:2] .<= [-0.4,-0.0])
X = solve(ODEProblem(g, [-1.0, -0.5, 0.0, 0.0], (0.0,20.0), dt=0.05), Euler())

D = load("out/demonstration.jld2", "u")

contourf(range(-1.1,-0.7,100),range(-0.7,-0.3,100),(s...)->v1(s...,0,0))
contourf!(range(-1.1,-0.7,100),range(-0.3,0.1,100),(s...)->v1(s...,0,0))
contourf!(range(-1.1,-0.7,100),range(-1.1,-0.7,100),(s...)->v1(s...,0,0))
contourf!(range(-0.7,-0.3,100),range(-1.1,-0.7,100),(s...)->v2(s...,0,0))
contourf!(range(-0.3,0.1,100),range(-1.1,-0.7,100),(s...)->v3(s...,0,0))
contourf!(range(-0.3,0.1,100),range(-0.7,-0.3,100),(s...)->v3(s...,0,0))
contourf!(range(-0.3,0.1,100),range(-0.3,0.1,100),(s...)->v3(s...,0,0))
contourf!(range(-0.7,-0.3,100),range(-0.3,0.1,100),(s...)->v4(s...,0,0))
scatter!(stack(X.u)[1,:],stack(X.u)[2,:],color=1,label="plan",grid=false,colorbar=false)
scatter!(stack(D)[1,:],stack(D)[2,:],color=2,alpha=0.5,label="demonstration",grid=false,colorbar=false)
xlims!(-1.1,0.1)
ylims!(-1.1,0.1)
savefig("out/regularized-wide.pdf")

plot!(ratio=:equal)
savefig("out/regularized.pdf")
