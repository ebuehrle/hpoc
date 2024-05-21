using MomentOpt, DynamicPolynomials, MosekTools
using Plots
using DifferentialEquations

eout(E,i) = E[:,1] .== i
einc(E,i) = E[:,2] .== i

@polyvar x[1:4] u[1:2]
k1 = @set -1.0 <= x[1] && x[1] <= -0.7 && -1.0 <= x[2] && x[2] <= -0.7
k2 = @set -0.7 <= x[1] && x[1] <= -0.3 && -1.0 <= x[2] && x[2] <= -0.7
k3 = @set -0.3 <= x[1] && x[1] <= -0.0 && -1.0 <= x[2] && x[2] <= -0.7
k4 = @set -1.0 <= x[1] && x[1] <= -0.7 && -0.7 <= x[2] && x[2] <= -0.3
#k5 = @set -0.7 <= x[1] && x[1] <= -0.3 && -0.7 <= x[2] && x[2] <= -0.3
k6 = @set -0.3 <= x[1] && x[1] <= -0.0 && -0.7 <= x[2] && x[2] <= -0.3
k7 = @set -1.0 <= x[1] && x[1] <= -0.7 && -0.3 <= x[2] && x[2] <= -0.0
k8 = @set -0.7 <= x[1] && x[1] <= -0.3 && -0.3 <= x[2] && x[2] <= -0.0
k9 = @set -0.3 <= x[1] && x[1] <= -0.0 && -0.3 <= x[2] && x[2] <= -0.0

E = [
    10 1;
    
    1 4;
    4 7;
    7 8;
    8 9;
    
    4 1;
    7 4;
    8 7;
    9 8;
    
    1 2;
    2 3;
    3 6;
    6 9;
    
    2 1;
    3 2;
    6 3;
    9 6;

    9 11;
]

K = [
    k1,
    k1, k4, k7, k8,
    k4, k7, k8, k9,
    k1, k2, k3, k6,
    k2, k3, k6, k9,
    k9,
]

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]
c = x'*x + u'*u

f = A*x + B*u
μ0 = DiracMeasure(x, x0)
μT = DiracMeasure(x, xT)

m = GMPModel(Mosek.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 3)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f

modes = collect(Set(E) ∩ Set(1:9))
@variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i])
@objective m Min sum(Mom.(c,μ[:,2]))
cn = @constraint m [i=1:length(K)] Mom.(dbdt,μ[i,2]) .== Mom.(b,μ[i,3]) - Mom.(b,μ[i,1])
@constraint m [i=modes] sum(μ[eout(E,i),1]) == sum(μ[einc(E,i),3])
@constraint m sum(μ[eout(E,10),1]) == μ0
@constraint m sum(μ[einc(E,11),3]) == μT
@constraint m Mom.(1,μ[:,1]) .<= 1
@constraint m Mom.(1,μ[:,3]) .<= 1

optimize!(m)

p = integrate.(1,μ[:,1])
@show objective_value(m)
@show p

v = [first.(-dual.(c))'*b for c in cn]

contourf(range(-1.0,-0.7,100),range(-1.0,-0.7,100),(x1,x2)->v[2](x1,x2,0.0,0.5))
contourf!(range(-1.0,-0.7,100),range(-0.7,-0.3,100),(x1,x2)->v[3](x1,x2,0.0,0.5))
contourf!(range(-1.0,-0.7,100),range(-0.3,-0.0,100),(x1,x2)->v[4](x1,x2,0.0,0.5))
contourf!(range(-0.7,-0.3,100),range(-0.3,-0.0,100),(x1,x2)->v[5](x1,x2,0.0,0.5))
contourf!(range(-0.3,-0.0,100),range(-0.3,-0.0,100),(x1,x2)->v[end](x1,x2,0.0,0.5))

contourf!(range(-0.7,-0.3,100),range(-1.0,-0.7,100),(x1,x2)->v[11](x1,x2,0.0,0.5))
contourf!(range(-0.3,-0.0,100),range(-1.0,-0.7,100),(x1,x2)->v[12](x1,x2,0.0,0.5))
contourf!(range(-0.3,-0.0,100),range(-0.7,-0.3,100),(x1,x2)->v[13](x1,x2,0.0,0.5))

plot!([-0.7, -0.3, -0.3, -0.7, -0.7], [-0.7, -0.7, -0.3, -0.3, -0.7], color=:orange, linestyle=:dash, label=false)

fc(x0, p, t) = let 
    v0 = v[2]
    if -1.0 <= x0[2] <= -0.7 v0 = v[2] end
    if -0.7 <= x0[2] <= -0.3 v0 = v[3] end
    if -0.3 <= x0[2] <= -0.0 && -1.0 <= x0[1] <= -0.7 v0 = v[4] end
    if -0.3 <= x0[2] <= -0.0 && -0.7 <= x0[1] <= -0.3 v0 = v[5] end
    if -0.3 <= x0[2] <= -0.0 && -0.3 <= x0[1] <= -0.0 v0 = v[end] end
    dv0 = differentiate(v0, x)
    dvx0 = [d(x0) for d in dv0]
    u0 = -B'*dvx0
    A*x0 + B*u0
end
prob = ODEProblem(fc, x0, (0.0, 10.0), dt=0.05)
sol = solve(prob, Euler())
scatter!([Tuple(u[1:2]) for u in sol.u], color=:violet, label=false)
savefig("test/gmp2.pdf")
