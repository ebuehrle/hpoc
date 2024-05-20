using MomentOpt, DynamicPolynomials, Clarabel
using Plots

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u
x0 = [-0.8, -0.8, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]

@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [x0; 0.0; 0.0])
μT = DiracMeasure([x;u], [xT; 0.0; 0.0])

K1 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && x[1] <= -0.7 && x[2] <= -0.7
K2 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && -0.7 <= x[1] && x[1] <= -0.3 && x[2] <= -0.7
K3 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && -0.3 <= x[1] && x[2] <= -0.7
K4 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && x[1] <= -0.7 && -0.7 <= x[2] && x[2] <= -0.3
K5 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && -0.3 <= x[1] && -0.7 <= x[2] && x[2] <= -0.3
K6 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && x[1] <= -0.7 && -0.3 <= x[2]
K7 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && -0.7 <= x[1] && x[1] <= -0.3 && -0.3 <= x[2]
K8 = @set -1.0 <= x[1] && x[1] <= 0.0 && -1.0 <= x[2] && x[2] <= 0.0 && -0.3 <= x[1] && -0.3 <= x[2]

K12 = @set K1 && K2
K23 = @set K2 && K3
K14 = @set K1 && K4
K35 = @set K3 && K5
K46 = @set K4 && K6
K58 = @set K5 && K8
K67 = @set K6 && K7
K78 = @set K7 && K8

E = [
    1 2;
    2 1;
    
    2 3;
    3 2;
    
    3 4;
    4 3;
    
    4 5;
    5 4;
    
    5 6;
    6 5;
    
    6 7;
    7 6;
    
    7 8;
    8 7;
    
    8 1;
    1 8;
    
    9 1;
    9 8;
    
    4 10;
    5 10;
]
eout(i) = E[:,1] .== i
einc(i) = E[:,2] .== i

K = [
    [K12, K2, K23],
    [K23, K2, K12],

    [K23, K3, K35],
    [K35, K3, K23],

    [K35, K5, K58],
    [K58, K5, K35],

    [K58, K8, K78],
    [K78, K8, K58],

    [K78, K7, K67],
    [K67, K7, K78],

    [K67, K6, K46],
    [K46, K6, K67],

    [K46, K4, K14],
    [K14, K4, K46],

    [K14, K1, K12],
    [K12, K1, K14],

    [K1, K1, K12],
    [K1, K1, K14],

    [K58, K8, K8],
    [K78, K8, K8],
]

_modes = 8

m = GMPModel(Clarabel.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 2)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f(x,u)

@variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i][j])
@objective m Min sum(Mom.(c(x,u), μ[:,2])) + 0.01*sum(Mom.(1, μ))
cn = @constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
@constraint m [i=1:_modes] sum(μ[eout(i),1]) == sum(μ[einc(i),3])
@constraint m sum(μ[eout(_modes+1),1]) == μ0
@constraint m sum(μ[einc(_modes+2),3]) == μT
@constraint m Mom.(1, μ[:,1]) .<= 1
@constraint m Mom.(1, μ[:,3]) .<= 1

optimize!(m)

p = integrate.(1,μ[:,1])
@constraint m Mom.(1, μ[:,1]) .== round.(p)
optimize!(m)

println(sum(integrate.(c(x,u), μ[:,2])))
t = integrate.(1,μ[:,2])
o = integrate.(1,μ[:,3])
J = integrate.(c(x,u),μ[:,2])
r = collect(zip(Tuple.(eachrow(E)),p,t,o,J))
[println(l) for l in r]
#println(p)
#println(J)

v = [first.(-dual.(c))'*b for c in cn]
contourf(range(-1,0,100),range(-1,0,100),(x1,x2)->v[1](x1,x2,0,0.5))
