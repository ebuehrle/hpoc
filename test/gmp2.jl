using MomentOpt, DynamicPolynomials, MosekTools

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
    
    1 2;
    2 3;
    3 6;
    6 9;

    9 11;
]

K = [
    k1,
    k1, k4, k7, k8,
    k1, k2, k3, k6,
    k9,
]

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f = A*x + B*u
c = x'*x + u'*u
μ0 = DiracMeasure(x, [-1.0, -1.0, 0.0, 0.5])
μT = DiracMeasure(x, [-0.0, -0.0, 0.0, 0.0])

m = GMPModel(Mosek.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 3)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f

modes = collect(Set(E) ∩ Set(1:9))
@variable m μ1[i=1:length(K)] Meas([x;u], support=K[i])
@variable m μ2[i=1:length(K)] Meas([x;u], support=K[i])
@variable m μ3[i=1:length(K)] Meas([x;u], support=K[i])
@objective m Min sum(Mom.(c,μ2))
@constraint m [i=1:length(K)] Mom.(dbdt,μ2[i]) .== Mom.(b,μ3[i]) - Mom.(b,μ1[i])
@constraint m [i=modes] sum(μ1[eout(E,i)]) == sum(μ3[einc(E,i)])
@constraint m sum(μ1[eout(E,10)]) == μ0
@constraint m sum(μ3[einc(E,11)]) == μT

optimize!(m)

p = integrate.(1,μ1)
@show objective_value(m)
@show p
