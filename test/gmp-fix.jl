using MomentOpt, DynamicPolynomials, Clarabel

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u
x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]

@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [x0; 0.0; 0.0])
μT = DiracMeasure([x;u], [xT; 0.0; 0.0])

K1 = @set(x[1] <= -0.7 && x[2] <= -0.7)
K2 = @set(x[1] <= -0.7 && x[2] >= -0.7 && x[2] <= -0.3)
K12 = @set(x[1] <= -0.7 && x[2] <= -0.7 && x[2] >= -0.7)
K3 = @set(x[1] <= -0.7 && x[2] >= -0.3)
K23 = @set(x[1] <= -0.7 && x[2] <= -0.3 && x[2] >= -0.3)
K4 = @set(x[1] >= -0.7 && x[1] <= -0.3 && x[2] >= -0.3)
K34 = @set(x[1] <= -0.7 && x[1] >= -0.7 && x[2] >= -0.3)
K5 = @set(x[1] >= -0.3 && x[2] >= -0.3)
K45 = @set(x[1] <= -0.3 && x[1] >= -0.3 && x[2] >= -0.3)
K = [
    [K1, K1, K12],
    [K12, K2, K23],
    [K23, K3, K34],
    [K34, K4, K45],
    [K45, K5, K5],
]

m = GMPModel(Clarabel.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 2)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f(x,u)

@variable m μ[i=1:length(K),j=1:3] Meas([x;u], support=K[i][j])
@objective m Min sum(Mom.(c(x,u), μ[:,2]))
@constraint m [i=1:length(K)] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
@constraint m μ[1,3] == μ[2,1]
@constraint m μ[2,3] == μ[3,1]
@constraint m μ[3,3] == μ[4,1]
@constraint m μ[4,3] == μ[5,1]
@constraint m μ[1,1] == μ0
@constraint m μ[5,3] == μT

optimize!(m)
println(objective_value(m))
