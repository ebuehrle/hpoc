using MomentOpt, DynamicPolynomials, Hypatia

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u
@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [-1.0, -1.0, 0.0, 0.5, 0.0, 0.0])
μT = DiracMeasure([x;u], [-0.0, -0.0, 0.0, 0.0, 0.0, 0.0])

K1 = @set(x[2] <= -0.7)
K2 = @set(x[1] >= -0.3)
K3 = @set(x[2] >= -0.3)
K4 = @set(x[1] <= -0.7)
K = [K1, K2, K3, K4]
E = [
    1 2;
    2 3;
    4 3;
    1 4;
]
eout(i) = E[:,1] .== i
einc(i) = E[:,2] .== i

m = GMPModel(Hypatia.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 2)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f(x,u)

@variable m μ[i=1:4,j=1:3] Meas([x;u], support=K[i])
@objective m Min sum(Mom.(c(x,u),μ[:,2]))
@constraint m [i=1:4] Mom.(dbdt, μ[i,2]) .== Mom.(b, μ[i,3]) - Mom.(b, μ[i,1])
@constraint m sum(μ[eout(1),1]) == μ0
@constraint m sum(μ[einc(3),3]) == μT
@constraint m [i in [2,4]] sum(μ[eout(i),1]) == sum(μ[einc(i),3])

optimize!(m)
