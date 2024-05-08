using MomentOpt, DynamicPolynomials, Clarabel

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
f(x,u) = A*x + B*u
c(x,u) = x'*x + u'*u
@polyvar x[1:4] u[1:2]
μ0 = DiracMeasure([x;u], [-1.0, -1.0, 0.0, 0.5, 0.0, 0.0])
μT = DiracMeasure([x;u], [-0.0, -0.0, 0.0, 0.0, 0.0, 0.0])

S1 = @set(x[1] <= -0.7)
S2 = @set(-0.3 <= x[2])

m = GMPModel(Clarabel.Optimizer)
set_approximation_mode(m, PRIMAL_RELAXATION_MODE())
set_approximation_degree(m, 2)
b = monomials(x, 0:approximation_degree(m))
dbdt = differentiate(b, x) * f(x,u)

@variable m μ10 Meas([x;u], support=S1)
@variable m μ1f Meas([x;u], support=S1)
@variable m μ1T Meas([x;u], support=S1)
@variable m μ20 Meas([x;u], support=S2)
@variable m μ2f Meas([x;u], support=S2)
@variable m μ2T Meas([x;u], support=S2)

@objective m Min Mom(c(x,u),μ1f) + Mom(c(x,u),μ2f)
@constraint m μ10 == μ0
@constraint m Mom.(dbdt, μ1f) .== Mom.(b, μ1T) - Mom.(b, μ10)
@constraint m μ1T == μ20
@constraint m Mom.(dbdt, μ2f) .== Mom.(b, μ2T) - Mom.(b, μ20)
@constraint m μ2T == μT

optimize!(m)
write("img/gmp2.txt","$(objective_value(m))")
