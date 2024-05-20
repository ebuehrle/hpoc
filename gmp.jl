using Symbolics, LazySets
using HybridSystems
using MosekTools, Gurobi
using Plots
include("pwa/product.jl")
include("pwa/gmp.jl")
include("pwa/simulate.jl")

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    G(!(HalfSpace(x[1] >= -0.7, x) & HalfSpace(x[1] <= -0.3, x) 
    & HalfSpace(x[2] >= -0.7, x) & HalfSpace(x[2] <= -0.3, x)))
end

s, q0, qT = PPWA(A, B, l)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

c(x,u) = x'*x + u'*u
policy = GMPPolicy(s, c, optimizer=Mosek.Optimizer)

x0 = [-1.0, -1.0, 0.0, 0.5]
xT = [-0.0, -0.0, 0.0, 0.0]
_, m, p = action(policy, (q0,x0), (qT,xT))

q0 = [(i ∈ q0) && (x0 ∈ HybridSystems.mode(s,i).X) for i in 1:nmodes(s)]
qT = [(i ∈ qT) && (xT ∈ HybridSystems.mode(s,i).X) for i in 1:nmodes(s)]
x, q, u = simulate(
    EulerSimulator(0.05, 15, Gurobi.Optimizer),
    s,
    ((q0,x0),) -> let (dv, _, _) = action(policy, (argmax(q0), x0), (argmax(qT), xT)); -B'*dv end,
    (q0, x0)
)

scatter(x[:,1], x[:,2], label=objective_value(m))
savefig("img/gmp.pdf")
write("img/gmp.txt","$(objective_value(m))")
write("img/gmpp.txt","$(p)")
