include("PWA.jl/PWA.jl")
using Symbolics, LazySets
using HybridSystems
using JuMP
using MosekTools
using Ipopt
using Plots
using .PWA

c(x,u) = x'*x + u'*u
A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]
B = [0 0; 0 0; 1 0; 0 1]
l = let x = Symbolics.variables(:x, 1:4)
    o1 = And_([HalfSpace(-0.4 <= x[1], x), HalfSpace(x[1] <= -0.3, x), HalfSpace(-0.4 <= x[2], x), HalfSpace(x[2] <= -0.3, x)])
    o2 = And_([HalfSpace( 0.3 <= x[1], x), HalfSpace(x[1] <=  0.4, x), HalfSpace(-0.4 <= x[2], x), HalfSpace(x[2] <= -0.3, x)])
    o3 = And_([HalfSpace( 0.3 <= x[1], x), HalfSpace(x[1] <=  0.4, x), HalfSpace( 0.3 <= x[2], x), HalfSpace(x[2] <=  0.4, x)])
    o4 = And_([HalfSpace(-0.4 <= x[1], x), HalfSpace(x[1] <= -0.3, x), HalfSpace( 0.3 <= x[2], x), HalfSpace(x[2] <=  0.4, x)])

    t1 = And_([HalfSpace(-0.4 <= x[1], x), HalfSpace(x[1] <= -0.3, x), HalfSpace(-0.05 <= x[2], x), HalfSpace(x[2] <= 0.05, x)])
    t2 = And_([HalfSpace( 0.3 <= x[1], x), HalfSpace(x[1] <=  0.4, x), HalfSpace(-0.05 <= x[2], x), HalfSpace(x[2] <= 0.05, x)])

    c = And_([HalfSpace(-0.2 <= x[1], x), HalfSpace(x[1] <= 0.2, x), HalfSpace(-0.2 <= x[2], x), HalfSpace(x[2] <= 0.2, x)])

    F(o1) & F(o2) & F(o3) & F(o4) & F(t1) & F(t2)
end
x0 = [0.0, 0.0, 0.0, 0.0]
xT = [0.0, 0.0, 0.0, 0.0]

s, q0, qT = PPWA(A, B, l, merge_modes=false)
println(HybridSystems.nmodes(s), " modes")
println(HybridSystems.ntransitions(s), " transitions")
println(q0)
println(qT)

policy = GMPPolicy(s, c; optimizer=Mosek.Optimizer)
uq, xq, qq, mq, m = extract(policy, (q0, x0), (qT, xT); T=20, optimizer=Ipopt.Optimizer)

scatter(xq[:,1],xq[:,2],label="J = $(round(objective_value(mq), digits=2)) ($(round(objective_value(m), digits=2)))")

o1 = plot([-0.4, -0.3, -0.3, -0.4, -0.4], [-0.4, -0.4, -0.3, -0.3, -0.4], color=:green, fill=true, fillalpha=0.2)
o2 = plot([ 0.3,  0.4,  0.4,  0.3,  0.3], [-0.4, -0.4, -0.3, -0.3, -0.4], color=:green, fill=true, fillalpha=0.2)
o3 = plot([ 0.3,  0.4,  0.4,  0.3,  0.3], [ 0.3,  0.3,  0.4,  0.4,  0.3], color=:green, fill=true, fillalpha=0.2)
o4 = plot([-0.4, -0.3, -0.3, -0.4, -0.4], [ 0.3,  0.3,  0.4,  0.4,  0.3], color=:green, fill=true, fillalpha=0.2)

t1 = plot([-0.4, -0.3, -0.3, -0.4, -0.4], [-0.05, -0.05, 0.05, 0.05, -0.05], color=:yellow, fill=true, fillalpha=0.2)
t2 = plot([ 0.3,  0.4,  0.4,  0.3,  0.3], [-0.05, -0.05, 0.05, 0.05, -0.05], color=:yellow, fill=true, fillalpha=0.2)

c = plot([-0.2, 0.2, 0.2, -0.2, -0.2], [-0.2, -0.2, 0.2, 0.2, -0.2], color=:blue, fill=true, fillalpha=0.2)

savefig("img/rover-1.pdf")
