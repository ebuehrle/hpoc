using Spot
using LazySets
import Polyhedra: vrep, polyhedron, npoints, nrays, nlines
using HybridSystems, MathematicalSystems
using LinearAlgebra
include("formula.jl")

function split_edge_disjunctions(E, L)
    split_disjunctions(x) = strip.(split(x, "|"))
    remove_braces(x) = (x[1] == '(') ? x[2:end-1] : x
    split_conjunctions(x) = strip.(split(x, "&"))
    remove_negation(x) = x[2:end]
    group_by_signs(x) = (filter(l->!startswith(l,'!'),x), remove_negation.(filter(l->startswith(l,'!'),x)))

    split_label(x) = group_by_signs.(split_conjunctions.(remove_braces.(split_disjunctions(x))))
    S = split_label.(L)

    T = [[e for _ in 1:length(s)] for (e,s) in zip(E,S)]
    T = vcat(T...)
    S = vcat(S...)

    return T, S
end

function ltla(translator::LTLTranslator, f::Formula)
    println("constructing automaton")
    l, d = translate(translator, f)

    V = collect(1:num_states(l))
    E = get_edges(l)
    L = string.(get_labels(l))

    E, L = split_edge_disjunctions(E, L)

    K = [(Set(p), Set(n)) for (p,n) in L]
    K = [([d[l] for l in p], [d[l] for l in n]) for (p,n) in K]
    K = [(p, [HalfSpace(-h.a, -h.b) for h in n]) for (p,n) in K]
    K = [[Vector{HalfSpace{Float64, Vector{Float64}}}(p);Vector{HalfSpace{Float64, Vector{Float64}}}(n)] for (p,n) in K]
    K = [HPolyhedron(k) for k in K]
    
    println("removing empty modes")
    zerovol(p) = let n = length(first(p.constraints).a); let v = vrep(polyhedron(p ∩ Hyperrectangle(zeros(n),1000*ones(n)))); rank(v.V .- v.V[1,:]') <= n-1 end end
    e = [!isempty(k) && !zerovol(k) for k in K]
    E = E[e]
    K = K[e]
    
    println("removing redundant modes")
    r = [!any((e == xe) && issubset(k,xk) for (xk,xe) in zip(K[1:i-1],E[1:i-1])) for (i,(k,e)) in enumerate(zip(K,E))]
    K = K[r]
    E = E[r]
    
    println("removing redundant constraints")
    K = [remove_redundant_constraints(k) for k in K]

    K = [k.constraints for k in K]
    K = [(stack(c.a for c in h)', [c.b for c in h]) for h in K]

    q0 = [get_init_state_number(l)]
    qT = collect(reduce(∪, get_rabin_acceptance(l)[1]))

    return V, E, K, q0, qT
end

function pwa(A, B, f::Formula, x0, xT, translator::LTLTranslator)
    V, E, K, q0, qT = ltla(translator, f)
    println("computing pwa")
    M = [@system(x' = A*x + B*u, x ∈ HPolytope(a,b), u ∈ Universe(size(B,2))) for (a,b) in K]
    Σ = [e1[2] == e2[1] for e1 in E, e2 in E]
    sgl(m1,m2) = let n = length(first(m1).a); let v = vrep(polyhedron(HPolyhedron([m1;m2]) ∩ Hyperrectangle(zeros(n),1000*ones(n)))); rank(v.V .- v.V[1,:]') <= n-2 end end
    adj(m1,m2) = !isempty(HPolytope(m1) ∩ HPolytope(m2)) && !sgl(m1, m2)
    Σ = [(i != j) && Σ[i,j] && adj(M[i].X.constraints, M[j].X.constraints) for i=1:size(Σ,1), j=1:size(Σ,2)]

    automaton = GraphAutomaton(length(M))
    [add_transition!(automaton, i, j, 1) for i=1:size(Σ,1) for j=1:size(Σ,2) if Σ[i,j]]
    transitions = [@system(x' = x) for i=1:size(Σ,1) for j=1:size(Σ,2) if Σ[i,j]]
    switch = [AutonomousSwitching() for i=1:size(Σ,1) for j=1:size(Σ,2) if Σ[i,j]]
    h = HybridSystem(automaton, M, transitions, switch)

    q0 = [i for (i,e) in enumerate(E) if (e[1] in q0) && (x0 ∈ M[i].X)]
    qT = [i for (i,e) in enumerate(E) if (e[2] in qT) && (xT ∈ M[i].X)]

    return h, q0, qT
end

function fpwa(A, B, l::Formula, x0, xT, translator::LTLTranslator)
    hs, hq0, hqT = pwa(A, B, l, x0, xT, translator)
    
    println("computing fpwa")
    VV = [(source(hs,t),target(hs,t)) for t in HybridSystems.transitions(hs)]
    V = Tuple.(Set([Set(m) for m in VV]))
    M = [remove_redundant_constraints(HPolyhedron([
        HybridSystems.mode(hs,i).X.constraints; 
        HybridSystems.mode(hs,j).X.constraints])) for (i,j) in V]
    Ex = [((i1,Set(v1)),(i2,Set(v2))) for (i1,v1) in enumerate(V) for (i2,v2) in enumerate(V) if (i1 != i2) && !isempty(Set(v1) ∩ Set(v2))]
    E = [(i1,i2) for ((i1,_),(i2,_)) in Ex]
    K = [HybridSystems.mode(hs, first(v1 ∩ v2)).X for ((_,v1),(_,v2)) in Ex]
    qm0 = [(i,first(Set(hq0) ∩ Set(v))) for (i,v) in enumerate(V) if !isempty(Set(hq0) ∩ Set(v)) && x0 in HybridSystems.mode(hs,first(Set(hq0) ∩ Set(v))).X]
    q0 = [i for (i,_) in qm0]
    K0 = [HybridSystems.mode(hs,k).X for (_,k) in qm0]
    qmT = [(i,first(Set(hqT) ∩ Set(v))) for (i,v) in enumerate(V) if !isempty(Set(hqT) ∩ Set(v)) && xT in HybridSystems.mode(hs,first(Set(hqT) ∩ Set(v))).X]
    qT = [i for (i,_) in qmT]
    KT = [HybridSystems.mode(hs,k).X for (_,k) in qmT]

    automaton = GraphAutomaton(length(M))
    [add_transition!(automaton, i, j, k) for (k,(i,j)) in enumerate(E)]
    modes = [@system(x' = 0, x ∈ m) for m in M]
    transitions = [@system(x' = A*x + B*u, x ∈ k, u ∈ Universe(size(B,2))) for k in K]
    switchings = [AutonomousSwitching() for _ in transitions]
    h = HybridSystem(automaton, modes, transitions, switchings)

    return h, q0, qT, K0, KT
end
