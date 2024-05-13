using Spot
using LazySets, Polyhedra
using HybridSystems, MathematicalSystems
include("formula.jl")

function ltla(translator::LTLTranslator, f::Formula)
    l, d = translate(translator, f)
    l = split_edges(l)

    V = collect(1:num_states(l))
    E = get_edges(l)

    K = [(Set(label_to_array(s)), setdiff(Set(atomic_propositions(l)), Set(label_to_array(s)))) for s in get_labels(l)]
    K = [([d[string(l)] for l in p], [d[string(l)] for l in n]) for (p,n) in K]
    K = [(p, [HalfSpace(-h.a, -h.b) for h in n]) for (p,n) in K]
    K = [[p;n] for (p,n) in K]
    K = [remove_redundant_constraints(HPolyhedron(k)) for k in K]
    
    e = [!isempty(k) for k in K]
    K = K[e]
    E = E[e]

    K = [k.constraints for k in K]
    K = [(stack(c.a for c in h)', [c.b for c in h]) for h in K]

    q0 = [get_init_state_number(l)]
    qT = collect(reduce(∪, get_rabin_acceptance(l)[1]))

    return V, E, K, q0, qT
end

function pwa(A, B, f::Formula, x0, xT, translator::LTLTranslator)
    V, E, K, q0, qT = ltla(translator, f)
    M = [@system(x' = A*x + B*u, x ∈ HPolytope(a,b), u ∈ Universe(size(B,2))) for (a,b) in K]
    Σ = [e1[2] == e2[1] for e1 in E, e2 in E]
    sgl(m1,m2) = let n = length(first(m1).a); let v = vrep(polyhedron(HPolyhedron([m1;m2]))); !(npoints(v) >= n || nrays(v) + nlines(v) >= n-1) end end
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
    hs, q0, qT = pwa(A, B, l, x0, xT, translator)

    M = [remove_redundant_constraints(HPolyhedron([
        HybridSystems.mode(hs,source(hs,r)).X.constraints; 
        HybridSystems.mode(hs,target(hs,r)).X.constraints])) for r in HybridSystems.transitions(hs)]
    E = [(i1,i2) for (i1,r1) in enumerate(HybridSystems.transitions(hs)) for (i2,r2) in enumerate(HybridSystems.transitions(hs)) if target(hs,r1) == source(hs,r2) && source(hs,r1) != target(hs,r2)]
    K = [HybridSystems.mode(hs,target(hs,r1)).X for (i1,r1) in enumerate(HybridSystems.transitions(hs)) for (i2,r2) in enumerate(HybridSystems.transitions(hs)) if (i1,i2) in E]
    q0 = [i for (i,e) in enumerate(HybridSystems.transitions(hs)) if source(hs,e) in q0]
    K0 = [HybridSystems.mode(hs,source(hs,t)).X for (i,t) in enumerate(HybridSystems.transitions(hs)) if i in q0]
    qT = [i for (i,e) in enumerate(HybridSystems.transitions(hs)) if target(hs,e) in qT]
    KT = [HybridSystems.mode(hs,target(hs,t)).X for (i,t) in enumerate(HybridSystems.transitions(hs)) if i in qT]

    automaton = GraphAutomaton(length(M))
    [add_transition!(automaton, i, j, k) for (k,(i,j)) in enumerate(E)]
    modes = [@system(x' = 0, x ∈ m) for m in M]
    transitions = [@system(x' = A*x + B*u, x ∈ k, u ∈ Universe(size(B,2))) for k in K]
    switchings = [AutonomousSwitching() for _ in transitions]
    h = HybridSystem(automaton, modes, transitions, switchings)

    return h, q0, qT, K0, KT
end
