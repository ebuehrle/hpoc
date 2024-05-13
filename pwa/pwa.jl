using Spot
using LazySets
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
    K = [(stack(c.a for c in h)', [c.b for c in h]) for h in K]

    e = [!isempty(HPolytope(k[1],k[2])) for k in K]
    K = K[e]
    E = E[e]

    q0 = [get_init_state_number(l)]
    qT = collect(reduce(∪, get_rabin_acceptance(l)[1]))

    return V, E, K, q0, qT
end

function pwa(A, B, f::Formula, translator::LTLTranslator)
    V, E, K, q0, qT = ltla(translator, f)
    M = [@system(x' = A*x + B*u, x ∈ HPolytope(a,b), u ∈ Universe(size(B,2))) for (a,b) in K]
    Σ = [e1[2] == e2[1] for e1 in E, e2 in E]

    automaton = GraphAutomaton(length(M))
    [add_transition!(automaton, i, j, 1) for i=1:size(Σ,1) for j=1:size(Σ,2) if Σ[i,j]]
    transitions = [@system(x' = x) for i=1:size(Σ,1) for j=1:size(Σ,2) if Σ[i,j]]
    switch = [AutonomousSwitching() for i=1:size(Σ,1) for j=1:size(Σ,2) if Σ[i,j]]
    h = HybridSystem(automaton, M, transitions, switch)

    q0 = [i for (i,e) in enumerate(E) if e[1] in q0]
    qT = [i for (i,e) in enumerate(E) if e[2] in qT]

    return h, q0, qT
end
