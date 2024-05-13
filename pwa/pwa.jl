using Spot
using LazySets
include("formula.jl")

function pwa(translator::LTLTranslator, f::Formula)
    l, d = translate(translator, f)
    l = split_edges(l)

    V = collect(1:num_states(l))
    E = get_edges(l)

    K = [(Set(label_to_array(s)), setdiff(Set(atomic_propositions(l)), Set(label_to_array(s)))) for s in get_labels(l)]
    K = [([d[string(l)] for l in p], [d[string(l)] for l in n]) for (p,n) in K]
    K = [(p, [HalfSpace(-h.a, -h.b) for h in n]) for (p,n) in K]
    K = [[p;n] for (p,n) in K]
    K = [(stack(c.a for c in h)', [c.b for c in h]) for h in K]

    q0 = [get_init_state_number(l)]
    qT = collect(reduce(∪, get_rabin_acceptance(l)[1]))

    return V, E, K, q0, qT
end
