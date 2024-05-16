using Spot
using LinearAlgebra
import Polyhedra
using HybridSystems, MathematicalSystems
include("formula.jl")

amb_dim(p) = let n = length(first(p.constraints).a);
    let v = tovrep(intersection(p, Hyperrectangle(zeros(n),1000*ones(n))));
        let V = stack(vertices_list(v))';
            rank(V .- V[1,:]')
        end 
    end 
end
fulldim(p) = LazySets.dim(p)
zerovolume(p) = amb_dim(p) <= fulldim(p) - 1
zerosurface(p) = amb_dim(p) <= fulldim(p) - 2

function partition(d)
    V = []
    for (k,v) in d
        if isempty(V)
            V = [
                (([k],[]), HPolyhedron([v])),
                (([],[k]), HPolyhedron([HalfSpace(-v.a,-v.b)])),
            ]
            continue
        end
        V1 = []
        for (xk,xv) in V
            xkp = ([xk[1];k], xk[2])
            xvp = intersection(xv, v)
            xkn = (xk[1], [xk[2];k])
            xvn = intersection(xv, HalfSpace(-v.a,-v.b))
            if !isempty(xvp) && !zerovolume(xvp) push!(V1, (xkp, xvp)) end
            if !isempty(xvn) && !zerovolume(xvn) push!(V1, (xkn, xvn)) end
        end
        V = V1
    end
    return V
end

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

function ppwa(A, B, f::Formula, x0, xT, translator::LTLTranslator)
    l, d = translate(translator, f)

    Vp = partition(d)
    O1 = [o for (o,_) in Vp]
    V1 = [v for (_,v) in Vp]
    E1 = [(i,j) for (i,(_,v1)) = enumerate(Vp) for (j,(_,v2)) = enumerate(Vp) if (i != j) && !isempty(intersection(v1, v2)) && !zerosurface(intersection(v1, v2))]
    q10 = [i for (i,v) in enumerate(V1) if x0 ∈ v]
    q1T = [i for (i,v) in enumerate(V1) if xT ∈ v]

    V2 = collect(1:num_states(l))
    E2 = get_edges(l)
    L2 = string.(get_labels(l))
    E2, L2 = split_edge_disjunctions(E2, L2)
    q20 = [get_init_state_number(l),]
    q2T = collect(reduce(∪, get_rabin_acceptance(l)[1]))

    V = [(v1,v2) for (i1,v1) in enumerate(V1) for (i2,v2) in enumerate(V2)]
    O = [O1[i1] for (i1,v1) in enumerate(V1) for (i2,v2) in enumerate(V2)]
    Ix = [(i1,i2) for (i1,v1) in enumerate(V1) for (i2,v2) in enumerate(V2)]
    E = [(i1,i2) for (i1,((i11,i12),v1)) in enumerate(zip(Ix,V)) for (i2,((i21,i22),v2)) in enumerate(zip(Ix,V)) if (i11,i21) ∈ E1 && (i12,i22) ∈ E2]

    E = filter(((i1,i2),) -> any(
        all(Set.(l) .⊆ Set.(O1[Ix[i2][1]])) && 
        e == (Ix[i1][2],Ix[i2][2]) for (e,l) in zip(E2,L2)), E)
    q0 = [i for (i,(q1,q2)) in enumerate(Ix) if q1 in q10 && q2 in q20]
    qT = [i for (i,(q1,q2)) in enumerate(Ix) if q1 in q1T && q2 in q2T]

    V = [(remove_redundant_constraints(k),i) for (k,i) in V]

    K = [k.constraints for (k,_) in V]
    K = [(stack(c.a for c in h)', [c.b for c in h]) for h in K]

    println("constructing pwa")
    modes = [@system(x' = A*x + B*u, x ∈ HPolytope(a,b), u ∈ Universe(size(B,2))) for (a,b) in K]
    automaton = GraphAutomaton(length(modes))
    [add_transition!(automaton, i, j, k) for (k,(i,j)) in enumerate(E)]
    transitions = [@system(x' = x) for _ in E]
    switchings = [AutonomousSwitching() for _ in E]
    h = HybridSystem(automaton, modes, transitions, switchings)

    return h, q0, qT

end

function fpwa(A, B, l::Formula, x0, xT, translator::LTLTranslator)
    hs, hq0, hqT = ppwa(A, B, l, x0, xT, translator)
    
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
