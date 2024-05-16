using Spot
using LinearAlgebra
import Polyhedra
using HybridSystems, MathematicalSystems
include("formula.jl")
include("merge.jl")

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

function pwa(A, B, f::Formula, x0, xT, translator::LTLTranslator)
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

    println("removing empty modes")
    e = [!isempty(k) && !zerovolume(k) for (k,_) in V]
    Vx = collect(1:length(V))[e]
    V = V[e]
    O = O[e]
    Ix = Ix[e]
    E = [(findfirst(x -> x == i, Vx), findfirst(x -> x == j, Vx)) for (i,j) in E if e[i] && e[j]]
    
    println("removing redundant modes")
    r = [!any((i == ik) && issubset(k,xk) for (xk,ik) in V[1:j-1]) for (j,(k,i)) in enumerate(V)]
    Vx = collect(1:length(V))[r]
    V = V[r]
    O = O[r]
    Ix = Ix[r]
    E = [(findfirst(x -> x == i, Vx), findfirst(x -> x == j, Vx)) for (i,j) in E if r[i] && r[j]]

    println("merging equivalent modes")
    for _ in 1:length(V)-1
        groups = [[i] for i in 1:length(V)]
        
        function can_merge_modes(i,j)
            if (i,j) ∉ E return end
            if (j,i) ∉ E return end
            p1,i1 = V[i]
            p2,i2 = V[j]
            if i1 != i2 return end
            p = can_merge_polyhedra(p1, p2)
            if !isnothing(p) return (p,i1) end
        end
        groups, Vm = _merge_equivalents(groups, V, can_merge_modes)
        if length(Vm) == length(V) break end
            
        update_mode(i) = first(findall(x -> (i in x), groups))
        Em = collect(Set([update_mode.(e) for e in E]))
        Em = filter(((i,j),) -> (i != j), Em)
        V = Vm
        E = Em
    end

    q0 = [i for (i,(v,k)) in enumerate(V) if x0 in v && k in q20]
    qT = [i for (i,(v,k)) in enumerate(V) if xT in v && k in q2T]

    println("remove unreachable or dead modes")
    T = [(i,j) in E for i in 1:length(V), j in 1:length(V)]
    T = sum(T^i for i in 0:length(V))
    reachable = dropdims(any(T[q0,:] .>= 1, dims=1), dims=1)
    alive = dropdims(any(T[:,qT] .>= 1, dims=2), dims=2)
    keep = reachable .&& alive
    println("pruning $(sum(1 .- keep)) modes ($(sum(1 .- reachable)) unreachable, $(sum(1 .- alive)) dead)")

    # Vx = collect(1:length(V))[keep]
    # V = V[keep]
    # E = [(findfirst(x -> x == i, Vx), findfirst(x -> x == j, Vx)) for (i,j) in E if keep[i] && keep[j]]
    # q0 = [findfirst(x -> x == q, Vx) for q in q0 if keep[q]]
    # qT = [findfirst(x -> x == q, Vx) for q in qT if keep[q]]

    println("removing redundant constraints")
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
