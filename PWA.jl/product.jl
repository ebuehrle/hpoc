using Spot
using LinearAlgebra
import Polyhedra
using HybridSystems, MathematicalSystems
using JuMP, HiGHS
using ProgressMeter
include("formula.jl")
include("graph.jl")

ambient_dim(p::HPolytope) = let v = tovrep(p)
    let V = stack(vertices_list(v))'
        rank(V .- V[1,:]')
    end 
end 
fulldim(p::Union{HPolytope, HPolyhedron}) = length(first(p.constraints).a)
zerovolume(p::HPolyhedron) = let n = fulldim(p); let p = intersection(p, Hyperrectangle(zeros(n),1000*ones(n))); ambient_dim(p) <= n - 1 end end
zerosurface(p::HPolyhedron) = let n = fulldim(p); let p = intersection(p, Hyperrectangle(zeros(n),1000*ones(n))); ambient_dim(p) <= n - 2 end end

function partition(d)
    V = []
    println("partitioning")
    @showprogress for (k,v) in d
        if isempty(V)
            V = [
                (([k],[]), HPolyhedron([v])),
                (([],[k]), HPolyhedron([LazySets.HalfSpace(-v.a,-v.b)])),
            ]
            continue
        end

        w = LazySets.HalfSpace(-v.a,-v.b)

        vp = collect((([xk[1];k], xk[2]), intersection(xv, v)) for (xk,xv) in V)
        vpne = map(((x,v),)->!isempty(v) && !zerovolume(v), vp)
        
        vn = collect(((xk[1], [xk[2];k]), intersection(xv, w)) for (xk,xv) in V)
        vnne = map(((x,v),)->!isempty(v) && !zerovolume(v), vn)
        
        V = [vp[vpne]; vn[vnne]]
    end
    return collect(Set(V))
end

function split_edge_disjunctions(E, L)
    split_disjunctions(x) = strip.(split(x, "|"))
    remove_braces(x) = (x[1] == '(') ? x[2:end-1] : x
    split_conjunctions(x) = strip.(split(x, "&"))
    remove_negation(x) = x[2:end]
    remove_1(x) = filter(l -> l != "1", x)
    group_by_signs(x) = (filter(l->!startswith(l,'!'),x), remove_negation.(filter(l->startswith(l,'!'),x)))

    split_label(x) = group_by_signs.(remove_1.(split_conjunctions.(remove_braces.(split_disjunctions(x)))))
    S = split_label.(L)

    T = [[e for _ in 1:length(s)] for (e,s) in zip(E,S)]
    T = vcat(T...)
    S = vcat(S...)

    return T, S
end

function _union_convex(v1, v2; tol=1e-3, M=1e3)
    if v1 ⊆ v2 return v2 end
    if v2 ⊆ v1 return v1 end

    @assert fulldim(v1) == fulldim(v2)
    n = fulldim(v1)

    h = Polyhedra.hrep(Polyhedra.convexhull(
        Polyhedra.polyhedron(v1),
        Polyhedra.polyhedron(v2)
    ))

    m = Model(HiGHS.Optimizer)
    set_silent(m)
    @variable m x[1:n]
    @variable m q1[1:length(v1.constraints)] Bin
    @variable m q2[1:length(v2.constraints)] Bin
    @constraint m h.A * x .<= h.b
    @constraint m [k=1:length(v1.constraints)] v1.constraints[k].a' * x / norm(v1.constraints[k].a) >= v1.constraints[k].b / norm(v1.constraints[k].a) + tol - M*(1-q1[k])
    @constraint m [k=1:length(v2.constraints)] v2.constraints[k].a' * x / norm(v2.constraints[k].a) >= v2.constraints[k].b / norm(v2.constraints[k].a) + tol - M*(1-q2[k])
    @constraint m sum(q1) >= 1
    @constraint m sum(q2) >= 1
    optimize!(m)

    if is_solved_and_feasible(m) return nothing end

    return HPolyhedron(h)
end

function _merge_modes(V, E)
    for (i1,v1) in enumerate(V)
        for (i2,v2) in enumerate(V[1:i1-1])
            if (i1,i2) ∉ E continue end
            if (i2,i1) ∉ E continue end
            if v1[2] != v2[2] continue end
            lout1 = Set([V[e[2]][2] for e in E[eout(stack(E)',i1)]])
            lout2 = Set([V[e[2]][2] for e in E[eout(stack(E)',i2)]])
            linc1 = Set([V[e[1]][2] for e in E[einc(stack(E)',i1)]])
            linc2 = Set([V[e[1]][2] for e in E[einc(stack(E)',i2)]])
            lconn = lout1 ∪ lout2 ∪ linc1 ∪ linc2
            if lconn != Set([v1[2]]) continue end
            v = _union_convex(v1[1], v2[1])
            if !isnothing(v)
                return (i2, i1), (v, v1[2])
            end
        end
    end
    return nothing
end

function merge_eq_modes(V, E, q0, qT)
    @showprogress for _ in 1:length(V)
        r = _merge_modes(V, E)
        if isnothing(r) break end

        (i1, i2), v = r
        Vx = collect(1:length(V))
        V = [V[1:i1-1]; v; V[i1+1:i2-1]; V[i2+1:end]]
        Vx = [Vx[1:i2-1]; i1; Vx[i2+1:end] .- 1]
        ix = x -> Vx[x]
        E = [ix.(e) for e in E]
        E = [e for e in E if e[1] != e[2]]
        q0 = ix.(q0)
        qT = ix.(qT)
    end

    E = collect(Set(E))
    q0 = collect(Set(q0))
    qT = collect(Set(qT))

    return V, E, q0, qT
end

function PPWA(A::Matrix, B::Union{Vector,Matrix}, f::Union{Formula, String}, translator = LTLTranslator(deterministic=true); merge_modes=true, remove_redundant=true, V=nothing, O=nothing)
    
    if isnothing(V)
        l, d = translate(translator, f)
        h = reduce(merge, x.f for x in values(d))
        Vp = partition(h)
        O1 = [o for (o,_) in Vp]
        O1 = map(o -> let k = collect(keys(d)); satisfied = [issatisfied(d[x], o[1], o[2]) for x in k]; (k[satisfied], k[(!).(satisfied)]) end, O1)
        V1 = [v for (_,v) in Vp]
    else
        @assert length(V) == length(O)
        l = Spot.translate(translator, Spot.SpotFormula(f))
        allO = reduce(vcat, O)
        O1 = map(o -> (o, setdiff(allO, o)), O)
        V1 = V
    end

    E1 = [(i,j) for (i,v1) = enumerate(V1) for (j,v2) = enumerate(V1) if (i != j) && !isempty(intersection(v1, v2)) && !zerosurface(intersection(v1, v2))]

    V2 = collect(1:num_states(l))
    E2 = get_edges(l)
    L2 = string.(get_labels(l))
    E2, L2 = split_edge_disjunctions(E2, L2)
    q20 = [get_init_state_number(l),]
    q2T = collect(reduce(∪, get_rabin_acceptance(l)[1]))

    println("constructing product automaton vertices")
    V = [(v1,v2) for (i1,v1) in enumerate(V1) for (i2,v2) in enumerate(V2)]
    O = [O1[i1] for (i1,v1) in enumerate(V1) for (i2,v2) in enumerate(V2)]
    Ix = [(i1,i2) for (i1,v1) in enumerate(V1) for (i2,v2) in enumerate(V2)]
    E = [(i1,i2) for (i1,((i11,i12),v1)) in enumerate(zip(Ix,V)) for (i2,((i21,i22),v2)) in enumerate(zip(Ix,V)) if (i11,i21) ∈ E1 && (i12,i22) ∈ E2]
    
    println("constructing product automaton edges")
    E = filter(((i1,i2),) -> any(
        e == (Ix[i1][2],Ix[i2][2]) &&
        all(Set.(l) .⊆ Set.(O1[Ix[i2][1]])) for (e,l) in zip(E2,L2)), E)
    q0 = [i for (i,(q1,q2)) in enumerate(Ix) if q2 in q20]
    qT = [i for (i,(q1,q2)) in enumerate(Ix) if q2 in q2T]

    if remove_redundant
        println("removing redundant constraints")
        V = @showprogress [(remove_redundant_constraints(k),i) for (k,i) in V]
    end

    if merge_modes
        println("merging equivalent modes")
        V, E, q0, qT = merge_eq_modes(V, E, q0, qT)
    end
    
    K = [k.constraints for (k,_) in V]
    K = [(stack(c.a for c in h)', [c.b for c in h]) for h in K]

    modes = [@system(x' = A*x + B*u, x ∈ HPolytope(a,b), u ∈ Universe(size(B,2))) for (a,b) in K]
    automaton = GraphAutomaton(length(modes))
    [add_transition!(automaton, i, j, k) for (k,(i,j)) in enumerate(E)]
    transitions = [@system(x' = x) for _ in E]
    switchings = [AutonomousSwitching() for _ in E]
    h = HybridSystem(automaton, modes, transitions, switchings)

    return h, q0, qT

end
