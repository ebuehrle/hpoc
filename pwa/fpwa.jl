using HybridSystems, MathematicalSystems

struct FPWA
    h::HybridSystem
    fh::HybridSystem
    fq0::Vector{Int}
    fqT::Vector{Int}
    V::Vector{Tuple}
end

function FPWA(hs, hq0, hqT)
    A = HybridSystems.mode(hs,1).A
    B = HybridSystems.mode(hs,1).B

    VV = [(source(hs,t),target(hs,t)) for t in HybridSystems.transitions(hs)]
    V = Tuple.(Set([Set(m) for m in VV]))
    M = [remove_redundant_constraints(HPolyhedron([
        HybridSystems.mode(hs,i).X.constraints; 
        HybridSystems.mode(hs,j).X.constraints])) for (i,j) in V]
    Ex = [((i1,Set(v1)),(i2,Set(v2))) for (i1,v1) in enumerate(V) for (i2,v2) in enumerate(V) if (i1 != i2) && !isempty(Set(v1) ∩ Set(v2))]
    E = [(i1,i2) for ((i1,_),(i2,_)) in Ex]
    K = [HybridSystems.mode(hs, first(Set(V[i]) ∩ Set(V[j]))).X for (i,j) in E]
    q0 = [i for (i,(j,k)) in enumerate(V) if ((j in hq0) || (k in hq0))]
    qT = [i for (i,(j,k)) in enumerate(V) if ((j in hqT) || (k in hqT))]

    automaton = GraphAutomaton(length(M))
    [add_transition!(automaton, i, j, k) for (k,(i,j)) in enumerate(E)]
    modes = [@system(x' = 0, x ∈ m) for m in M]
    transitions = [@system(x' = A*x + B*u, x ∈ k, u ∈ Universe(size(B,2))) for k in K]
    switchings = [AutonomousSwitching() for _ in transitions]
    h = HybridSystem(automaton, modes, transitions, switchings)

    return FPWA(hs, h, q0, qT, V)
end
