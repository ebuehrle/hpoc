function can_merge_polyhedra(p1, p2)
    p1r = remove_redundant_constraints(p1)
    p2r = remove_redundant_constraints(p2)
    if p1r == p2r return p1r end

    for i = 1:length(p1.constraints)
        p1f = HPolyhedron([
            p1.constraints[1:i-1];
            HalfSpace(-p1.constraints[i].a, -p1.constraints[i].b);
            p1.constraints[i+1:end]
        ])
        p1m = HPolyhedron([
            p1.constraints[1:i-1];
            p1.constraints[i+1:end]
        ])
        p1r = remove_redundant_constraints(p1f)
        p2r = remove_redundant_constraints(p2)
        if p1r == p2r return remove_redundant_constraints(p1m) end
    end
    
    for i = 1:length(p2.constraints)
        p2f = HPolyhedron([
            p2.constraints[1:i-1];
            HalfSpace(-p2.constraints[i].a, -p2.constraints[i].b);
            p2.constraints[i+1:end]
        ])
        p2m = HPolyhedron([
            p2.constraints[1:i-1];
            p2.constraints[i+1:end]
        ])
        p1r = remove_redundant_constraints(p1)
        p2r = remove_redundant_constraints(p2f)
        if p1r == p2r return remove_redundant_constraints(p2m) end
    end
end

function _merge_equivalents(L, K, can_merge)
    for (j,(lj,kj)) in enumerate(zip(L,K))
        for (i,(li,ki)) in enumerate(zip(L[1:j-1],K[1:j-1]))
            p = can_merge(i,j)
            if !isnothing(p)
                L1 = [[[li;lj]]; L[1:i-1]; L[i+1:j-1]; L[j+1:end]]
                K1 = [p; K[1:i-1]; K[i+1:j-1]; K[j+1:end]]
                return L1, K1
            end
        end
    end
    return L, K
end

function merge_equivalents(K, equivalent)
    L = [[i] for i in 1:length(K)]
    for _ in 1:length(L)-1
        n = length(L)
        L, K = _merge_equivalents(L, K, equivalent)
        if length(L) == n break end
    end
    return L, K
end
