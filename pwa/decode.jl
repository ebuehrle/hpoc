include("graph.jl")

function decode(E, p, s, t)

    v = s
    P = [s]
    for _ in 1:length(p)

        e = argmax(p .* eout(E, v))
        v = E[e, 2]
        push!(P, v)

        if v == t break end

    end

    @assert P[end] == t

    return P

end
