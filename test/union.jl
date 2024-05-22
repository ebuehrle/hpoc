include("../pwa/product.jl")
using Symbolics, LazySets

Symbolics.@variables x[1:4]
v1 = HPolytope([
    HalfSpace(x[1] <= 0, x),
    HalfSpace(x[2] <=  1, x),
    HalfSpace(x[2] >= -1, x),
])
v2 = HPolytope([
    HalfSpace(x[1] >= 0, x),
    HalfSpace(x[2] <=  1, x),
    HalfSpace(x[2] >= -1, x),
])
_union_convex(v1, v2)
