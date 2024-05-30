using LazySets
using Symbolics
using Spot
import Base: !, &, |

abstract type Formula end
Atom = Union{Formula, LazySets.HalfSpace}
struct Not <: Formula
    f::Atom
end
struct And <: Formula
    f1::Atom
    f2::Atom
end
struct And_ <: Formula
    f::Vector{HalfSpace}
end
struct Or <: Formula
    f1::Atom
    f2::Atom
end
struct G <: Formula
    f::Atom
end
struct F <: Formula
    f::Atom
end
struct U <: Formula
    f1::Atom
    f2::Atom
end
(!)(f::Atom) = Not(f)
(&)(f1::Atom, f2::Atom) = And(f1, f2)
(|)(f1::Atom, f2::Atom) = Or(f1, f2)

struct Conj f::Dict end

function issatisfied(f::Conj, op::Array, on::Array)
    return keys(f.f) ⊆ op
end

_to_lit(x) = "o$(x)"

function _ltl(f::LazySets.HalfSpace)
    k = _to_lit(hash(f))
    return k, Dict(k => Conj(Dict(k => f)))
end
function _ltl(f::And_)
    k = _to_lit(hash(f))
    return k, Dict(k => Conj(Dict(_to_lit(hash(h)) => h for h in f.f)))
end
function _ltl(f::Not)
    s1, d1 = _ltl(f.f)
    return "!($(s1))", d1
end
function _ltl(f::G)
    s1, d1 = _ltl(f.f)
    return "G($(s1))", d1
end
function _ltl(f::F)
    s1, d1 = _ltl(f.f)
    return "F($(s1))", d1
end
function _ltl(f::And)
    s1, d1 = _ltl(f.f1)
    s2, d2 = _ltl(f.f2)
    return "($(s1)) & ($(s2))", merge(d1, d2)
end
function _ltl(f::Or)
    s1, d1 = _ltl(f.f1)
    s2, d2 = _ltl(f.f2)
    return "($(s1)) | ($(s2))", merge(d1, d2)
end
function _ltl(f::U)
    s1, d1 = _ltl(f.f1)
    s2, d2 = _ltl(f.f2)
    return "($(s1)) U ($(s2))", merge(d1, d2)
end

translate(translator::LTLTranslator, f::Formula) = let (f, d) = _ltl(f)
    Spot.translate(translator, SpotFormula(f)), d
end
