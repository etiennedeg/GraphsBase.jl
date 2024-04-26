import Base: Pair, Tuple, show, ==, hash, isless
import GraphsBase: AbstractEdge, src, dst

struct Edge{V} <: AbstractEdge{V, Int}
    src::V
    dst::V
end

struct WeightedEdge{V, U} <: AbstractEdge{V, U}
    src::V
    dst::V
    weight::U
end

Edge(t::Tuple) = Edge(t[1], t[2])
WeightedEdge(t::Tuple, w) = WeightedEdge(t[1], t[2], w)
Edge(p::Pair) = Edge(p.first, p.second)
WeightedEdge(p::Pair, w) = WeightedEdge(p.first, p.second, w)
# Edge{V}(p::Pair) where {V<:AbstractVertex} = Edge(V(p.first), V(p.second))
# WeightedEdge{V}(p::Pair, w) where {V<:Vertices} = WeightedEdge(V(p.first), V(p.second), w)
# Edge{V}(t::Tuple) where {V<:AbstractVertex} = Edge(V(t[1]), V(t[2]))
# WeightedEdge{V}(t::Tuple, w) where {V<:AbstractVertex} = WeightedEdge(V(t[1]), V(t[2]), w)

# Accessors
src(e::Edge) = e.src
dst(e::Edge) = e.dst
src(e::WeightedEdge) = e.src
dst(e::WeightedEdge) = e.dst

# I/O
show(io::IO, e::Edge) = print(io, "Edge $(e.src) => $(e.dst)")
show(io::IO, e::WeightedEdge) = print(io, "Edge $(e.src) => $(e.dst) with weight $(e.weight)")

Edge(e::AbstractEdge) where {T<:Integer} = Edge(e.src, e.dst)

# Convenience functions
reverse(e::T) where {T<:AbstractSimpleEdge} = T(dst(e), src(e))
function ==(e1::AbstractSimpleEdge, e2::AbstractSimpleEdge)
    return (src(e1) == src(e2) && dst(e1) == dst(e2))
end
hash(e::AbstractSimpleEdge, h::UInt) = hash(src(e), hash(dst(e), h))
isless(e1::AbstractSimpleEdge, e2::AbstractSimpleEdge) = (src(e1) < src(e2)) || ((src(e1) == src(e2)) && (dst(e1) < dst(e2)))
