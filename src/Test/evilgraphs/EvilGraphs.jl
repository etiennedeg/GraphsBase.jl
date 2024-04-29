module EvilGraphs

using GraphsBase
using SimpleGraphs
using SimpleTraits
using Random


import Base:
    eltype, show, ==, Pair, Tuple, copy, length, issubset, reverse, zero, in, iterate

# import GraphsBase:
#     _NI, AbstractGraph, AbstractEdge, AbstractEdgeIter,
#     src, dst, edgetype, nv, ne, vertices, edges, outedges, inedges, is_directed,
# 	is_simply_mutable, is_range_based,
#     has_vertex, has_edge, inneighbors, outneighbors, all_neighbors,
# 	get_vertex_container, get_edge_container,
#     deepcopy_adjlist, indegree, outdegree, degree, has_self_loops,
#     insorted

# export AbstractSimpleGraph, AbstractSimpleEdge,
#     SimpleEdge, SimpleGraph, SimpleGraphFromIterator, SimpleGraphEdge,
#     SimpleDiGraph, SimpleDiGraphFromIterator, SimpleDiGraphEdge,
#     add_vertex!, add_edge!, rem_vertex!, rem_vertices!, rem_edge!

rng_gen(g::AbstractGraph) = Xoshiro(hash(g))

"""
    EvilEdge <: AbstractEdge

An adversarial edge type that can be used in to tests functions that relay on the jl interface.

"""
struct EvilEdge{T} <: AbstractEdge{T, Int}
    e::SimpleEdge{T}
end

src(e::EvilEdge) = src(e.e)

dst(e::EvilEdge) = dst(e.e)

Base.reverse(e::EvilEdge) = EvilEdge(reverse(e.e))


"""
    EvilGraph

An adversarial graph structure that can be used in to tests functions that relay on the jl interface.
"""
abstract type AbstractSimpleGraph{T<:Integer} <: AbstractGraph{T, AbstractSimpleEdge{T}} end

struct EvilGraph{T} <: AbstractGraph{T}
    g::SimpleGraph{T}
end

function EvilGraph(elist::Vector{SimpleGraphEdge{T}}) where {T<:Integer}
    return EvilGraph{T}(SimpleGraph(elist))
end

"""
    EvilDiGraph{T} <: AbstractGraph{T}

A directed graph type that can  be used to tests functions that relay on the jl interface.

"""
struct EvilDiGraph{T} <: AbstractGraph{T}
    g::SimpleDiGraph{T}
end

"""
    evil_graph(g::Union{SimpleGraph, SimpleDiGraph})

Return either a EvilGraph or EvilDiGraph that wraps a copy of g.
"""
function evil_graph(g::Union{SimpleGraph,SimpleDiGraph})
    g = copy(g)
    return is_directed(g) ? EvilDiGraph(g) : EvilGraph(g)
end

function EvilDiGraph(elist::Vector{SimpleDiGraphEdge{T}}) where {T<:Integer}
    return EvilDiGraph{T}(SimpleDiGraph(elist))
end

is_directed(::Type{<:EvilGraph}) = false
is_directed(::Type{<:EvilDiGraph}) = true

edges(g::EvilGraph) = (EvilEdge(e) for e in shuffle(rng_gen(g), edges(g.g)))
edges(g::EvilDiGraph) = (EvilEdge(e) for e in shuffle(rng_gen(g), edges(g.g)))

edgetype(g::EvilGraph) = EvilEdge{eltype(g)}
edgetype(g::EvilDiGraph) = EvilEdge{eltype(g)}

has_edge(g::EvilGraph, s, d) = has_edge(g.g, s, d)
has_edge(g::EvilDiGraph, s, d) = has_edge(g.g, s, d)

has_vertex(g::EvilGraph, v) = has_vertex(g.g, v)
has_vertex(g::EvilDiGraph, v) = has_vertex(g.g, v)

inneighbors(g::EvilGraph, v) = (u for u in shuffle(rng_gen(g), inneighbors(g.g, v)))
inneighbors(g::EvilDiGraph, v) = (u for u in shuffle(rng_gen(g), inneighbors(g.g, v)))

outneighbors(g::EvilGraph, v) = (u for u in shuffle(rng_gen(g), outneighbors(g.g, v)))
outneighbors(g::EvilDiGraph, v) = (u for u in shuffle(rng_gen(g), outneighbors(g.g, v)))

ne(g::EvilGraph) = ne(g.g)
ne(g::EvilDiGraph) = ne(g.g)

nv(g::EvilGraph) = nv(g.g)
nv(g::EvilDiGraph) = nv(g.g)

vertices(g::EvilGraph) = (v for v in vertices(g.g))
vertices(g::EvilDiGraph) = (v for v in vertices(g.g))

end