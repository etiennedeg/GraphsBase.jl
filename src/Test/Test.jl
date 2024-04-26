
"""
    GraphsBase.Test

A module that provides utilities for testing functions that should work with any `GraphsBase.AbstractGraph`.
"""
module Test

using GraphsBase

export GenericEdge, GenericGraph, GenericDiGraph, generic_graph

"""
    GenericEdge <: GraphsBase.AbstractEdge

An edge type that can  be used to tests functions that relay on the GraphsBase.jl interface.

"""
struct GenericEdge{V} <: GraphsBase.AbstractEdge{V, U}
    e::GraphsBase.SimpleEdge{V}
end

GraphsBase.src(e::GenericEdge) = GraphsBase.src(e.e)

GraphsBase.dst(e::GenericEdge) = GraphsBase.dst(e.e)

Base.reverse(e::GenericEdge) = GenericEdge(reverse(e.e))

"""
    GenericGraph{T} <: GraphsBase.AbstractGraph{T}

An undirected graph type that can  be used to tests functions that relay on the GraphsBase.jl interface.

"""
struct GenericGraph{T} <: GraphsBase.AbstractGraph{T}
    g::SimpleGraph{T}
end

"""
    GenericDiGraph{T} <: GraphsBase.AbstractGraph{T}

A directed graph type that can  be used to tests functions that relay on the GraphsBase.jl interface.

"""
struct GenericDiGraph{T} <: GraphsBase.AbstractGraph{T}
    g::SimpleDiGraph{T}
end

"""
    generic_graph(g::Union{SimpleGraph, SimpleDiGraph})

Return either a GenericGraph or GenericDiGraph that wraps a copy of g.
"""
function generic_graph(g::Union{SimpleGraph,SimpleDiGraph})
    g = copy(g)
    return is_directed(g) ? GenericDiGraph(g) : GenericGraph(g)
end

function GenericDiGraph(elist::Vector{GraphsBase.SimpleDiGraphEdge{T}}) where {T<:Integer}
    return GenericDiGraph{T}(SimpleDiGraph(elist))
end

GraphsBase.is_directed(::Type{<:GenericGraph}) = false
GraphsBase.is_directed(::Type{<:GenericDiGraph}) = true

GraphsBase.edges(g::GenericGraph) = (GenericEdge(e) for e in GraphsBase.edges(g.g))
GraphsBase.edges(g::GenericDiGraph) = (GenericEdge(e) for e in GraphsBase.edges(g.g))

GraphsBase.edgetype(g::GenericGraph) = GenericEdge{eltype(g)}
GraphsBase.edgetype(g::GenericDiGraph) = GenericEdge{eltype(g)}

GraphsBase.has_edge(g::GenericGraph, s, d) = GraphsBase.has_edge(g.g, s, d)
GraphsBase.has_edge(g::GenericDiGraph, s, d) = GraphsBase.has_edge(g.g, s, d)

GraphsBase.has_vertex(g::GenericGraph, v) = GraphsBase.has_vertex(g.g, v)
GraphsBase.has_vertex(g::GenericDiGraph, v) = GraphsBase.has_vertex(g.g, v)

GraphsBase.inneighbors(g::GenericGraph, v) = (u for u in GraphsBase.inneighbors(g.g, v))
GraphsBase.inneighbors(g::GenericDiGraph, v) = (u for u in GraphsBase.inneighbors(g.g, v))

GraphsBase.outneighbors(g::GenericGraph, v) = (u for u in GraphsBase.outneighbors(g.g, v))
GraphsBase.outneighbors(g::GenericDiGraph, v) = (u for u in GraphsBase.outneighbors(g.g, v))

GraphsBase.ne(g::GenericGraph) = GraphsBase.ne(g.g)
GraphsBase.ne(g::GenericDiGraph) = GraphsBase.ne(g.g)

GraphsBase.nv(g::GenericGraph) = GraphsBase.nv(g.g)
GraphsBase.nv(g::GenericDiGraph) = GraphsBase.nv(g.g)

GraphsBase.vertices(g::GenericGraph) = (v for v in GraphsBase.vertices(g.g))
GraphsBase.vertices(g::GenericDiGraph) = (v for v in GraphsBase.vertices(g.g))

end # module