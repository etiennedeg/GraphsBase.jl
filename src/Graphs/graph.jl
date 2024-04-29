"""
    Graph{V, U}

A type representing an undirected graph.
"""
mutable struct Graph{V} <: AbstractGraph{V, E}
    ne::Int
    fadjlist::Vector{Vector{V}}
    vertex_to_index::Dict{V, Int}

    function Graph{V}(ne::Int, fadjlist::Vector{Vector{V}}, vertex_to_index::Dict{V, Int}) where {V}
        return new{V}(ne, fadjlist, vertex_to_index)
    end
end

function Graph(ne, fadjlist::Vector{Vector{V}}, vertex_to_index::Dict{V, Int}) where {V}
    return Graph{V}(ne, fadjlist, vertex_to_index)
end

Graph(n::V) where {V} = Graph{V}(n)

Graph() = Graph{Int}()

"""
    Graph(::Type{V})

Construct an empty `Graph{V}` with 0 vertices and 0 edges.

## Examples
```jldoctest
julia> using Graphs

julia> Graph(UInt8)
{0, 0} undirected simple UInt8 graph
```
"""
Graph(::Type{V}) where {V} = Graph{V}(0, Vector{V}(), Dict{V, Int}())

"""
    Graph{V}(adjm::AbstractMatrix)

Construct a `Graph{V}` from the adjacency matrix `adjm`.
If `adjm[i][j] != 0`, an edge `(i, j)` is inserted. `adjm` must be a square and symmetric matrix.
The element type `V` can be omitted.

## Examples
```jldoctest
julia> using Graphs

julia> A1 = [false true; true false];

julia> Graph(A1)
{2, 1} undirected simple Int64 graph

julia> A2 = [2 7; 7 0];

julia> Graph{Int16}(A2)
{2, 2} undirected simple Int16 graph
```
"""
Graph(adjmx::AbstractMatrix) = Graph{Int}(adjmx)

function Graph{V}(adjmx::AbstractMatrix) where {V}
    dima, dimb = size(adjmx)
    isequal(dima, dimb) ||
        throw(ArgumentError("Adjacency / distance matrices must be square"))
    issymmetric(adjmx) ||
        throw(ArgumentError("Adjacency / distance matrices must be symmetric"))

    g = Graph(V(dima))
    @inbounds for i in findall(triu(adjmx) .!= 0)
        add_edge!(g, i[1], i[2])
    end
    return g
end

"""
    Graph{V}(g::Graph)

Construct a copy of g.
If the element type `V` is specified, the vertices of `g` are converted to this type.
Otherwise the element type is the same as for `g`.

## Examples
```jldoctest
julia> using Graphs

julia> g = complete_graph(5)
{5, 10} undirected simple Int64 graph

julia> Graph{UInt8}(g)
{5, 10} undirected simple UInt8 graph
```
"""
Graph(g::Graph) = copy(g)


function Graph{V}(g::Graph) where {V}
    h_fadj = [Vector{V}(x) for x in fadj(g)]
    return Graph(ne(g), h_fadj)
end

"""
    Graph(g::SimpleDiGraph)

Construct an undirected `Graph` from a directed `SimpleDiGraph`.
Every directed edge in `g` is added as an undirected edge.
The element type is the same as for `g`.

## Examples
```jldoctest
julia> using Graphs

julia> g = path_digraph(Int8(5))
{5, 4} directed simple Int8 graph

julia> Graph(g)
{5, 4} undirected simple Int8 graph
```
"""
function Graph(g::SimpleDiGraph)
    gnv = nv(g)
    edgect = 0
    newfadj = deepcopy_adjlist(g.fadjlist)
    @inbounds for i in vertices(g)
        for j in badj(g, i)
            index = searchsortedfirst(newfadj[i], j)
            if index <= length(newfadj[i]) && newfadj[i][index] == j
                edgect += 1     # this is an existing edge - we already have it
                if i == j
                    edgect += 1 # need to count self loops
                end
            else
                insert!(newfadj[i], index, j)
                edgect += 2      # this is a new edge only in badjlist
            end
        end
    end
    iseven(edgect) ||
        throw(AssertionError("invalid edgect in graph creation - please file bug report"))
    return Graph(edgect ÷ 2, newfadj)
end

@inbounds function cleanupedges!(fadjlist::Vector{Vector{V}}) where {V}
    neg = 0
    for v in 1:length(fadjlist)
        if !issorted(fadjlist[v])
            sort!(fadjlist[v])
        end
        unique!(fadjlist[v])
        neg += length(fadjlist[v])
        # self-loops should count as one edge
        for w in fadjlist[v]
            if w == v
                neg += 1
                break
            end
        end
    end
    return neg ÷ 2
end

"""
    Graph(edge_list::Vector)

Construct a `Graph` from a vector of edges.
The element type is taken from the edges in `edge_list`.
The number of vertices is the highest that is used in an edge in `edge_list`.

### Implementation Notes
This constructor works the fastest when `edge_list` is sorted
by the lexical ordering and does not contain any duplicates.

### See also
[`GraphFromIterator`](@ref)

## Examples
```jldoctest
julia> using Graphs

julia> el = Edge.([ (1, 2), (1, 5) ])
2-element Vector{Graphs.Graphs.Edge{Int64}}:
 Edge 1 => 2
 Edge 1 => 5

julia> Graph(el)
{5, 2} undirected simple Int64 graph
```
"""
function Graph(edge_list::Vector{Edge{V}}) where {V}
    nvg = zero(V)
    @inbounds(
        for e in edge_list
            nvg = max(nvg, src(e), dst(e))
        end
    )

    list_sizes = ones(Int, nvg)
    degs = zeros(Int, nvg)
    @inbounds(
        for e in edge_list
            s, d = src(e), dst(e)
            (s >= 1 && d >= 1) || continue
            degs[s] += 1
            if s != d
                degs[d] += 1
            end
        end
    )

    fadjlist = Vector{Vector{V}}(undef, nvg)
    @inbounds(
        for v in 1:nvg
            fadjlist[v] = Vector{V}(undef, degs[v])
        end
    )

    @inbounds(
        for e in edge_list
            s, d = src(e), dst(e)
            (s >= 1 && d >= 1) || continue
            fadjlist[s][list_sizes[s]] = d
            list_sizes[s] += 1
            if s != d
                fadjlist[d][list_sizes[d]] = s
                list_sizes[d] += 1
            end
        end
    )

    neg = cleanupedges!(fadjlist)
    g = Graph{V}()
    g.fadjlist = fadjlist
    g.ne = neg

    return g
end

@inbounds function add_to_fadjlist!(
    fadjlist::Vector{Vector{V}}, s::V, d::V
) where {V}
    nvg = length(fadjlist)
    nvg_new = max(nvg, s, d)
    for v in (nvg + 1):nvg_new
        push!(fadjlist, Vector{V}())
    end

    push!(fadjlist[s], d)
    if s != d
        push!(fadjlist[d], s)
    end
end

# Try to get the eltype from the first element
function _GraphFromIterator(iter)::Graph
    next = iterate(iter)
    if (next === nothing)
        return Graph(0)
    end

    e = first(next)
    E = typeof(e)
    if !(E <: Edge{})
        throw(DomainError(iter, "Edges must be of type Edge{V <: Integer}"))
    end

    V = eltype(e)
    g = Graph{V}()
    fadjlist = Vector{Vector{V}}()

    while next != nothing
        (e, state) = next

        if !(e isa E)
            throw(DomainError(iter, "Edges must all have the same type."))
        end
        s, d = src(e), dst(e)
        if ((s >= 1) & (d >= 1))
            add_to_fadjlist!(fadjlist, s, d)
        end

        next = iterate(iter, state)
    end

    neg = cleanupedges!(fadjlist)
    g.fadjlist = fadjlist
    g.ne = neg

    return g
end

function _GraphFromIterator(iter, ::Type{V}) where {V}
    g = Graph{V}()
    fadjlist = Vector{Vector{V}}()

    @inbounds(
        for e in iter
            s, d = src(e), dst(e)
            (s >= 1 && d >= 1) || continue
            add_to_fadjlist!(fadjlist, s, d)
        end
    )

    neg = cleanupedges!(fadjlist)
    g.fadjlist = fadjlist
    g.ne = neg

    return g
end

"""
    GraphFromIterator(iter)

Create a [`Graph`](@ref) from an iterator `iter`. The elements in iter must
be of `type <: Edge`.

# Examples
```jldoctest
julia> using Graphs

julia> g = Graph(3);

julia> add_edge!(g, 1, 2);

julia> add_edge!(g, 2, 3);

julia> h = GraphFromIterator(edges(g));

julia> collect(edges(h))
2-element Vector{Graphs.Graphs.Edge{Int64}}:
 Edge 1 => 2
 Edge 2 => 3
```
"""
function GraphFromIterator(iter)::Graph
    if Base.IteratorEltype(iter) == Base.HasEltype()
        E = eltype(iter)
        if (E <: Edge{} && isconcretetype(E))
            V = eltype(E)
            if isconcretetype(V)
                return _GraphFromIterator(iter, V)
            end
        end
    end

    return _GraphFromIterator(iter)
end

edgetype(::Graph{V}) where {V} = Edge{V}

"""
    badj(g::Graph[, v::Integer])

Return the backwards adjacency list of a graph. If `v` is specified,
return only the adjacency list for that vertex.

### Implementation Notes
Returns a reference to the current graph's internal structures, not a copy.
Do not modify result. If the graph is modified, the behavior is undefined:
the array behind this reference may be modified too, but this is not guaranteed.
"""
badj(g::Graph) = fadj(g)
badj(g::Graph, v::Integer) = fadj(g, v)

"""
    adj(g[, v])

Return the adjacency list of a graph. If `v` is specified, return only the
adjacency list for that vertex.

### Implementation Notes
Returns a reference to the current graph's internal structures, not a copy.
Do not modify result. If the graph is modified, the behavior is undefined:
the array behind this reference may be modified too, but this is not guaranteed.
"""
adj(g::Graph) = fadj(g)
adj(g::Graph, v::Integer) = fadj(g, v)

copy(g::Graph) = Graph(g.ne, deepcopy_adjlist(g.fadjlist))

function ==(g::Graph, h::Graph)
    return vertices(g) == vertices(h) && ne(g) == ne(h) && fadj(g) == fadj(h)
end

"""
    is_directed(g)

Return `true` if `g` is a directed graph.
"""
is_directed(::Type{<:Graph}) = false

function has_edge(g::Graph{V}, s, d) where {V}
    verts = vertices(g)
    (s in verts && d in verts) || return false  # edge out of bounds
    @inbounds list_s = g.fadjlist[s]
    @inbounds list_d = g.fadjlist[d]
    if length(list_s) > length(list_d)
        d = s
        list_s = list_d
    end
    return insorted(d, list_s)
end

function has_edge(g::Graph{V}, e::Edge{V}) where {V}
    s, d = V.(Tuple(e))
    return has_edge(g, s, d)
end

"""
    add_edge!(g, e)

Add an edge `e` to graph `g`. Return `true` if edge was added successfully,
otherwise return `false`.

# Examples
```jldoctest
julia> using Graphs

julia> g = Graph(2);

julia> add_edge!(g, 1, 2)
true

julia> add_edge!(g, 2, 3)
false
```
"""
function add_edge!(g::Graph{V}, e::Edge{V}) where {V}
    s, d = T.(Tuple(e))
    verts = vertices(g)
    (s in verts && d in verts) || return false  # edge out of bounds
    @inbounds list = g.fadjlist[s]
    index = searchsortedfirst(list, d)
    @inbounds (index <= length(list) && list[index] == d) && return false  # edge already in graph
    insert!(list, index, d)

    g.ne += 1
    s == d && return true  # selfloop

    @inbounds list = g.fadjlist[d]
    index = searchsortedfirst(list, s)
    insert!(list, index, s)
    return true  # edge successfully added
end

"""
    rem_edge!(g, e)

Remove an edge `e` from graph `g`. Return `true` if edge was removed successfully,
otherwise return `false`.

### Implementation Notes
If `rem_edge!` returns `false`, the graph may be in an indeterminate state, as
there are multiple points where the function can exit with `false`.

# Examples
```jldoctest
julia> using Graphs

julia> g = Graph(2);

julia> add_edge!(g, 1, 2);

julia> rem_edge!(g, 1, 2)
true

julia> rem_edge!(g, 1, 2)
false
```
"""
function rem_edge!(g::Graph{V}, e::Edge{V}) where {V}
    s, d = T.(Tuple(e))
    verts = vertices(g)
    (s in verts && d in verts) || return false  # edge out of bounds
    @inbounds list = g.fadjlist[s]
    index = searchsortedfirst(list, d)
    @inbounds (index <= length(list) && list[index] == d) || return false  # edge not in graph
    deleteat!(list, index)

    g.ne -= 1
    s == d && return true  # selfloop

    @inbounds list = g.fadjlist[d]
    index = searchsortedfirst(list, s)
    deleteat!(list, index)
    return true  # edge successfully removed
end

fd
"""
    add_vertex!(g)

Add a new vertex to the graph `g`. Return `true` if addition was successful.

# Examples
```jldoctest
julia> using Graphs

julia> g = Graph(Int8(typemax(Int8) - 1))
{126, 0} undirected simple Int8 graph

julia> add_vertex!(g)
true

julia> add_vertex!(g)
false
```
"""
function add_vertex!(g::Graph{V}) where {V}
    (nv(g) + one(V) <= nv(g)) && return false       # test for overflow
    push!(g.fadjlist, Vector{V}())
    return true
end

"""
    rem_vertices!(g, vs, keep_order=false) -> vmap

Remove all vertices in `vs` from `g`.
Return a vector `vmap` that maps the vertices in the modified graph to the ones in
the unmodified graph.
If `keep_order` is `true`, the vertices in the modified graph appear in the same
order as they did in the unmodified graph. This might be slower.

### Implementation Notes
This function is not part of the official Graphs API and is subject to change/removal between major versions.

# Examples
```jldoctest
julia> using Graphs

julia> g = complete_graph(5)
{5, 10} undirected simple Int64 graph

julia> vmap = rem_vertices!(g, [2, 4], keep_order=true);

julia> vmap
3-element Vector{Int64}:
 1
 3
 5

julia> g
{3, 3} undirected simple Int64 graph
```
"""
function rem_vertices!(
    g::Graph{V}, vs::AbstractVector{}; keep_order::Bool=false
) where {V}
    # TODO There might be some room for performance improvements.
    # At the moment, we check for all edges if they stay in the graph.
    # If some vertices keep their position, this might be unnecessary.

    n = nv(g)
    isempty(vs) && return collect(Base.OneTo(n))

    # Sort and filter the vertices that we want to remove
    remove = sort(vs)
    unique!(remove)
    (1 <= remove[1] && remove[end] <= n) ||
        throw(ArgumentError("Vertices to be removed must be in the range 1:nv(g)."))

    # Create a vmap that maps vertices to their new position
    # vertices that get removed are mapped to 0
    vmap = Vector{V}(undef, n)

    if keep_order
        # traverse the vertex list and shift if a vertex gets removed
        i = 1
        @inbounds for u in vertices(g)
            if i <= length(remove) && u == remove[i]
                vmap[u] = 0
                i += 1
            else
                vmap[u] = u - (i - 1)
            end
        end
    else
        # traverse the vertex list and replace vertices that get removed
        # with the furthest one to the back that does not get removed
        i = 1
        j = length(remove)
        v = n
        @inbounds for u in vertices(g)
            u > v && break
            if i <= length(remove) && u == remove[i]
                while v == remove[j] && v > u
                    vmap[v] = 0
                    v -= one(V)
                    j -= 1
                end
                # v > remove[j] || u == v
                vmap[v] = u
                vmap[u] = 0
                v -= one(V)
                i += 1
            else
                vmap[u] = u
            end
        end
    end

    fadjlist = g.fadjlist

    # count the number of edges that will be removed
    # for an edge that gets removed we have to ensure that
    # such an edge does not get counted twice when both endpoints
    # get removed. That's why we relay on the ordering >= on the vertices.
    num_removed_edges = 0
    @inbounds for u in remove
        for v in fadjlist[u]
            if v >= u || vmap[v] != 0
                num_removed_edges += 1
            end
        end
    end
    g.ne -= num_removed_edges

    # move the lists in the adjacency list to their new position
    # The order of traversal is very important here, as otherwise we
    # could overwrite lists, that we want to keep!
    @inbounds for u in (keep_order ? (one(V):1:n) : (n:-1:one(V)))
        if vmap[u] != 0
            fadjlist[vmap[u]] = fadjlist[u]
        end
    end
    resize!(fadjlist, n - length(remove))

    # remove vertices from the lists in fadjlist
    @inbounds for list in fadjlist
        Δ = 0
        for (i, v) in enumerate(list)
            if vmap[v] == 0
                Δ += 1
            else
                list[i - Δ] = vmap[v]
            end
        end
        resize!(list, length(list) - Δ)
        if !keep_order
            sort!(list)
        end
    end

    # we create a reverse vmap, that maps vertices in the result graph
    # to the ones in the original graph. This resembles the output of
    # induced_subgraph
    reverse_vmap = Vector{V}(undef, nv(g))
    @inbounds for (i, u) in enumerate(vmap)
        if u != 0
            reverse_vmap[u] = i
        end
    end

    return reverse_vmap
end
