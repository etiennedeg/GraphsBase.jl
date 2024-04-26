using Aqua
using Documenter
using GraphsBase
using JET
using JuliaFormatter
using Test
# using Graphs

DocMeta.setdocmeta!(GraphsBase, :DocTestSetup, :(using GraphsBase); recursive=true)

# using Graphs.SimpleGraphs
# using Graphs.Experimental
# using Graphs.Test
# using SparseArrays
# using LinearAlgebra
# using Compat
# using DelimitedFiles
# using Base64
# using Random
# using Statistics: mean, std
# using StableRNGs
# using Pkg

const testdir = dirname(@__FILE__)

function testgraphs(g)
    return if is_directed(g)
        [g, DiGraph{UInt8}(g), DiGraph{Int16}(g)]
    else
        [g, Graph{UInt8}(g), Graph{Int16}(g)]
    end
end
testgraphs(gs...) = vcat((testgraphs(g) for g in gs)...)
testdigraphs = testgraphs

# function test_generic_graphs(g; eltypes=[UInt8, Int16], skip_if_too_large::Bool=false)
#     SG = is_directed(g) ? SimpleDiGraph : SimpleGraph
#     GG = is_directed(g) ? GenericDiGraph : GenericGraph
#     result = GG[]
#     for T in eltypes
#         if skip_if_too_large && nv(g) > typemax(T)
#             continue
#         end
#         push!(result, GG(SG{T}(g)))
#     end
#     return result
# end

# function test_generic_graphs(gs...; kwargs...)
#     return vcat((test_generic_graphs(g; kwargs...) for g in gs)...)
# end

# function test_large_generic_graphs(g; skip_if_too_large::Bool=false)
#     return test_generic_graphs(
#         g; eltypes=[UInt16, Int32], skip_if_too_large=skip_if_too_large
#     )
# end

tests = [
    "simplegraphs/runtests",
    "interface",
    "core",
    "utils",
]


@testset verbose = true "GraphsBase.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(GraphsBase)
    end
    @testset "Formatting (JuliaFormatter.jl)" begin
        @test format(GraphsBase; verbose=false, overwrite=false)
    end
    # @testset "Doctests (Documenter.jl)" begin
    #     doctest(GraphsBase)
    # end
    @testset "Static checking (JET.jl)" begin
        if VERSION >= v"1.9"
            JET.test_package(GraphsBase; target_modules=(GraphsBase,))
        end
    end
    
    @testset verbose = true "Actual tests" begin
        for t in tests
            tp = joinpath(testdir, "$(t).jl")
            include(tp)
        end
    end
end


