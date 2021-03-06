module LibXCTests
using LibXC
using Base.Test
using Unitful
using DataFrames: DataFrame
using DFTShims
using Documenter

const DHA = Dispatch.Hartree.Scalars
macro lintpragma(s) end

include("fixtures.jl")
using .Fixtures

@testset "> Internal API" begin
    @test LibXC.LIB_VERSION == v"3.0.0"
    @test LibXC.FUNCTIONALS[:lda_x] == 1
    @test LibXC.FUNCTIONALS[:mgga_x_ms2] == 223
end

@testset "> Array manipulations" begin
    include("ArrayManips.jl")
end

@testset "> Internal Functional API" begin
    include("internals.jl")
end

@testset "> Output tuples" begin
    include("OutputTuples.jl")
end

@testset "> LDA" begin
    include("lda.jl")
end

@testset "> GGA" begin
    include("gga.jl")
end

@testset "Docstests" begin
    mktempdir() do path
        makedocs(modules = [LibXC], clean = true,
             build = joinpath(path, "build"),
             pages = LibXC._doc_pages(),
             root=joinpath(dirname(@__DIR__), "docs"))
    end
end
end
