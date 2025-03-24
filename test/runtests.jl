#test written by Sungsik Kong 2025
using CSV
using PhyloNetworks
using DataFrames
using FileCmp
using SymbolicQuartetNetworkCoal
using Test

@testset "SymbolicQuartetNetworkCoal tests" begin
    include("test_symbolicQNC.jl")
end