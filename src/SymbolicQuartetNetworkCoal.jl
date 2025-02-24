#test written by Sungsik Kong 2025

module SymbolicQuartetNetworkCoal

    using Combinatorics
    using CSV
    using DataFrames
    using Dates
    using PhyloNetworks
    using StaticArrays
    import Random

    const dpoints=7 #decimal points for all parameters when randomly generated
    const eLab="t_"
    const PN = PhyloNetworks

    export
    readTopologyrand,
    network_expectedCF,
    makeEdgeLabel
    #plot_ntwk_with_Symbolic_Names,

    include("symbolicQNC.jl")

end # module
