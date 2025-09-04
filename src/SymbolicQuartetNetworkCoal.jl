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
    aloha,
    readTopologyrand,
    network_expectedCF_formulas,
    makeEdgeLabel,
    assignBinaryEdgeLengths,
    binaryValToString

    include("misc.jl")
    include("topology.jl")
    include("symbolicQNC.jl")

end # module