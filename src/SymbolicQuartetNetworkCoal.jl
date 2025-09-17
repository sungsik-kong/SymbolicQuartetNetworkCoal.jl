module SymbolicQuartetNetworkCoal
    using Combinatorics
    using CSV
    using DataFrames
    using Dates
    using PhyloNetworks
    using StaticArrays
    import Random
    
    using Logging        
    #global_logger(ConsoleLogger(stderr, Debug))

    const dpoints=7 #decimal points for all parameters when randomly generated
    const eLab="t_"
    const PN = PhyloNetworks

    export
    aloha,
    readTopologyrand,
    network_expectedCF_formulas,
    makeEdgeLabel,
    assignBinaryEdgeLengths
        
    include("misc.jl")
    include("symbolicQNC.jl")
    include("topology.jl")    
end # module