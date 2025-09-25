module SymbolicQuartetNetworkCoal
    using Combinatorics
    using CSV
    using DataFrames
    using Dates
    using PhyloNetworks
    using StaticArrays
    import Random
    using Distributions
    using Logging        
    #global_logger(ConsoleLogger(stderr, Debug))

    const dpoints=3 #decimal points for all parameters when randomly generated
    const eLab="t_"
    const PN = PhyloNetworks

    export
    aloha,
    export_csv,
    readTopologyrand,
    network_expectedCF_formulas,
    makeEdgeLabel,
    assignBinaryEdgeLengths
        
    include("misc.jl")
    include("symbolicQNC.jl")
    include("topology.jl")
    include("reformatting.jl")    

end # module