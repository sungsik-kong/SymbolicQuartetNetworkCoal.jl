module SymbolicQuartetNetworkCoal
    using CSV
    using DataFrames
    #using Dates
    using PhyloNetworks
    using StaticArrays
    import Random

    const dpoints=10 #decimal points for all parameters when randomly generated
    const eLab="t_"
    const rLab="r_"
    const PN = PhyloNetworks

    export
    aloha,
    export_csv,
    reformat_export,
    readTopologyrand,
    network_expectedCF_formulas,
    makeEdgeLabel,
    assignBinaryEdgeLengths
        
    include("misc.jl")
    include("symbolicQNC.jl")
    include("reformatting.jl")    

end # module