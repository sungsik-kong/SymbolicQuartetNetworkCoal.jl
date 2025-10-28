module SymbolicQuartetNetworkCoal
using CSV
using DataFrames
#using Dates
using PhyloNetworks
using StaticArrays
import Random

const dpoints=10 # decimal points for all parameters when randomly generated
const eLab="t_"
const gLab="g_"
const PN = PhyloNetworks

export
aloha,
export_csv,
format_export,
readTopologyrand,
network_expectedCF_formulas,
# makeEdgeLabel_v3,
makeEdgeLabel,
assignBinaryEdgeLengths

include("misc.jl")
include("symbolicCF.jl")
include("formatting.jl")    

end # module
