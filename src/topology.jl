"""
    readTopologyrand(net; scaleparameter::Float64=1.0, decimalpoints::Integer=dpoints)

Assigns randomized edge lengths and inheritance probabilities to a phylogenetic network.

# Description
The input can be provided either as:
- a `PhyloNetworks.HybridNetwork` object, or
- a Newick/Extended Newick string (which will be imported into a `PhyloNetworks.HybridNetwork`).

The function assumes a **binary network** and produces a **non-ultrametric network** with dummy parameters
that are suitable for testing or downstream processing, but not meant to represent realistic biological values.

## Edge lengths
- Each edge length is assigned as: `round(scaleparameter * rand(), digits=decimalpoints)`
- By default, all edge lengths < scaleparameter.

## Inheritance probabilities (gamma values)
- For each hybrid node, two incoming edges are expected.
- One edge is assigned a gamma value sampled uniformly from: `Uniform(1/nedge, 0.499)` where nedge is the number of edges in the network.
- The other incoming edge is assigned the complementary probability: `1 - gamma`. Values are rounded to the specified number of decimalpoints.

## Arguments
- `net` A PhyloNetworks.HybridNetwork object or a Newick/Extended Newick string.
- `scaleparameter::Float64=1.0` Multiplier applied to generated edge lengths.
- `decimalpoints::Integer=dpoints` Number of decimal places for rounding (defaults to global dpoints).

## Returns
- A `PhyloNetworks.HybridNetwork` object with randomized edge lengths and inheritance probabilities.

## Notes
- This function will throw an error if any hybrid node does not have exactly two incoming edges.
- The parameters generated are arbitrary placeholders and should not be interpreted as biologically realistic.
"""
function readTopologyrand(net;scaleparameter::Float64=1.0,decimalpoints::Integer=dpoints)
    #---------read in topology: input is either newick string or HybridNetwork object---------#
    if net isa PhyloNetworks.HybridNetwork
    else net=PhyloNetworks.readTopology(net) end

    #--------generate arbitrary edge lengths--------#
    for e in net.edge e.length=round(scaleparameter*(rand()),digits=decimalpoints) end #edge lengths <1, may be more realistic values
    #for e in net.edge e.length=round((scaleparameter*(e.number+0.01)),digits=dpoints) end #this option makes it easier to keep track of elengths during debugging

    #--------generaete arbitrary inheritance probabilities--------#
    #----preambles----#
    reticulatenodeindex=Int[]
    nedge=length(net.edge)
    nreticulate=net.numHybrids
    gammavec=zeros(nreticulate)
    #getting hybrid node index numbers
    for n in net.node n.hybrid && push!(reticulatenodeindex,n.number) end
    #check the number of hybrid nodes are counted correctly
    length(reticulatenodeindex)==nreticulate || @error "Inheritance probability generation failed. Retry."
    #generate arbitrary gamma values n=number of reticulation nodes (that will be assigned to one of the incoming edges)
    gammavec .= round.(rand(Uniform(parse(Float64,"0.$nedge"),0.499),nreticulate), digits=decimalpoints)
    #assign inheritance probabilities to all reticulate edges
    for (nthgamma, nodeidx) in enumerate(reticulatenodeindex)
        # collect the incoming hybrid edges to this node
        incoming = [e for e in net.edge if e.hybrid && PhyloNetworks.getchild(e).number == nodeidx]

        if length(incoming) != 2
            error("Hybrid node $nodeidx has $(length(incoming)) incoming edges (expected 2).")
        end

        incoming[1].gamma = gammavec[nthgamma]
        incoming[2].gamma = round(1 - gammavec[nthgamma], digits=decimalpoints)
    end
    
    #=
    for nthgamma in 1:nreticulate
        visits = 0
        for e in net.edge
            if e.hybrid && PhyloNetworks.getchild(e).number == reticulatenodeindex[nthgamma]
                visits += 1
                if visits == 1
                    e.gamma = gammavec[nthgamma]
                elseif visits == 2
                    e.gamma = round(1 - gammavec[nthgamma], digits=dpoints)
                else
                    error("Hybrid node $(reticulatenodeindex[nthgamma]) has more than 2 incoming edges.")
                end
            end
        end
    end
    =#

    return net
end