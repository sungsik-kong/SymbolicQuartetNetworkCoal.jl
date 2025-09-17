"""
    readTopologyrand(net; scaleparameter=1.0, dpoints=7)

Generates randomized values for all edge lengths and inheritance probabilities in the input network.

## Description
The input network can be provided as either a Newick/Extended Newick string or a PhyloNetworks `HybridNetwork` object, 
with or without parameters specified. The function assumes a **binary network** and produces a **non-ultrametric network**.

- Edge lengths are determined using:  
  `scaleparameter * (edge.number + random value in [0,1])`  
  The `scaleparameter` defaults to `1.0`.
- Inheritance probabilities (gamma values) are assigned random values in the range (0,1).

**Caution:**  
This function does **not** simulate realistic parameters in an empty network topology. 
This function is not intended to simulate parameters with biological background.
Instead, it assigns **dummy** values to allow the topology to be processed by downstream functions.

## Arguments
- `net`: The input network (Newick string, file containing Newick, or a `HybridNetwork` object).
- `scaleparameter`: Multiplier applied to the generated edge lengths. Defaults to `1.0`.
- `dpoints`: Number of decimal places for rounding. Defaults to `7`.

## Returns
- A PhyloNetworks `HybridNetwork` object with randomly assigned edge lengths and inheritance probabilities.
"""
function readTopologyrand(net;scaleparameter::Float64=1.0,dpoints::Integer=dpoints)
    #---------read in topology: input is either newick string or HybridNetwork object---------#
    if net isa PhyloNetworks.HybridNetwork
    else net=PhyloNetworks.readTopology(net) end

    #--------generate arbitrary edge lengths--------#
    #for e in net.edge e.length=round((scaleparameter*(e.number+rand())),digits=dpoints) end
    for e in net.edge e.length=round((scaleparameter*(e.number+0.1)),digits=dpoints) end

    #--------generaete arbitrary inheritance probabilities--------#
    #----preambles----#
    reticulatenodeindex=Int[]
    nreticulate=net.numHybrids
    gammavec=zeros(nreticulate)
    #getting hybrid node index numbers
    for n in net.node n.hybrid && push!(reticulatenodeindex,n.number) end
    #check the number of hybrid nodes are counted correctly
    length(reticulatenodeindex)==nreticulate || @error "Inheritance probability generation failed. Retry."
    #generate arbitrary gamma values n=number of reticulation nodes (that will be assigned to one of the incoming edges)
    gammavec .= round.(rand(nreticulate), digits=dpoints)
    #assign inheritance probabilities to all reticulate edges
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

    return net
end