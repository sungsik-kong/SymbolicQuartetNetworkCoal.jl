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
    for e in net.edge e.length=round(scaleparameter*(1+rand()),digits=decimalpoints) end #edge lengths <1, may be more realistic values
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
    gammavec .= round.(rand(), digits=decimalpoints)
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
    return net
end

"""
    aloha()

Prints an ASCII art representation along with the text `Hawai'i-Five-O`.
If you see the Hawaiian "shaka" hand gesture, relax and take it easy, 
because your `SymbolicQuartetNetworkCoal.jl` is installed correctly.
"""
function aloha(;scale::Int=1)
    ascii = raw"""
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣠⣤⣴⠂⢀⡀⠀⢀⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⣄⣠⣀⠀⠘⠋⠉⠉⠁⠀⠺⣿⡷⣿⣿⣿⡿⠀⢀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣿⣿⣛⠛⠉⠀⠀⠀⠀⠺⣷⣦⠀⠀⠀⠙⠛⠉⠀⠀⠈⣿⣦⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⢀⣴⣿⣆⠀⠈⠉⠁⠀⠀⠀⠀⠀⠀⠀⠙⠉⠀⠀⢸⣦⠀⠀⠀⢀⣼⣿⣿⣿⣿⣷⡄⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⢺⣿⣿⡿⠁⠀⠀⠀ ⠀⠀⠀⠀⠀⠀  ⠀⠀⠀⠀⠀⠀⠀⠀⢻⣿⣿⣿⣿⣿⣿⣧⡀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⢀⡆⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣤⣤⣀⠀⠀⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣿⣿⣿⣷⠄⠀⠀⠀⠀⠀
⢠⣾⣷⣦⡀⠘⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣿⣿⣿⣿⣿⣷⣦⡀⠀⠀⠀⠀⢠⣿⣿⣿⡿⠟⠛⠋⠁⣀⣠⣤⣄⣀⠀
⠘⣿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⢠⣴⣾⣶⣿⣿⣿⣿⣿⣿⣟⠘⣿⣷⡀⠀⠀⠘⠿⡿⠉⠀⠀⣀⣴⣾⣿⣿⣿⣿⡿⡂
⠀⠈⠿⢟⣿⣿⣆⠀⠀⠀⠀⢀⣤⣤⣿⣿⣿⣿⣿⣎⠛⢫⣿⣿⣿⣷⡘⢿⣿⣆⠀⠀⠀⠀⠀⢀⣾⣿⣿⣿⣿⣿⡿⠛⠋⠁
⠀⠀⠀⢺⣿⣿⣿⣷⡄⠀⢰⣿⣿⣿⣯⠹⣿⣿⣿⣷⣶⡜⢿⣿⣿⣿⣷⡄⠹⣿⣷⣄⠀⠀⣴⣼⣿⣿⣿⣿⠟⠉⠀⠀⠀⠀
⠀⠀⠀⠀⠙⢻⣿⣿⣿⣦⡘⢿⣿⣿⣿⣃⠉⢿⣿⣿⣿⣿⡌⢻⣝⠻⠿⢃⡀⣿⣿⣿⣷⣶⣿⣿⣿⣿⡿⠃⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢻⣿⣿⣿⣷⣄⠛⣿⣿⣿⣿⣄⠻⣿⣿⣿⣿⡆⠙⠷⠶⠟⢠⣿⣿⣿⣿⣿⣿⣿⣿⡟⠁⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⡇⠈⠻⣿⣿⣿⣷⠘⠧⣉⣁⡴⠀⢠⣤⣶⣿⣿⣿⣿⣿⣿⣿⣿⡟⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⣿⣿⣿⣷⣄⠙⠧⣍⣩⡜⢀⣀⣀⠄⣴⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢿⣿⣿⣿⣿⣷⣦⣄⣀⣤⣾⣿⣿⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠛⠛⠻⠛⠛⠁⠉⠙⠛⠉⠉⠉⠀⠀⠀⠀ 
               Hawai'i-Five-O"""

    function scale_ascii_art(art::String; scale::Int=scale)
        lines = split(art, '\n')
        new_lines = String[]
        for (i, line) in enumerate(lines)
            if (i - 1) % scale == 0   
                reduced = String([c for (j, c) in enumerate(line) if (j - 1) % scale == 0]) 
                push!(new_lines, reduced)
            end
        end
        return join(new_lines, "\n")
    end

    println(scale_ascii_art(ascii, scale=scale))
end

"""
    binary_to_tstring(binary_str::String) -> String

Convert a binary string (e.g., "10101") into a string of positions prefixed by the global constant eLab (`t_` by default) and joined by `-`.
Each `1` in the binary string indicates a position, counted from the rightmost bit as position 1.

## Arguments:
- binary_str::String : A string containing only '0' and '1'.

## Returns:
- String : A string representing positions of '1's in the form "t_1-t_3-...".

## Examples:
- binary_to_tstring("10101") # returns "t_5-t_3-t_1"
- binary_to_tstring("0100") # returns "t_3"
"""
function binary_to_tstring(binary_str::String; edge_label=eLab::String)
    # Check that it's binary
    if any(c -> c != '0' && c != '1', binary_str)
        error("Input must represent a binary number (e.g., 10101101).")
    end    
    # Find positions with '1'
    t = length(binary_str)
    positions = [t - i + 1 for (i, c) in enumerate(binary_str) if c !== '0']
    # Join as required
    return join(["$edge_label{$i}" for i in positions], "-")
end
   
"""
    gettingSymbolicInput(net::HybridNetwork, df, inheritancecorrelation)

Transforms numerical values in the concordance factor (CF) equations into symbolic representations.

## Description
This function updates a dataframe (`df`) containing CF equations by replacing:
- The inheritance correlation `"rho"` with its numerical value.
- Edge-related terms `"-t_{e}"` with `"X{e}"`.
- Hybrid-related terms `"r_{e}"` with `"R{e}"`.
- Exponential terms `"exp(-t_{e})"` into `"(X{e}"`.
- Other symbolic replacements for consistency.

Additionally, it extracts and returns a **sorted list of unique symbolic parameters** used in the CF equations.

## Arguments
- `net`: A `HybridNetwork` object.
- `df`: A dataframe containing CF equations.
- `inheritancecorrelation`: The numerical value of the inheritance correlation.

## Returns
- A sorted vector of unique symbolic parameters used in the CF equations.
"""   
function gettingSymbolicInput(net::HybridNetwork, df, inheritancecorrelation)
    edgenumber = length(net.edge)
    retnumber = length(net.hybrid)
    numCFs = size(df, 1)
    
    params = String[]

    for i in 1:numCFs
        df[i, 2] = replace(df[i, 2], "rho" => "$inheritancecorrelation")
        record = false
        expressions = String[]

        for e in 1:edgenumber
            if occursin("-t_{$e}", df[i, 2])
                push!(expressions, "X$e")
            end
        end

        for e in 1:retnumber
            if occursin("r_{$e}", df[i, 2])
                push!(expressions, "R$e")
            end
        end    

        append!(params, expressions)
    end

    params = unique(params)

    # Replace symbolic expressions in CF equations
    for cf in 1:numCFs
        df[cf, 2] = replace(df[cf, 2], "r_{" => "R")   # Replace hybrid parameter notation
        df[cf, 2] = replace(df[cf, 2], "exp(-t_{" => "(X")  # Replace exponential notation
        df[cf, 2] = replace(df[cf, 2], "-t_{" => "*X")  # Replace edge length notation
        df[cf, 2] = replace(df[cf, 2], "})" => ")")     # Close parentheses
        df[cf, 2] = replace(df[cf, 2], "}" => "")       # Remove extra braces
    end

    return sort!(params)
end
    
"""
    makeEdgeLabel(net; showTerminalEdgeLabels=false)

Generates a dataframe mapping edge numbers to their symbolic labels.

## Description
This function creates labels for the edges of a `HybridNetwork` in the format `"t_{e}"`, where `e` is the edge number.  
By default, labels are only assigned to **non-terminal edges** (i.e., edges that do not end at leaf nodes).  
The output dataframe is used as input for PhyloPlots' option `edgelabel=``.
Setting `showTerminalEdgeLabels=true` includes labels for terminal edges as well.

## Arguments
- `net`: A `HybridNetwork` object.
- `showTerminalEdgeLabels`: A boolean flag (default = `false`).  
   - `false`: Excludes terminal edges.  
   - `true`: Includes all edges.  

## Returns
- A `DataFrame` with columns:
  - `number`: Edge numbers.
  - `label`: Corresponding symbolic labels (`"t_{e}"`).
"""
function makeEdgeLabel(net; showTerminalEdgeLabels=false)
    edges_to_include = [e for e in net.edge if !PhyloNetworks.getchild(e).leaf || showTerminalEdgeLabels]
    df = DataFrame(
        number = [e.number for e in edges_to_include],
        label = ["t_{$(e.number)}" for e in edges_to_include]
    )
    return df
end

function sqnc_removedegree2nodes!(net::HybridNetwork, dict, keeproot::Bool=false)
    rootnode = getroot(net)
    # caution: the root and its incident edges may change when degree-2 nodes
    #          are removed. Indices of nodes to be removed would change too.
    rootin2cycle(nn) = isrootof(nn,net) && all(e.hybrid for e in nn.edge)
    toberemoved(nn) = (keeproot ? length(nn.edge) == 2 && nn !== rootnode :
                                  length(nn.edge) == 2 && !rootin2cycle(nn))
    ndegree2nodes = sum(toberemoved.(net.node))
    for _ in 1:ndegree2nodes # empty if 0 degree-2 nodes
        i = findfirst(toberemoved, net.node)
        # i may be 'nothing' if the root was initially thought to be removed
        # but later its edges turned to be hybrids, so should not be removed
        isnothing(i) || Qfuseedgesat!(i, net, dict)
    end
    return net,dict
end





function Qdeleteleaf!(net::HybridNetwork, node::PhyloNetworks.Node, synth_e_dict; kwargs...)
    node.leaf || error("node number $(node.number) is not a leaf.")
    Qdeleteleaf!(net, node.number,synth_e_dict; kwargs..., index=false)
end

function Qdeleteleaf!(net::HybridNetwork, nodeName::AbstractString, synth_e_dict; kwargs...)
    tmp = findall(n -> n.name == nodeName, net.node)
    if length(tmp)==0
        error("node named $nodeName was not found in the network.")
    elseif length(tmp)>1
        error("several nodes were found with name $nodeName.")
    end
    Qdeleteleaf!(net, tmp[1],synth_e_dict; kwargs..., index=true)
end

# recursive algorithm. nodes previously removed are all necessaily
# *younger* than the current node to remove. Stated otherwise:
# edges previously removed all go "down" in time towards current node:
# - tree edge down to an original leaf,
# - 2 hybrid edges down to a hybrid node.
# hybrid edges from node to another node are not removed. fused instead.
# consequence: node having 2 hybrid edges away from node should not occur.
function Qdeleteleaf!(net::HybridNetwork, nodeNumber::Integer, synth_e_dict;
                     index=false::Bool, nofuse=false::Bool,
                     simplify=true::Bool, unroot=false::Bool,
                     multgammas=false::Bool, keeporiginalroot=false::Bool)

    i = nodeNumber # good if index=true
    if !index
        i = findfirst(n -> n.number == nodeNumber, net.node)
        i !== nothing ||
            error("cannot delete leaf number $(nodeNumber) because it is not part of net")
    elseif i > length(net.node)
        error("node index $i too large: the network only has $(length(net.node)) nodes.")
    end
    nodei = net.node[i]
    nodeidegree = length(nodei.edge)
    if nodeidegree == 0
        length(net.node)==1 || error("leaf $(nodei.name) has no edge but network has $(length(net.node)) nodes (instead of 1).")
        @warn "Only 1 node. Removing it: the network will be empty"
        PhyloNetworks.deleteNode!(net,nodei) # empties the network
    elseif nodeidegree == 1
        pe = nodei.edge[1]
        pn = PhyloNetworks.getOtherNode(pe, nodei) # parent node of leaf
        if net.root == i && keeporiginalroot
            return nothing
        end
        # keep nodei if pn is a leaf: keep 1 edge for the single remaining leaf
        if pn.leaf
            net.root = i # it should have been i before anyway
            length(net.edge)==1 || error("neighbor of degree-1 node $(nodei.name) is a leaf, but network had $(length(net.edge)) edges (instead of 1).")
            length(pn.edge)==1 || error("neighbor of $(nodei.name) is a leaf, incident to $(length(pn.edge)) edges (instead of 1)")
            return nothing
        end
        # remove nodei and pe.
        PhyloNetworks.removeNode!(pn,pe)  # perhaps useless. in case gc() on pe affects pn
        PhyloNetworks.removeEdge!(pn,pe)
        PhyloNetworks.deleteEdge!(net,pe,part=false)
        if net.root==i # if node was the root, new root = pn
            net.root = findfirst(x -> x===pn, net.node)
        end
        PhyloNetworks.deleteNode!(net,nodei) # this updates the index net.root
        Qdeleteleaf!(net, pn.number, synth_e_dict; nofuse = nofuse, simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
        return nothing
    elseif nodeidegree > 2
        # do nothing: nodei has degree 3+ (through recursive calling)
        return nothing
    end
    # if we get to here, nodei has degree 2 exactly: --e1-- nodei --e2--
    if i==net.root && (keeporiginalroot || !unroot)
        return nothing # node = root of degree 2 and we want to keep it
    end
    e1 = nodei.edge[1]
    e2 = nodei.edge[2]
    if e1.hybrid && e2.hybrid
        cn  = PhyloNetworks.getchild(e1)
        cn2 = PhyloNetworks.getchild(e2)
        if !(nodei ≡ cn && nodei ≡ cn2) # nodei *not* the child of both e1 and e2
            # possible at the root, in which case e1,e2 should have same child
            (i==net.root && cn ≡ cn2) ||
                error("after removing descendants, node $(nodei.number) has 2 hybrid edges but is not the child of both.")
            # delete e1,e2,nodei and move the root to their child cn
            cn.hybrid || error("child node $(cn.number) of hybrid edges $(e1.number) and $(e2.number) should be a hybrid.")
            # check that cn doesn't have any other parent than e1 and e2
            any(PhyloNetworks.getchild(e) ≡ cn && e !== e1 && e !==e2 for e in cn.edge) &&
                error("root has 2 hybrid edges, but their common child has an extra parent")
            PhyloNetworks.removeEdge!(cn,e1); PhyloNetworks.removeEdge!(cn,e2)
            PhyloNetworks.removeHybrid!(net,cn) # removes cn from net.hybrid, updates net.numHybrids
            cn.hybrid = false # !! allowrootbelow! not called: would require correct isChild1
            empty!(e1.node); empty!(e2.node)
            PhyloNetworks.deleteEdge!(net,e1,part=false); PhyloNetworks.deleteEdge!(net,e2,part=false)
            empty!(nodei.edge)
            PhyloNetworks.deleteNode!(net,nodei)
            net.root = findfirst(x -> x ≡ cn, net.node)
            Qdeleteleaf!(net, net.root, synth_e_dict; index=true, nofuse=nofuse, simplify=simplify,
                unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
            return nothing
        end
        # by now, nodei is the child of both e1 and e2
        p1 = PhyloNetworks.getparent(e1) # find both parents of hybrid leaf
        p2 = PhyloNetworks.getparent(e2)
        # remove node1 and both e1, e2
        sameparent = (p1≡p2) # 2-cycle
        PhyloNetworks.removeNode!(p1,e1);  PhyloNetworks.removeNode!(p2,e2) # perhaps useless
        PhyloNetworks.removeEdge!(p1,e1);  PhyloNetworks.removeEdge!(p2,e2)
        PhyloNetworks.deleteEdge!(net,e1,part=false); PhyloNetworks.deleteEdge!(net,e2,part=false)
        if net.root==i net.root=getIndex(p1,net); end # should never occur though.
        PhyloNetworks.deleteNode!(net,nodei)
        # recursive call on both p1 and p2.
        Qdeleteleaf!(net, p1.number, synth_e_dict; nofuse = nofuse, simplify=simplify, unroot=unroot,
                    multgammas=multgammas, keeporiginalroot=keeporiginalroot)
        # p2 may have already been deleted: e.g. if sameparent, or other scenarios
        if !sameparent
          p2idx = findfirst(n -> n.number == p2.number, net.node)
          isnothing(p2idx) ||
            Qdeleteleaf!(net, p2idx,synth_e_dict; index=true, nofuse=nofuse, simplify=simplify,
                        unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
        end
    elseif !nofuse
        e1 = Qfuseedgesat!(i,net, synth_e_dict,multgammas) # fused edge
        if simplify && e1.hybrid # check for 2-cycle at new hybrid edge
            cn = PhyloNetworks.getchild(e1)
            e2 = PhyloNetworks.getpartneredge(e1, cn) # companion hybrid edge
            pn  = PhyloNetworks.getparent(e1)
            if pn ≡ PhyloNetworks.getparent(e2)
                # e1 and e2 have same child and same parent. Remove e1.
                e2.hybrid = false # assumes bicombining at cn: no third hybrid parent
                e2.isMajor = true
                e2.gamma = PhyloNetworks.addBL(e1.gamma, e2.gamma)
                PhyloNetworks.removeEdge!(pn,e1); PhyloNetworks.removeEdge!(cn,e1)
                PhyloNetworks.deleteEdge!(net,e1,part=false)
                PhyloNetworks.removeHybrid!(net,cn) # removes cn from net.hybrid, updates net.numHybrids
                cn.hybrid = false # !! allowrootbelow! not called: would require correct isChild1
                # call recursion again because pn and/or cn might be of degree 2 (or even 1).
                Qdeleteleaf!(net, cn.number, synth_e_dict; nofuse = nofuse, simplify=simplify, unroot=unroot,
                            multgammas=multgammas, keeporiginalroot=keeporiginalroot)
                pnidx = findfirst(n -> n.number == pn.number, net.node)
                isnothing(pnidx) ||
                Qdeleteleaf!(net, pnidx, synth_e_dict; index=true, nofuse=nofuse, simplify=simplify,
                    unroot=unroot, multgammas=multgammas, keeporiginalroot=keeporiginalroot)
            end
        end
    end
    return nothing
end

function Qfuseedgesat!(i::Integer, net::HybridNetwork, synth_e_dict, multgammas=false::Bool)
    i <= length(net.node) ||
      error("node index $i too large: only $(length(net.node)) nodes in the network.")
    nodei = net.node[i]
    length(nodei.edge) == 2 ||
      error("can't fuse edges at node number $(nodei.number): connected to $(length(nodei.edge)) edges.")
    !(nodei.edge[1].hybrid && nodei.edge[2].hybrid) ||
      error("can't fuse edges at node number $(nodei.number): connected to exactly 2 hybrid edges")
    j = argmax([e.number for e in nodei.edge])
    pe = nodei.edge[j] # edge to remove: pe.number > ce.number
pen=pe.number
pep=PN.getparent(pe).number
pec=PN.getchild(pe).number
#pee=deepcopy(pe)
    ce = nodei.edge[j==1 ? 2 : 1]
cee=deepcopy(ce)
    if pe.hybrid       # unless it's a hybrid: should be --tree--> node i --hybrid-->
        (ce,pe) = (pe,ce) # keep the hybrid edge: keep its isMajor
    end
    isnodeiparent = (nodei ≡ getparent(ce))
    (!ce.hybrid || isnodeiparent) ||
      error("node $(nodei.number) has 1 tree edge ($(pe.number)) and 1 hybrid edge ($(ce.number)), but is child of the hybrid edge.")
    pn = PN.getOtherNode(pe,nodei)
    PN.removeEdge!(nodei,ce) # perhaps useless. in case gc() on ith node affects its edges.
    PN.removeNode!(nodei,ce)
    PN.removeEdge!(pn,pe)
    PN.removeNode!(pn,pe)    # perhaps useless. in case gc() on pe affects its nodes.
    PN.setEdge!(pn,ce)
    PN.setNode!(ce,pn)       # pn comes 2nd in ce now: ce.node is: [original, pn]
    ce.isChild1 = isnodeiparent # to retain same direction as before.
    ce.length = PN.addBL(ce.length, pe.length)
    if multgammas
        ce.gamma = PN.multiplygammas(ce.gamma, pe.gamma)
    end
    if net.root==i # isnodeiparent should be true, unless the root and ce's direction were not in sync
        newroot = pn
        if newroot.leaf && !ce.hybrid # then reverse ce's direction. pn.leaf and ce.hybrid should never both occur!
            newroot = ce.node[1] # getOtherNode(ce, pn)
            ce.isChild1 = false
        end
        net.root = findfirst(isequal(newroot), net.node)
    end
    PN.deleteNode!(net,nodei)
    PN.deleteEdge!(net,pe,part=false) # do not update partitions. irrelevant for networks of level>1.

if !(PN.getchild(ce).leaf) 
p=PN.getparent(ce).number
c=PN.getchild(ce).number
f1=parse(BigInt,synth_e_dict[(cee.number,PN.getparent(cee).number,PN.getchild(cee).number)])
f2=parse(BigInt,synth_e_dict[(pen,pep,pec)])
synth_e_dict[(ce.number,p,c)]=string(f1+f2)#binary(ce.number,synth_e_dict[pep,pec])
synth_e_dict[(ce.number,c,p)]=string(f1+f2)#binary(ce.number,synth_e_dict[pep,pec])
end    
    return ce
end

function Qdeletehybridedge!(
    net::HybridNetwork,
    edge::PhyloNetworks.Edge,
    synth_e_dict,
    nofuse::Bool=false,
    unroot::Bool=false,
    multgammas::Bool=false,
    simplify::Bool=true,
    keeporiginalroot::Bool=false
)
    edge.hybrid || error("edge $(edge.number) has to be hybrid for deletehybridedge!")
    n1 = getchild(edge)  # child of edge, to be deleted unless nofuse
    n1.hybrid || error("child node $(n1.number) of hybrid edge $(edge.number) should be a hybrid.")
    n1degree = length(n1.edge)
    n2 = getparent(edge)  # parent of edge, to be deleted too
    n2degree = length(n2.edge)
    # next: keep hybrid node n1 if it has 4+ edges or if keepNode.
    #       otherwise: detach n1, then delete recursively
    delete_n1_recursively = false
    if n1degree < 3
        error("node $(n1.number) has $(length(n1.edge)) edges instead of 3+");
    # alternatively: error if degree < 2 or leaf,
    #   warning if degree=2 and internal node, then
    #   delete_n1_recursively = true # n1 doesn't need to be detached first
    elseif n1degree == 3 && !nofuse # then fuse 2 of the edges and detach n1
        delete_n1_recursively = true
        pe = nothing # will be other parent (hybrid) edge of n1
        ce = nothing # will be child edge of n1, to be merged with pe
        for e in n1.edge
            if e.hybrid && e!==edge && n1===getchild(e) pe = e; end
            if !e.hybrid || n1===getparent(e)  ce = e; end # does *not* assume correct ischild1 for tree edges :)
        end
        pn = PhyloNetworks.getparent(pe); # parent node of n1, other than n2
        atRoot = (net.node[net.root] ≡ n1) # n1 should not be root, but if so, pn will be new root
        # if pe may contain the root, then allow the root on ce and below
        if pe.containRoot
            PhyloNetworks.allowrootbelow!(ce) # warning: assumes correct `ischild1` for ce and below
        end
        # next: replace ce by pe+ce, detach n1 from pe & ce, remove pe from network.
        ce.length =PhyloNetworks.addBL(ce.length, pe.length)
        if multgammas
            ce.gamma = PhyloNetworks.multiplygammas(ce.gamma, pe.gamma)
        end

ce1=deepcopy(ce)
pe1=deepcopy(pe)

        PhyloNetworks.removeNode!(n1,ce) # ce now has 1 single node cn
        PhyloNetworks.setNode!(ce,pn)    # ce now has 2 nodes in this order: cn, pn
        ce.isChild1 = true
        PhyloNetworks.setEdge!(pn,ce)

if !(PN.getchild(ce).leaf) 
f1=parse(BigInt,synth_e_dict[(pe1.number,PN.getparent(pe1).number, PN.getchild(pe1).number)])
f2=parse(BigInt,synth_e_dict[(ce1.number,PN.getparent(ce1).number, PN.getchild(ce1).number)])
synth_e_dict[(ce.number,PN.getparent(pe).number,PN.getchild(ce).number)]=string(f1+f2)
synth_e_dict[(ce.number,PN.getchild(ce).number,PN.getparent(pe).number)]=string(f1+f2)
end

        PhyloNetworks.removeEdge!(pn,pe)
        # if (pe.number<ce.number) ce.number = pe.number; end # bad to match edges between networks
        PhyloNetworks.removeEdge!(n1,pe); PhyloNetworks.removeEdge!(n1,ce) # now n1 attached to edge only
        PhyloNetworks.deleteEdge!(net,pe,part=false) # decreases net.numedges   by 1
        PhyloNetworks.removeHybrid!(net,n1) # decreases net.numhybrids by 1
        n1.hybrid = false
        edge.hybrid = false; edge.isMajor = true
        # n1 not leaf, and not added to net.leaf, net.numtaxa unchanged
        if atRoot
            i = findfirst(x -> x===pn, net.node)
            i !== nothing || error("node $(pn.number) not in net!")
            net.rooti = i
        end
        # below: we will need to delete n1 recursively (hence edge)
    else # n1 has 4+ edges (polytomy) or 3 edges but we want to keep it anyway:
        # keep n1 but detach it from 'edge', set its remaining parent to major tree edge
        pe = getpartneredge(edge, n1) # partner edge: keep it this time
        if !pe.isMajor pe.isMajor=true; end
        pe.hybrid = false
        # note: pe.gamma *not* set to 1.0 here
        PhyloNetworks.removeEdge!(n1,edge) # does not update n1.hybrid at this time
        PhyloNetworks.removeHybrid!(net,n1) # removes n1 from net.hybrid, updates net.numhybrids
        n1.hybrid = false
        if pe.containRoot
            PhyloNetworks.allowrootbelow!(pe) # warning: assumes correct `ischild1` for pe and below
        end
        # below: won't delete n1, delete edge instead
    end

    formernumhyb = net.numHybrids
    # next: delete n1 recursively, or delete edge and delete n2 recursively.
    # keep n2 if it has 4+ edges (or if nofuse). 1 edge should never occur.
    #       If root, would have no parent: treat network as unrooted and change the root.
    if delete_n1_recursively
        Qdeleteleaf!(net, n1.number,synth_e_dict; index=false, nofuse=nofuse,
                    simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
    # else: delete "edge" then n2 as appropriate
    elseif n2degree == 1
        error("node $(n2.number) (parent of hybrid edge $(edge.number) to be deleted) has 1 edge only!")
    else

        # fixit: if n2degree == 2 && n2 === net.node[net.rooti] and
        #        if we want to keep original root: then delete edge but keep n2
        # detach n2 from edge, remove hybrid 'edge' from network
        PhyloNetworks.removeEdge!(n2,edge)

if n2.hybrid
    p1=PN.getparent(n2.edge[1])
    c1=PN.getchild(n2.edge[1])
    p2=PN.getparent(n2.edge[2])
    c2=PN.getchild(n2.edge[2])
    if !(c1.leaf) && !(c2.leaf)
    f1=parse(BigInt,synth_e_dict[(n2.edge[1].number,p1.number,c1.number)])
    f2=parse(BigInt,synth_e_dict[(n2.edge[2].number,p2.number,c2.number)])
    synth_e_dict[(n2.edge[1].number,p1.number,p2.number)]=string(f1+f2)
    synth_e_dict[(n2.edge[1].number,p2.number,p1.number)]=string(f1+f2)
    end
else
    p1=PN.getparent(n2.edge[1])
    c1=PN.getchild(n2.edge[1])
    p2=PN.getparent(n2.edge[2])
    c2=PN.getchild(n2.edge[2])
    if !(c1.leaf) && !(c2.leaf)
    f1=parse(BigInt,synth_e_dict[(n2.edge[1].number,p1.number,c1.number)])
    f2=parse(BigInt,synth_e_dict[(n2.edge[2].number,p2.number,c2.number)])
        if (p1,c1)==(p2,c2)
            synth_e_dict[(n2.edge[1].number,p1.number,c2.number)]=string(f1)
            synth_e_dict[(n2.edge[1].number,c2.number,p1.number)]=string(f1)
        else
            synth_e_dict[(n2.edge[1].number,p1.number,c2.number)]=string(f1+f2)
            synth_e_dict[(n2.edge[1].number,c2.number,p1.number)]=string(f1+f2)
        end
    end
end
        PhyloNetworks.deleteEdge!(net,edge,part=false)
        # remove n2 as appropriate later (recursively)   
        Qdeleteleaf!(net, n2.number,synth_e_dict; index=false, nofuse=nofuse,
                    simplify=simplify, unroot=unroot, multgammas=multgammas,
                    keeporiginalroot=keeporiginalroot)
    end
    if net.numHybrids != formernumhyb # deleteleaf! does not update containroot
        PhyloNetworks.allowrootbelow!(net)
    end
    return net
end

#=
    """
        binary(n::Int, binstr::String) -> String

    Set the nth bit of a binary string to 1, returning the updated binary string.
    If the input string is empty, returns a string of length `n` with the first bit set to 1 and the rest zeros.

    Arguments:
    - n::Int : The bit position to set (1-based, from the rightmost bit).
    - binstr::String : A binary string. Can be empty.

    Returns:
    - String : The binary string with the nth bit set.

    Examples:
    - binary(3, "10") # returns "110"
    - binary(5, "") # returns "10000"
    """
    function bina11ry(n::Int, binstr::String)
        if isempty(binstr)
            # If empty, create 1 followed by (n-1) zeros
            return "1" * repeat("0", n-1)
        end    
        # Convert binary string to integer
        original = parse(Int, binstr)#; base=2)    
        # Set the nth bit to 1
        new_value = original | (1 << (n-1))    
        # Determine the minimum length: keep at least n bits
        min_len = max(length(binstr), n)  
        # Convert to binary string
        new_binstr = string(new_value)#, base=2)
        # Pad with leading zeros if necessary
        if length(new_binstr) < min_len
            new_binstr = lpad(new_binstr, min_len, '0')
        end
        return new_binstr
    end
=#

#=
    """
        parameterDictionary1(net, inheritancecorrelation; gammaSymbol::String="r")

    Creates a dictionary mapping inheritance probabilities (γ) and inheritance correlation (ρ) 
    in a phylogenetic network to symbolic labels.

    # Description
    This function generates a dictionary of symbolic parameter labels for use in 
    algebraic or likelihood-based calculations.

    - **Inheritance probabilities (γ):**  
    For each hybrid node, the two incoming edges are assigned labels:
    - `"r_{j}"` for the first edge  
    - `"(1-r_{j})"` for the second edge  
    where `j` is the index of the hybrid node (1-based).  
    Each `γ` value is rounded to `dpoints` decimal places before being used as a dictionary key.

    - **Inheritance correlation (ρ):**  
    The provided `inheritancecorrelation` value is assigned the label `"&rho"`.  
    Its complement, `1 - ρ`, is rounded and assigned the label `"1-&rho"`.

    This function assumes each hybrid node has exactly two incoming edges, and will 
    throw an error otherwise.

    # Arguments
    - `net::PhyloNetworks.HybridNetwork`: The input network object.
    - `inheritancecorrelation::Float64`: The inheritance correlation value.
    - `gammaSymbol::String="r"`: Optional prefix for inheritance probability labels.

    # Returns
    - `Dict{Any, String}`: A dictionary mapping numerical parameter values to symbolic labels.
    """
    function paramesssterDictionary1(net, inheritancecorrelation; gammaSymbol="r_"::String)
        dict = Dict()

        # Dictionary for inheritance probabilities (γ)
        hybridNodeNumbers = [n.number for n in net.node if n.hybrid]
        for (j, hybNode) in enumerate(hybridNodeNumbers)
            incoming = [e for e in net.edge if PhyloNetworks.getchild(e).number == hybNode]
            if length(incoming) != 2
                error("Hybrid node $hybNode has $(length(incoming)) incoming edges (expected 2).")
            end
            for (k, e) in enumerate(incoming)
                e.gamma = round(e.gamma, digits=dpoints)
                dict[e.gamma] = k == 1 ? "$gammaSymbol{$j}" : "(1-$gammaSymbol{$j})"
            end
        end

        # Inheritance correlation (ρ)
        dict[inheritancecorrelation] = "rho"
        dict[round(1 - inheritancecorrelation, digits=dpoints)] = "1-rho"
        
        return dict
    end
=#