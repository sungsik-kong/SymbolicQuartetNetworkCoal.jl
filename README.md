# SymbolicQuartetNetworkCoal.jl

`SymbolicQuartetNetworkCoal.jl` is a Julia package that provides useful functions to study identifiability of phylogenetic trees and networks using quartet concordance factors with techniques from algebraic statistics. This package is constructed by extending another Julia package `QuartetNetworkGoodnessFit`, aka "Quarnet GoF" or simply "QGoF" (https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl).

## Main functions
### Parametrized topology with random values
```@julia
julia> ik1=SymbolicQuartetNetworkCoal.readTopologyrand("((C,A),(((G,H),(((E,F))#H2)#H1),((#H2,(B,D)),#H1)));")
PhyloNetworks.HybridNetwork, Rooted Network
20 edges
19 nodes: 8 tips, 2 hybrid nodes, 9 internal tree nodes.
tip labels: C, A, G, H, ...
((C:1.851,A:2.477):3.573,(((G:4.285,H:5.184):6.99,(((E:7.091,F:8.515):9.043)#H2:10.453::0.102)#H1:11.065::0.798):12.244,((#H2:13.007::0.898,(B:14.831,D:15.746):16.095):17.455,#H1:18.34::0.202):19.582):20.523);
```
### Obtain formulas for CF computation for each quartet
#### Numerical formulas

#### Symbolic formulas
```@julia
julia> using Pkg
julia> Pkg.add("PhyNEST")
```
### Creating Macaulay2 and Matlab input file
```@julia
julia> using Pkg
julia> Pkg.add("PhyNEST")
```
### Visualizing network with parameter names

## Citation