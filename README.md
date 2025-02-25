# SymbolicQuartetNetworkCoal.jl

Aloha!

`SymbolicQuartetNetworkCoal.jl` is a Julia package that provides useful functions to study identifiability of phylogenetic trees and networks using quartet concordance factors with techniques from algebraic statistics. This package is constructed by extending another Julia package `QuartetNetworkGoodnessFit`, aka "Quarnet GoF" or simply "QGoF" (https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl).

## Installation
TBA

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
```@julia
julia> julia> network_expectedCF(ik1)
(PhyloNetworks.QuartetT{StaticArraysCore.MVector{3, Float64}}[4-taxon set number 1; taxon numbers: 1,2,3,4
data: [9.618622528537664e-35, 1.0, 9.618622528537664e-35], 4-taxon set number 2; taxon numbers: 1,2,3,5
data: [9.321874190057184e-13, 0.9999999999981356, 9.321874190057184e-13], 4-taxon set number 3; taxon numbers: 1,2,4,5
data: [3.0628420794597945e-8, 3.0628420794597945e-8, 0.9999999387431584], 4-taxon set number 4; taxon numbers: 1,3,4,5
data: [0.9999999999981356, 9.321874190057184e-13, 9.321874190057184e-13], 4-taxon set number 5; taxon numbers: 2,3,4,5
data: [3.0628420794597945e-8, 0.9999999387431584, 3.0628420794597945e-8], 4-taxon set number 6; taxon numbers: 1,2,3,6
data: [9.321874190057184e-13, 0.9999999999981356, 9.321874190057184e-13], 4-taxon set number 7; taxon numbers: 1,2,4,6
data: [3.0628420794597945e-8, 3.0628420794597945e-8, 0.9999999387431584], 4-taxon set number 8; taxon numbers: 1,3,4,6
data: [0.9999999999981356, 9.321874190057184e-13, 9.321874190057184e-13], 4-taxon set number 9; taxon numbers: 2,3,4,6
data: [3.0628420794597945e-8, 0.9999999387431584, 3.0628420794597945e-8], 4-taxon set number 10; taxon numbers: 1,2,5,6
data: [0.9999782998462678, 1.0800076806077543e-5, 1.0800076806077543e-5]  …  4-taxon set number 61; taxon numbers: 3,4,7,8
data: [0.9999999970451936, 1.4774031943872099e-9, 1.4774031943872099e-9], 4-taxon set number 62; taxon numbers: 1,5,7,8
data: [0.9999499336595082, 2.503317024587204e-5, 2.503317024587204e-5], 4-taxon set number 63; taxon numbers: 2,5,7,8
data: [0.9999499363732899, 2.5031813355072422e-5, 2.5031813355072422e-5], 4-taxon set number 64; taxon numbers: 3,5,7,8
data: [0.9999499336595082, 2.503317024587204e-5, 2.503317024587204e-5], 4-taxon set number 65; taxon numbers: 4,5,7,8
data: [0.9999499363732899, 2.5031813355072422e-5, 2.5031813355072422e-5], 4-taxon set number 66; taxon numbers: 1,6,7,8
data: [0.9999499336595082, 2.503317024587204e-5, 2.503317024587204e-5], 4-taxon set number 67; taxon numbers: 2,6,7,8
data: [0.9999499363732899, 2.5031813355072422e-5, 2.5031813355072422e-5], 4-taxon set number 68; taxon numbers: 3,6,7,8
data: [0.9999499336595082, 2.503317024587204e-5, 2.503317024587204e-5], 4-taxon set number 69; taxon numbers: 4,6,7,8
data: [0.9999499363732899, 2.5031813355072422e-5, 2.5031813355072422e-5], 4-taxon set number 70; taxon numbers: 5,6,7,8
data: [0.9999998894163947, 5.291742700465496e-9, 5.291742700465496e-9]], ["A", "B", "C", "D", "E", "F", "G", "H"])
```

```
SymbolicQuartetNetworkCoal.jl log
Timestamp: 2025-02-24T18:45:18.762
------------------------
General setting: 
Symbolic option: off
Store output in .csv file: off
Write Matlab file: off
Write Macaulay2 file: off
------------------------
Topology:
((C:1.85113,A:2.4774289):3.5732234,(((G:4.2852874,H:5.1839555):6.9904621,(((E:7.0906201,F:8.5153955):9.0431873)#H2:10.4526595::0.1021757)#H1:11.0652352::0.7983348):12.2439055,((#H2:13.0069595::0.8978243,(B:14.8312632,D:15.7456604):16.0949443):17.4547141,#H1:18.3397273::0.2016652):19.5821649):20.5231182);
((C:t_{1},A:t_{2}):t_{3},(((G:t_{4},H:t_{5}):t_{6},(((E:t_{7},F:t_{8}):t_{9})#H2:t_{10}::r_{1})#H1:11-&rho652352::r_{2}):t_{12},((#H2:t_{13}::(1-r_{1}),(B:t_{14},D:t_{15}):t_{16}):t_{17},#H1:t_{18}::(1-r_{2})):t_{19}):t_{20});
------------------------
Parameters:
Parameter		Value
t_{1}		1.85113
t_{2}		2.4774289
t_{3}		3.5732234
t_{4}		4.2852874
t_{5}		5.1839555
t_{6}		6.9904621
t_{7}		7.0906201
t_{8}		8.5153955
t_{9}		9.0431873
t_{10}		10.4526595
r_{1}		0.1021757
t_{11}		11.0652352
r_{2}		0.7983348
t_{12}		12.2439055
t_{13}		13.0069595
(1-r_{1})		0.8978243
t_{14}		14.8312632
t_{15}		15.7456604
t_{16}		16.0949443
t_{17}		17.4547141
t_{18}		18.3397273
(1-r_{2})		0.2016652
t_{19}		19.5821649
t_{20}		20.5231182
&rho		0
------------------------
Concordance factor:
Quartet		Formula
AB|CD		(exp(-77.2281649)/3)
AC|BD		(1-2*exp(-77.2281649)/3)
AD|BC		(exp(-77.2281649)/3)
.
.
.
```

```@julia
julia> network_expectedCF(ik1,savecsv=true)
```

#### Symbolic formulas
```@julia
julia> julia> network_expectedCF(ik1,savecsv=true,symbolic=true)
```

```
SymbolicQuartetNetworkCoal.jl log
Timestamp: 2025-02-24T18:46:51.858
------------------------
General setting: 
Symbolic option: on
Store output in .csv file: on
Write Matlab file: off
Write Macaulay2 file: off
------------------------
Topology:
((C:1.7566785,A:2.9934842):3.1356123,(((G:4.551512,H:5.1216977):6.8029044,(((E:7.3012429,F:8.3767442):9.6975132)#H2:10.7946048::0.2514125)#H1:11.5440906::0.7356495):12.5929871,((#H2:13.9026236::0.7485875,(B:14.9021465,D:15.328871):16.1521944):17.0057032,#H1:18.8169413::0.2643505):19.5466342):20.6892273);
((C:t_{1},A:t_{2}):t_{3},(((G:t_{4},H:t_{5}):t_{6},(((E:t_{7},F:t_{8}):t_{9})#H2:t_{10}::r_{1})#H1:t_{11}::r_{2}):t_{12},((#H2:t_{13}::(1-r_{1}),(B:t_{14},D:t_{15}):t_{16}):t_{17},#H1:t_{18}::(1-r_{2})):t_{19}):t_{20});
------------------------
Parameters:
Parameter		Value
t_{1}		1.7566785
t_{2}		2.9934842
t_{3}		3.1356123
t_{4}		4.551512
t_{5}		5.1216977
t_{6}		6.8029044
t_{7}		7.3012429
t_{8}		8.3767442
t_{9}		9.6975132
t_{10}		10.7946048
r_{1}		0.2514125
t_{11}		11.5440906
r_{2}		0.7356495
t_{12}		12.5929871
t_{13}		13.9026236
(1-r_{1})		0.7485875
t_{14}		14.9021465
t_{15}		15.328871
t_{16}		16.1521944
t_{17}		17.0057032
t_{18}		18.8169413
(1-r_{2})		0.2643505
t_{19}		19.5466342
t_{20}		20.6892273
&rho		0
------------------------
Concordance factor:
Quartet		Formula
AB|CD		(exp(-t_{3}-t_{16}-t_{17}-t_{19}-t_{20})/3)
AC|BD		(1-2*exp(-t_{3}-t_{16}-t_{17}-t_{19}-t_{20})/3)
AD|BC		(exp(-t_{3}-t_{16}-t_{17}-t_{19}-t_{20})/3)
AB|CE		(((exp(-t_{3}-t_{20})/3)*r_{2}+(exp(-t_{3}-t_{19}-t_{20})/3)*(1
.
.
.
```
### Creating Macaulay2 and Matlab input file
```@julia
julia> network_expectedCF(ik1,savecsv=true,symbolic=true,macaulay=true,matlab=true)
```

```@julia
julia> network_expectedCF(ik1,savecsv=true,symbolic=false,macaulay=true,matlab=true)
ERROR: symbolic must be set to true.
Stacktrace:
 [1] error(s::String)
   @ Base ./error.jl:33
```

### Visualizing network with parameter names
```@julia
julia> el=makeEdgeLabel(ik1)
20×2 DataFrame
 Row │ number   label
     │ Integer  String
─────┼─────────────────
   1 │       1  t_{1}
   2 │       2  t_{2}
   3 │       3  t_{3}
   4 │       4  t_{4}
   5 │       5  t_{5}
   6 │       6  t_{6}
   7 │       7  t_{7}
   8 │       8  t_{8}
   9 │       9  t_{9}
  10 │      10  t_{10}
  11 │      11  t_{11}
  12 │      12  t_{12}
  13 │      13  t_{13}
  14 │      14  t_{14}
  15 │      15  t_{15}
  16 │      16  t_{16}
  17 │      17  t_{17}
  18 │      18  t_{18}
  19 │      19  t_{19}
  20 │      20  t_{20}
```

```@julia  
julia> using PhyloPlots
julia> plot(ik1,edgelabel=elabels)
```
![Alt text](example/edgelabeled-ik1.png)


## Citation