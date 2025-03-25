    #test written by Sungsik Kong 2025

    @testset begin
        threshold=0.05 #we want the absolute difference between the true and computed values to be <threshold
        ih=0#inheritancecorrelation-works when 0/1 but not works for other values, like 0.9, why?
        nreps=1#number of retplicates for m2,matlab
        count=0

        #true vs numeric
        #filename="ik3.txt"
        filename="topologies_n5_l2.txt"
#        filename="sim_nets.txt" #2000 simulated network from Morley et al. 

        
        open(filename, "r") do file
            for i in eachline(file)

                df_true=DataFrame(Split=String[], CF_true=Float64[])
                df_numeric=DataFrame(Split=String[], CF_numeric=String[], CF_val_from_numeric=Float64[])
                df_symbolic=DataFrame(Split=String[], CF_symbolic=String[], CF_sym_substituted=String[], CF_val_from_symbolic=Float64[])
                
                count+=1
                println("Working on case $count...")
                println("Currenet tree:\n$i")
                #ik1=SymbolicQuartetNetworkCoal.readTopologyrand("((C,A),(((G,H),(((E,F))#H2)#H1),((#H2,(B,D)),#H1)));")
                ik1=SymbolicQuartetNetworkCoal.readTopologyrand(i)
                #Q0,T0=SymbolicQuartetNetworkCoal.network_expectedCF_formulas(ik1,filename="temp")

                
                #numeric res stored in df_numeric
                Q0,T0=SymbolicQuartetNetworkCoal.network_expectedCF_formulas(ik1,symbolic=false,savecsv=true,filename="temp-num",inheritancecorrelation=ih)
                df_num=CSV.read("temp-num.csv",DataFrame,header=false) #make sure CSV file is reliable
                for nthrow in 1:nrow(df_num) 
                    push!(df_numeric,(df_num[nthrow,1],df_num[nthrow,2],eval(Meta.parse(df_num[nthrow,2])))) 
                end #compute numeric formula and record the value

                #true res stored in df_true
                for nthq in 1:binomial(length(T0),4)
                    push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[2]])|$(T0[Q0[nthq].taxonnumber[3]])$(T0[Q0[nthq].taxonnumber[4]])",Q0[nthq].data[1]))
                    push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[3]])|$(T0[Q0[nthq].taxonnumber[2]])$(T0[Q0[nthq].taxonnumber[4]])",Q0[nthq].data[2]))
                    push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[4]])|$(T0[Q0[nthq].taxonnumber[2]])$(T0[Q0[nthq].taxonnumber[3]])",Q0[nthq].data[3]))
                end


                #symbolic res stored in df_symbolic
                Q2,T2=SymbolicQuartetNetworkCoal.network_expectedCF_formulas(ik1,symbolic=true,savecsv=true,filename="temp-sym",inheritancecorrelation=ih)
                df_sym=CSV.read("temp-sym.csv",DataFrame,header=false)
                dict=SymbolicQuartetNetworkCoal.parameterDictionary(ik1,ih)
                revdict=Dict(value => key for (key, value) in dict)
                
                
                for nthrow in 1:nrow(df_sym)
                    df_replace="$(df_sym[nthrow,2])"
                    for edgenum in 1:length(ik1.edge) df_replace=(replace(df_replace,"t_{$edgenum}"=>revdict["t_{$edgenum}"])) end #replace taus to edge lengths
                    for retnum in 1:length(ik1.hybrid) df_replace=(replace(df_replace,"r_{$retnum}"=>revdict["r_{$retnum}"])) end #replace gammas to inheritance probabilities
                    df_replace=(replace(df_replace,"rho"=>ih)) #replace rho by ih
                    push!(df_symbolic,(df_sym[nthrow,1],df_sym[nthrow,2],df_replace,(eval(Meta.parse(df_replace))))) #compute new numeric formula and record the value
                end

                #merge all results based on the :Split entry
                df_test1=outerjoin(df_true,df_numeric,df_symbolic,on=:Split)
                CSV.write("temp-$(filename)_$count.csv",df_test1,header=false) 
                #evaluate values
                for nthrow in 1:nrow(df_test1)
                    @test abs(df_test1[nthrow,2]-df_test1[nthrow,4]) < threshold || error("True and numerical values do not match") #checking if the true value and numerical computation matches
                    @test abs(df_test1[nthrow,2]-df_test1[nthrow,7]) < threshold || error("True and symbolic values do not match") #checking if the true value and symbolic computation matches
                    @test abs(df_test1[nthrow,4]-df_test1[nthrow,7]) < threshold || error("Numerical and symbolic values do not match") #checking if the numerical computation and symbolic computation matches
                end        

                #generate nreps macaluay and matlab files one a same topoplogy but randomized parameters
                for x in 1:nreps
                    #ik1=readTopologyrand("((C,A),(((G,H),(((E,F))#H2)#H1),((#H2,(B,D)),#H1)));")
                    ik1=SymbolicQuartetNetworkCoal.readTopologyrand(i)
                    QM,TM=network_expectedCF_formulas(ik1,symbolic=true,savecsv=true,filename="temp-macaulay-matlab-$x",inheritancecorrelation=ih,macaulay=true,matlab=true)
                end
                for y in 1:(nreps-1)
                    for z in 2:nreps
                        @test filecmp("temp-macaulay-matlab-$y.m2.txt", "temp-macaulay-matlab-$z.m2.txt")==true || error("Macaulay2 files do not match") 
                        @test filecmp("temp-macaulay-matlab-$y.matlab.txt", "temp-macaulay-matlab-$z.matlab.txt")==true || error("Matlab files do not match") 
                    end
                end
                foreach(rm, filter(startswith("temp"), readdir()))#remove all temp files
                #foreach(rm, filter(startswith("$filename"), readdir()))#remove all temp files

            end
        end
    end




@testset "makeEdgeLabel Tests" begin
    # Test network: ((A:1.0,B:2.0):3.0,C:4.0);
    net = readTopology("((A:1.0,B:2.0):3.0,C:4.0);")
    
    # Test 1: Default behavior (exclude terminal edges)
    df_default = makeEdgeLabel(net)
    @test nrow(df_default) == 1  # Only the internal edge should be included
    @test df_default.number[1] == 3  # Assuming edge 3 is the internal one (may vary)
    @test df_default.label[1] == "t_{3}"
    
    # Test 2: Include terminal edges
    df_all = makeEdgeLabel(net, showTerminalEdgeLabels=true)
    @test nrow(df_all) == 4  # All edges (A, B, and internal)
    @test Set(df_all.number) == Set([1, 2, 3, 4])  # Exact numbers may vary
    @test Set(df_all.label) == Set(["t_{1}", "t_{2}", "t_{3}", "t_{4}"])
    
    # Test 3: Empty network case
    net_empty = HybridNetwork()  # Empty network
    df_empty = makeEdgeLabel(net_empty)
    @test nrow(df_empty) == 0
    
    # Test 4: All leaves excluded
    net_single = readTopology("(A:1.0);")
    df_single = makeEdgeLabel(net_single)
    @test nrow(df_single) == 0  # No non-terminal edges
end