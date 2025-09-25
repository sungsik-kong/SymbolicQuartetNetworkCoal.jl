#test written by Sungsik Kong 2025
@testset begin
    threshold=0.001 #we want the absolute difference between the true and computed values to be <threshold
    ih=0.1 #inheritancecorrelation
    eLab="t_"
    gammaSymbol="r_"
    dpoints=3
    filename=["topologies_n5_l1.txt","topologies_n5_l2.txt"]
    count1=0

    for afile in filename
        count1+=1
        open(afile, "r") do file
            count=0
            for i in eachline(file)
                df_true=DataFrame(Split=String[], CF_true=Float64[])
                df_numeric=DataFrame(Split=String[], CF_numeric=String[], CF_val_from_numeric=Float64[])
                df_symbolic=DataFrame(Split=String[], CF_symbolic=String[], CF_sym_substituted=String[], CF_val_from_symbolic=Float64[])
                
                ###ioprint###
                    count+=1
                    println("Case [$count1, $count]: $i")
                    testnet=SymbolicQuartetNetworkCoal.readTopologyrand(i)
                ###ioprint###
        
                #compute true CFs and numeric formulas
                Q0,T0,df_num=SymbolicQuartetNetworkCoal.network_expectedCF_formulas(testnet,symbolic=false,inheritancecorrelation=ih)

                #true CFs stored in df_true
                for nthq in 1:binomial(length(T0),4)
                    push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[2]])|$(T0[Q0[nthq].taxonnumber[3]])$(T0[Q0[nthq].taxonnumber[4]])",Q0[nthq].data[1]))
                    push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[3]])|$(T0[Q0[nthq].taxonnumber[2]])$(T0[Q0[nthq].taxonnumber[4]])",Q0[nthq].data[2]))
                    push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[4]])|$(T0[Q0[nthq].taxonnumber[2]])$(T0[Q0[nthq].taxonnumber[3]])",Q0[nthq].data[3]))
                end
                
                #numeric formulas stored in df_numeric
                for nthrow in 1:nrow(df_num) 
                    push!(df_numeric,(df_num[nthrow,1],df_num[nthrow,2],eval(Meta.parse(df_num[nthrow,2])))) 
                    #compute numeric formula and record the value
                end 

                #symbolic formulas stored in df_numeric
                Q2,T2,df_sym=SymbolicQuartetNetworkCoal.network_expectedCF_formulas(testnet,symbolic=true,inheritancecorrelation=ih)
                revdict=Dict()


 
    # Dictionary for inheritance probabilities (γ)
    hybridNodeNumbers = [n.number for n in testnet.node if n.hybrid]
    for j in 1:testnet.numHybrids
        hybNode = hybridNodeNumbers[j]
        visitCount = 1
        for e in testnet.edge
            if PhyloNetworks.getchild(e).number == hybNode
                e.gamma = round(e.gamma, digits=dpoints)
                if visitCount == 1
                    revdict["$gammaSymbol{$j}"] = e.gamma
                    visitCount += 1
                elseif visitCount == 2
                    revdict["(1-$gammaSymbol{$j})"] = e.gamma
                else
                    error("Hybrid node $(hybNode) has more than 2 incoming edges.")
                end
            end
        end
    end
    # Inheritance correlation (ρ)
    #inheritancecorrelation = round(inheritancecorrelation, digits=dpoints)
    revdict["&rho"] = ih
    revdict["1-&rho"] = round(1 - ih, digits=dpoints)



                    for e in testnet.edge 
                            revdict["$eLab{$(e.number)}"] = e.length
                    end

                for nthrow in 1:nrow(df_sym)
                    df_replace="$(df_sym[nthrow,2])"
                    for edgenum in 1:length(testnet.edge) df_replace=(replace(df_replace,"t_{$edgenum}"=>revdict["t_{$edgenum}"])) end #replace taus to edge lengths
                    for retnum in 1:length(testnet.hybrid) df_replace=(replace(df_replace,"r_{$retnum}"=>revdict["r_{$retnum}"])) end #replace gammas to inheritance probabilities
                    df_replace=(replace(df_replace,"rho"=>ih)) #replace rho by ih
                    push!(df_symbolic,(df_sym[nthrow,1],df_sym[nthrow,2],df_replace,(eval(Meta.parse(df_replace))))) #compute new numeric formula and record the value
                end

                #merge all results based on the :Split entry
                df_test1=outerjoin(df_true,df_numeric,df_symbolic,on=:Split)
                CSV.write("temp-$(filename)_$count.csv",df_test1,header=false) 
                havingerror=false
                #evaluate values
                for nthrow in 1:nrow(df_test1)
                    a=abs(df_test1[nthrow,2]-df_test1[nthrow,4]) 
                    b=abs(df_test1[nthrow,2]-df_test1[nthrow,7])
                    c=abs(df_test1[nthrow,4]-df_test1[nthrow,7])
                    @test a < threshold || error("True and numerical values do not match: $a") #checking if the true value and numerical computation matches
                    @test b < threshold || error("True and symbolic values do not match: $b") #checking if the true value and symbolic computation matches
                    @test c < threshold || error("Numerical and symbolic values do not match: $c") #checking if the numerical computation and symbolic computation matches
                    if !(a < threshold && b <threshold && c<threshold)
                        havingerror=true
                    end
                end    
                if !(havingerror) 
                    rm("temp-$(filename)_$count.csv")
                end
            end
        end
    end
end