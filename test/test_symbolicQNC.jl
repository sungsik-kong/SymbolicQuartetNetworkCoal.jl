#test written by Sungsik Kong 2025
const eLab="t_"
const rLab="r_"
const dpoints=3

@testset begin
    threshold=0.0000000001 #we want the absolute difference between the true and computed values to be <threshold
    ih=0.1 #inheritancecorrelation
    filename=["topologies_n5_l1.txt","topologies_n5_l2.txt"]
    count1=0

    function parameterDictionary1REV(testnet, ih)
        revdict = Dict()
        
        for e in testnet.edge 
            revdict["$eLab{$(e.number)}"] = e.length
        end
        # Dictionary for inheritance probabilities (γ)
        hybridNodeNumbers = [n.number for n in testnet.node if n.hybrid]
        for j in 1:testnet.numHybrids
            hybNode = hybridNodeNumbers[j]
            visitCount = 1
            for e in testnet.edge
                if PhyloNetworks.getchild(e).number == hybNode
                    e.gamma = round(e.gamma, digits=dpoints)
                    if visitCount == 1
                        revdict["$rLab{$j}"] = e.gamma
                        visitCount += 1
                    elseif visitCount == 2
                        revdict["(1-$rLab{$j})"] = e.gamma
                    else
                        error("Hybrid node $(hybNode) has more than 2 incoming edges.")
                    end
                end
            end
        end

        # Inheritance correlation (ρ)
        revdict["&rho"] = ih
        revdict["1-&rho"] = round(1 - ih, digits=dpoints)
        
        return revdict
    end

    for afile in filename
        count1+=1
        open(afile, "r") do file
            count=0
            for i in eachline(file)
                ###ioprint###
                count+=1
                println("[Level-$count1 network, $count]: $i")
                ###ioprint###
                df_true=DataFrame(Split=String[], CF_true=Float64[])
                df_numeric=DataFrame(Split=String[], CF_numeric=String[], CF_val_from_numeric=Float64[])
                df_symbolic=DataFrame(Split=String[], CF_symbolic=String[], CF_sym_substituted=String[], CF_val_from_symbolic=Float64[])
                
                testnet=SymbolicQuartetNetworkCoal.readTopologyrand(i)
                revdict=parameterDictionary1REV(testnet,ih)                    

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
                    #compute numeric formula and record the value
                    push!(df_numeric,(df_num[nthrow,1],df_num[nthrow,2],eval(Meta.parse(df_num[nthrow,2])))) 
                end 

                #symbolic formulas stored in df_numeric
                Q2,T2,df_sym=SymbolicQuartetNetworkCoal.network_expectedCF_formulas(testnet,symbolic=true,inheritancecorrelation=ih)
                for nthrow in 1:nrow(df_sym)
                    df_replace="$(df_sym[nthrow,2])"
                    for edgenum in 1:length(testnet.edge) df_replace=(replace(df_replace,"$eLab{$edgenum}"=>revdict["$eLab{$edgenum}"])) end #replace taus to edge lengths
                    for retnum in 1:length(testnet.hybrid) df_replace=(replace(df_replace,"$rLab{$retnum}"=>revdict["$rLab{$retnum}"])) end #replace gammas to inheritance probabilities
                    df_replace=(replace(df_replace,"rho"=>ih)) #replace rho by ih
                    push!(df_symbolic,(df_sym[nthrow,1],df_sym[nthrow,2],df_replace,(eval(Meta.parse(df_replace))))) #compute new numeric formula and record the value
                end

                #merge all results based on the :Split entry and evaluate
                havingerror=false
                df_test1=outerjoin(df_true,df_numeric,df_symbolic,on=:Split)
                df_test1[!,:truexnumeric] .= 1.0
                df_test1[!,:truexsymbolic] .= 1.0
                df_test1[!,:numericxsymbolic] .= 1.0
                for nthrow in 1:nrow(df_test1)
                    df_test1[nthrow,8]=abs(df_test1[nthrow,2]-df_test1[nthrow,4]) 
                    df_test1[nthrow,9]=abs(df_test1[nthrow,2]-df_test1[nthrow,7])
                    df_test1[nthrow,10]=abs(df_test1[nthrow,4]-df_test1[nthrow,7])
                    @test df_test1[nthrow,8] < threshold || error("True and numerical values do not match: $(df_test1[nthrow,8])") #checking if the true value and numerical computation matches
                    @test df_test1[nthrow,9] < threshold || error("True and symbolic values do not match: $(df_test1[nthrow,9])") #checking if the true value and symbolic computation matches
                    @test df_test1[nthrow,10] < threshold || error("Numerical and symbolic values do not match: $(df_test1[nthrow,10])") #checking if the numerical computation and symbolic computation matches
                    if !(df_test1[nthrow,8] < threshold && df_test1[nthrow,9] <threshold && df_test1[nthrow,10]<threshold)
                        havingerror=true
                    end
                end    

                if (havingerror) 
                    select!(df_test1, [:truexnumeric,:truexsymbolic,:numericxsymbolic], Not([:truexnumeric,:truexsymbolic,:numericxsymbolic]))
                    CSV.write("temp-$(afile)_$count.csv",df_test1,header=false) 
                else
                    isfile("temp-$(afile)_$count.csv") ? rm("temp-$(afile)_$count.csv") : continue
                end
            end
        end
    end
end