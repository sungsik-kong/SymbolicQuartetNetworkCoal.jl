    #test written by Sungsik Kong 2025

    @testset begin
        threshold=0.01 #we want the absolute difference between the true and computed values to be <threshold
        ih=0#inheritancecorrelation-works when 0/1 but not works for other values, like 0.9, why?
        nreps=30#number of retplicates for m2,matlab

        df_true=DataFrame(Split=String[], CF_true=Float64[])
        df_numeric=DataFrame(Split=String[], CF_numeric=String[], CF_val_from_numeric=Float64[])
        df_symbolic=DataFrame(Split=String[], CF_symbolic=String[], CF_sym_substituted=String[], CF_val_from_symbolic=Float64[])
        
        #true vs numeric
        ik1=SymbolicQuartetNetworkCoal.readTopologyrand("((C,A),(((G,H),(((E,F))#H2)#H1),((#H2,(B,D)),#H1)));")
        Q0,T0=SymbolicQuartetNetworkCoal.network_expectedCF(ik1,filename="temp")

        #true res stored in df_true
        for nthq in 1:binomial(length(T0),4)
            push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[2]])|$(T0[Q0[nthq].taxonnumber[3]])$(T0[Q0[nthq].taxonnumber[4]])",Q0[nthq].data[1]))
            push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[3]])|$(T0[Q0[nthq].taxonnumber[2]])$(T0[Q0[nthq].taxonnumber[4]])",Q0[nthq].data[2]))
            push!(df_true,("$(T0[Q0[nthq].taxonnumber[1]])$(T0[Q0[nthq].taxonnumber[4]])|$(T0[Q0[nthq].taxonnumber[2]])$(T0[Q0[nthq].taxonnumber[3]])",Q0[nthq].data[3]))
        end

        #numeric res stored in df_numeric
        Q1,T1=SymbolicQuartetNetworkCoal.network_expectedCF(ik1,symbolic=false,savecsv=true,filename="temp-num",inheritancecorrelation=ih)
        df_num=CSV.read("temp-num.csv",DataFrame,header=false) #make sure CSV file is reliable
        for nthrow in 1:nrow(df_num) push!(df_numeric,(df_num[nthrow,1],df_num[nthrow,2],eval(Meta.parse(df_num[nthrow,2])))) end #compute numeric formula and record the value

        #symbolic res stored in df_symbolic
        Q2,T2=SymbolicQuartetNetworkCoal.network_expectedCF(ik1,symbolic=true,savecsv=true,filename="temp-sym",inheritancecorrelation=ih)
        df_sym=CSV.read("temp-sym.csv",DataFrame,header=false)
        dict=SymbolicQuartetNetworkCoal.dictionary(ik1,ih)
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
        
        #evaluate values
        for nthrow in 1:nrow(df_test1)
            @test abs(df_test1[nthrow,2]-df_test1[nthrow,4]) < threshold || error("True and numerical values do not match") #checking if the true value and numerical computation matches
            @test abs(df_test1[nthrow,2]-df_test1[nthrow,7]) < threshold || error("True and symbolic values do not match") #checking if the true value and symbolic computation matches
            @test abs(df_test1[nthrow,4]-df_test1[nthrow,7]) < threshold || error("Numerical and symbolic values do not match") #checking if the numerical computation and symbolic computation matches
        end        

        #generate nreps macaluay and matlab files one a same topoplogy but randomized parameters
        for i in 1:nreps
            ik1=readTopologyrand("((C,A),(((G,H),(((E,F))#H2)#H1),((#H2,(B,D)),#H1)));")
            network_expectedCF(ik1,symbolic=true,savecsv=true,filename="temp-macaulay-matlab-$i",inheritancecorrelation=ih,macaulay=true,matlab=true)
        end
        for i in 1:(nreps-1)
            for j in 2:nreps
                @test filecmp("temp-macaulay-matlab-$i.m2.txt", "temp-macaulay-matlab-$j.m2.txt")==true || error("Macaulay2 files do not match") 
                @test filecmp("temp-macaulay-matlab-$i.matlab.txt", "temp-macaulay-matlab-$j.matlab.txt")==true || error("Matlab files do not match") 
            end
        end
    end

    foreach(rm, filter(startswith("temp"), readdir()))#remove all temp files