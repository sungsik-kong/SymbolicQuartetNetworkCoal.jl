        
    #macaulay output
    numCFs=size(df)[1]
    if(symbolic) 
        dataframe=deepcopy(df)
        params=gettingSymbolicInput(net, dataframe, inheritancecorrelation) 
    end


    if(macaulay) symbolic || error("symbolic must be set to true.") end
    str*="Write Macaulay2 file: "
    str*=(macaulay ? "on\n" : "off\n")

    if(macaulay)
        open("$filename.m2", "w") do file
        str="R = QQ["
        for par in params str=str*par*"," end
        str=str*"C_1..C_$numCFs]\n"
        str=str*"I = ideal(\n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2])-C_$i,\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])-C_$numCFs);\n"
        str=str*"G = eliminate (I, {"
        numparams=length(params)   
        for par in 1:(length(params)-1) str=str*params[par]*"," end
        str=str*"$(params[numparams])})\n"
        str=str*"S = QQ[C_1..C_$numCFs]\nJ = sub(G,S)\ndim J"
        write(file, str)
        end
    end

    #matlab output
    if(matlab) symbolic || error("symbolic must be set to true.") end
    str*="Write Matlab file: "
    str*=(matlab ? "on\n" : "off\n")

    if(matlab)
        open("$filename.m", "w") do file 
            str="% Declare variables\n"
            str=str*"syms "
            for par in params str=str*par*" " end
            for i in 1:numCFs str=str*"C_$i " end
            str=str*"\n\n% matrix of generating polynomials\n"
            str=str*"F=["
            i=1
            while i<numCFs
                str=str*"$(dataframe[i,2])-C_$i,\n"
                i+=1
            end
            str=str*"$(dataframe[numCFs,2])-C_$numCFs];\n"
            str=str*"\n% matrix of generating polynomials\n"
            str=str*"\n% Array of all variables\n"
            str=str*"V=["
            for par in params str=str*par*" " end
            for i in 1:numCFs-1 str=str*"C_$i " end
            str=str*"C_$numCFs]\n"
            str=str*"\n% Compute dimension\nCoalDim(F,V)"
        write(file, str)
        end
    end

    if(multigraded)
        open("$filename.im.m2.txt", "w") do file
        str="needsPackage \"MultigradedImplicitization\"\n"
        str*="R = QQ["
        for par in params str=str*par*"," end
        str*="T]\n"
        str*="S = QQ["
        str=str*"C_0..C_$numCFs]\n"
        str=str*"im = {\n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2]),\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])};\n"

        str=str*"im = {T} | apply(im, f -> f * T);\n"
        str=str*"phi = map(R, S, im)\n"
        str=str*"d=1\n"
        str=str*"I = time componentsOfKernel(d, phi)\n"
        str=str*"L = flatten values I"

        write(file, str)
        end
    end

    if(singular)
        open("$filename.sing", "w") do file
        str="ring R = 0, ("
        for par in params str=str*par*"," end
        for i in 1:numCFs str=str*"C$i," end
        str=chop(str)
        str=str*"), dp;\n"
        str=str*"ideal I = \n"
        i=1
        while i<numCFs
            str=str*"$(dataframe[i,2])-C$i,\n"
            i+=1
        end
        str=str*"$(dataframe[numCFs,2])-C$numCFs;\n"

        str=str*"ideal G = eliminate (I, "
        for par in params str=str*par*"*" end
        str=chop(str)
        str=str*");\n"

        str=str*"ring S = 0, ("
        for i in 1:numCFs str=str*"C$i," end
        str=chop(str)
        str=str*"), dp;\n"
        str=str*"ideal J = imap(R, G);\n"
        str=str*"int dimJ = dim(J);\n"
        str=str*"dimJ;"

        write(file, str)
    end
end

    if(savecsv) CSV.write("$filename.csv",df,header=false) end
    write(logfile,str)
    flush(logfile)