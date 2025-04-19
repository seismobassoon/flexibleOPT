include("../src/batchDiffTools.jl")

# I kind of use what Tibo did in ExactSolutions.jl

function timeMarchingScheme(opt, Nt, Δnum;sourceType="Ricker",t₀=50,f₀=0.03,initialCondition=0.0)

    #region reading data and source time functions

    operators = opt["numOperators"]
    costfunctions,fieldLHS,fieldRHS,champsLimité = operators
    #@show size(costfunctions),size(fieldLHS[1,1]),size(fieldRHS)

    itVec=collect(1:1:Nt)
    t=(itVec.-1).*Δnum[end] # time vector # if it's not time marching t will give you just 0.0 regardless of Δnum[end]

    #@show costfunctions[1,431] # check if we can see the source terms

    #specificication of parameters such as Nt (which can be 1 for no time-marching scheme) and source time function 


    if sourceType === "Ricker"
        # making a small source

        # if you want to plot:
        # myRicker(x)=Ricker(x,20,0.03)
        # scatter(0:1:100,myRicker) or lines 

       
        myRicker(x)=Ricker(x,t₀,f₀)
        
        sourceTime = myRicker.(t)
    end

    #endregion

    #region initial condition, 1Dvectorisation of knowns and unknowns

    #@show champsLimité[1,1], fieldLHS[1:end,1:end-1][1:end]

    knownField = fieldLHS[1:end,1:end-1][1:end]
    unknownField = fieldLHS[1:end,end][1:end]
    knownSource = similar(fieldRHS)
    @show costfunctions[1:end]

    


    #endregion

    
    #region

    for it in itVec

    end


    #endregion



end