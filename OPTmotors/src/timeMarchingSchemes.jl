include("../src/batchDiffTools.jl")

# I kind of use what Tibo did in ExactSolutions.jl

function timeMarchingScheme(opt, Nt, Δnum;sourceType="Ricker",t₀=50,f₀=0.03,initialCondition=0.0)

    #region reading data and source time functions

    operators = opt["numOperators"]
    costfunctions,fieldLHS,fieldRHS,champsLimité = operators
    costfunctions=reduce(vcat,costfunctions)
    #@show size(costfunctions),size(fieldLHS[1,1]),size(fieldRHS)

    timePointsUsedForOneStep = size(fieldLHS)[2]

    itVec=collect(1:1:Nt)
    t=(itVec.-1).*Δnum[end] # time vector # if it's not time marching t will give you just 0.0 regardless of Δnum[end]
    sourceTime = nothing
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

    #region 1Dvectorisation of knowns and unknowns

    #@show champsLimité[1,1], fieldLHS[1:end,1:end-1][1:end]

    symbKnownField = reduce(vcat,fieldLHS[1:end,1:end-1][1:end])
    symbUnknownField = reduce(vcat,fieldLHS[1:end,end][1:end])

    symbKnownForce = nothing
    if champsLimité === nothing
        symbKnownForce = reduce(vcat,fieldRHS[1:end,1:end][1:end])
    else
        symbKnownForce = reduce(vcat,champsLimité[1:end,1:end][1:end])
    end
   
 
    knownField = similar(symbKnownField)
    unknownField = copy(symbUnknownField)
    knownForce = similar(symbKnownForce)
    
     #endregion

    
    #region

    knownField .= initialCondition
    knownForce .= initialCondition

    sparseColouring(costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
    @show symbKnownForce


    # source time function will be shifted with timePointsUsedForOneStep - 1

    prepend!(sourceTime,zeros(timePointsUsedForOneStep))

    for it in itVec
        knownForce[1:timePointsUsedForOneStep] = sourceTime[it:it+timePointsUsedForOneStep-1]
        
    end


    #endregion



end