include("../src/batchDiffTools.jl")
using DrWatson,JLD2
using GLMakie: Figure, Axis, heatmap!, Colorbar, record
# I kind of use what Tibo did in ExactSolutions.jl

function timeMarchingScheme(opt, Nt, Δnum;sourceType="Ricker",t₀=50,f₀=0.03,initialCondition=0.0)

    #region reading data and source time functions
    if !isdir(datadir("fieldResults"))
        mkdir(datadir("fieldResults"))
    end
    @show sequentialFileName=datadir("fieldResults", savename((Nt,Δnum...,sourceType),"jld2"))
    # filename that can be given by DrWatson
    @show compactFileName=datadir("fieldResults", savename("compact",(Nt,Δnum...,sourceType),"jld2"))


    operators = opt["numOperators"]
    costfunctions,fieldLHS,fieldRHS,champsLimité = operators
    costfunctions=reduce(vcat,costfunctions)
    #@show size(costfunctions),size(fieldLHS[1,1]),size(fieldRHS)

    timePointsUsedForOneStep = size(fieldLHS)[2]
    NField = size(fieldLHS)[1]
    @show pointsField = size(fieldLHS[1,1])

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
    
    symbKnownField = Array{Any,1}(undef,timePointsUsedForOneStep-1)
    for iT in 1:timePointsUsedForOneStep-1
        symbKnownField[iT] = reduce(vcat,fieldLHS[1:end,iT][1:end])
    end
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
    knownField .= initialCondition
    knownForce .= initialCondition
    #endregion

    #region sparse matrix colouring 

    J,colors=sparseColouring(costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
 
    #endregion

    #region time marching

    # source time function will be shifted with timePointsUsedForOneStep - 1

    prepend!(sourceTime,zeros(timePointsUsedForOneStep))

    F = zeros(length(costfunctions))
    unknownField .= initialCondition
    fieldFile = jldopen(sequentialFileName, "w") # file open


    
    for it in itVec
        knownForce[1:timePointsUsedForOneStep] = sourceTime[it:it+timePointsUsedForOneStep-1]
        #field to be shifted from the past
        

        timeStepOptimisation!(F, costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce,J,colors)

        # update
        knownField[1:end-1] = knownField[2:end]
        knownField[end] = unknownField

        # reshape
        @show newField = reshape(unknownField,NField,pointsField...)       
        fieldFile["timestep_$it"] = newField
    end
    close(fieldFile)

    #endregion


    #region making the file compact

    # Reconstruct array from separate datasets
    file = jldopen(sequentialFileName, "r")
    a = [file["timestep_$it"] for it in itVec]
    close(file)


    #endregion

    #region compact file
    # Convert to single 3D array: a[N, L, K]
    #a = permutedims(reshape(reduce(hcat, a), NField,pointsField...,Nt), (3, 1, 2))
    
    a=reduce(hcat,a)

    # Now save as compact single dataset
    @save compactFileName a compress=true

    #endregion

end