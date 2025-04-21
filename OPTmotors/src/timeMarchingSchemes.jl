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
    costfunctions=reduce(vcat,costfunctions) # we know that costfunctions is a 2D array

    #@show size(costfunctions),size(fieldLHS[1,1]),size(fieldRHS)

    timePointsUsedForOneStep = size(fieldLHS)[2]
    NField = size(fieldLHS)[1]

    pointsFieldSpace = size(fieldLHS[1,1]) # the list of points in space which is useful for reshaping
    NpointsSpace = length(fieldLHS[1,1]) # the full number of points

    itVec=collect(1:1:Nt)
    t=(itVec.-1).*Δnum[end] # time vector # if it's not time marching t will give you just 0.0 regardless of Δnum[end]
    sourceTime = nothing
    #@show costfunctions[1,431] # check if we can see the source terms


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
    
    symbKnownField = Array{Num,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep-1)
    knownField = Array{Number,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep-1)
    Rcoord = CartesianIndices(pointsFieldSpace)
    for iT in 1:timePointsUsedForOneStep-1
        for iField in 1:NField
            for j in Rcoord
                linearJ = LinearIndices(Rcoord)[j]
                symbKnownField[linearJ,iField,iT]=fieldLHS[iField,iT][j]
                knownField[linearJ,iField,iT] = initialCondition
            end
        end
    end

    symbUnknownField = Array{Num,2}(undef,NpointsSpace,NField)
    unknownField = Array{Float64,2}(undef,NpointsSpace,NField)
    for iField in 1:NField
        for j in Rcoord
            linearJ = LinearIndices(Rcoord)[j]
            symbUnknownField[linearJ,iField] = fieldLHS[iField,end][j]
            unknownField[linearJ,iField] = initialCondition
        end
    end

    symbKnownForce = nothing
    if champsLimité === nothing
        symbKnownForce=Array{Num,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                for j in Rcoord
                    linearJ = LinearIndices(Rcoord)[j]
                    symbKnownField[linearJ,iField,iT]=fieldLHS[iField,iT][j]
                    knownField[linearJ,iField,iT] = initialCondition
                end
            end
        end
    else
        NpointsLimités = length(champsLimité[1,1])
        symbKnownForce=Array{Num,3}(undef,NpointsLimités,NField,timePointsUsedForOneStep)
        knownForce=Array{Number,3}(undef,NpointsLimités,NField,timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                for j in eachindex(NpointsLimités)
                    symbKnownForce[j,iField,iT]=champsLimité[iField,iT][j]
                    knownForce[j,iField,iT] = initialCondition
                end
            end
        end
    end
   
    #endregion

    #region numerical operators

    f= buildNumericalFunctions(costfunctions,symbUnknownField,symbKnownField,symbKnownForce)

    #end


    #region sparse matrix colouring 

    J,cache=sparseColouring(f,unknownField,knownField,knownForce)
 
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
        
        # this is not true!!! Just for debugging!
        knownForce[1:timePointsUsedForOneStep] .= 1.0

        timeStepOptimisation!(f,unknownField,knownField,knownForce,J,cache,pointsFieldSpace)

        # update
        knownField[:,:,1:end-1] = knownField[:,:,2:end]
        knownField[:,:,end] = unknownField[:,:]

        # reshape
        newField = reshape(unknownField,NField,pointsFieldSpace...)       
        fieldFile["timestep_$it"] = newField
        scene = heatmap(Float32.(newField), colormap =  :deep,colorrange=(-1.e-5,1.e-5))
        display(scene)
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