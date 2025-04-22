# Nobuaki Fuji (April 2025)
#
# Here in this project, I will try to benchmark OPT operators with different configurations for a Poisson 1D problem

using  Pkg, BenchmarkTools, Symbolics, UnPack

cd(@__DIR__)
Pkg.activate("../..")


include("../src/OPTwrappers.jl") 


famousEquationType="1DpoissonHetero"
exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂²=famousEquations(famousEquationType)

if length(coordinates)  !== 1
    @error "with this benchmark, we just look for 1D heterogeneous problems"
end


modelName="1D_for_Poisson"

logsOfHinverse = [0.2*i for i in 0:5]

numPointsX = [2:4]

cases=[]

# manufactured ExactSolutions 

cases = push!(cases,(name="samewavelength",u=cos(x),β=sin(x)+2))
cases = push!(cases,(name="halfwavelength",u=cos(x),β=sin(x/2) + 2))
cases = push!(cases,(name="samewavelength_shifted_thirdpi",u=cos(x),β=sin(x+π/3) + 2))
cases = push!(cases,(name="twicewavelength",u=cos(x),β=cos(x).^2 + 1))
cases = push!(cases,(name="parabols",u=cos(x),β=x^2 .+ 1))
cases = push!(cases,(name="homogeneous",u=cos(x),β=1.0))
#

L = 10*π # the length of the segment

misfit = Array{Float64,3}(undef,length(logsOfHinverse),length(cases),length(numPointsX))

for iPointsUsed in eachindex(numPointsX)
    for iCase in eachindex(cases)
        @unpack name,u,β = cases[iCase]
        q = -mySimplify(β*∂[1](u))
        @show qₓ = mySimplify(∂[1](q))
        for iH in eachindex(logsOfHinverse)
            
            Δx = 1.0/exp(logsOfHinverse[iH])
            Δnum = (Δx)

            Nx = Int(L÷Δx) + 2

            @show X = [Δx * (i-1) for i ∈ range(1,Nx)]

            @show model=[substitute(β,Dict(x=>X[i])) for i ∈ range(1,Nx)]
            models=[]
            models=push!(models, model)
            

            IneedExternalSources = false
            maskedRegionForSourcesInSpace = nothing

            #DrWatson configurations

            orderBtime=1
            orderBspace=1
            pointsInSpace=numPointsX[iPointsUsed]
            pointsInTime=0

            modelName = name*string(Nx)

            modelPoints = (Nx)
            maskedRegionForSourcesInSpace =nothing
            
            concreteModelParameters = @strdict famousEquationType Δnum orderBtime orderBspace pointsInSpace pointsInTime IneedExternalSources modelName models modelPoints maskedRegionForSourcesInSpace


            opt,file=@produce_or_load(makeCompleteCostFunctions,concreteModelParameters,datadir("numOperators");filename = config -> savename("quasiNum",concreteModelParameters))
            syntheticData=timeMarchingScheme(opt, Nt, Δnum,modelName;videoMode=false)

            anayliticalData = [substitute(qₓ,Dict(x=>X[i])) for i ∈ range(1,Nx)]

            misfit[iH,iCase,iPointsUsed] = norm(syntheticData-anayliticalData)



        end
    end
end


@show misfit