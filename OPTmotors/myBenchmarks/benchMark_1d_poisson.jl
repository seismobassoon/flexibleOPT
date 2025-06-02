# Nobuaki Fuji (April-june 2025)
#
# Here in this project, I will try to benchmark OPT operators with different configurations for a Poisson 1D problem

using  Pkg, BenchmarkTools, Symbolics, UnPack,Symbolics,JLD2

cd(@__DIR__)
Pkg.activate("../..")


include("../src/OPTwrappers.jl") 


# Choose backend depending on environment
if isdefined(Main, :IJulia) || get(ENV, "JULIA_EDITOR", "") == "code"
    # Use CairoMakie for inline if in notebook or VS Code
    using CairoMakie
else
    using GLMakie
end


famousEquationType="1DpoissonHetero"
exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂²=famousEquations(famousEquationType)

if length(coordinates)  !== 1
    @error "with this benchmark, we just look for 1D heterogeneous problems"
end


#modelName="1D_for_Poisson"


logsOfHinverse = [1.0*i for i in 0:3]

numPointsX = collect(2:3)
tmpOrderBtime=1
tmpOrderBspace=-1

cases=[]

# manufactured ExactSolutions 
prefix="B"*string(tmpOrderBspace)*"_"
#cases = push!(cases,(name=prefix*"sameλ",u=cos(x),β=sin(x)+2))
#cases = push!(cases,(name=prefix*"twiceλ",u=cos(x),β=sin(x/2) + 2))
#cases = push!(cases,(name=prefix*"sameλ_shifted_π_3",u=cos(x),β=sin(x+π/3) + 2))
cases = push!(cases,(name=prefix*"λ_2",u=cos(x),β=cos(x).^2 + 1))
#cases = push!(cases,(name=prefix*"parabols",u=cos(x),β=x^2+ 1))
cases = push!(cases,(name=prefix*"homogeneous",u=cos(x),β=1.0))
#

L = 10.0*π # the length of the segment

misfit = Array{Float64,3}(undef,length(logsOfHinverse),length(cases),length(numPointsX))

fileMisfit="misfit_B"*string(tmpOrderBspace)*string(numPointsX[end])*".jld2"

if !isfile(fileMisfit)
for iPointsUsed in eachindex(numPointsX)
    for iCase in eachindex(cases)
        @unpack name,u,β = cases[iCase]
        q = mySimplify(β*∂[1](u))
        @show qₓ = mySimplify(∂[1](q))
        for iH in eachindex(logsOfHinverse)
            iExperiment = (iH=iH,iCase=iCase,iPointsUsed=iPointsUsed)
            ΔxTry = 1.0/exp(logsOfHinverse[iH])
            Nx = Int(L÷ΔxTry) +1
            Δx = L/(Nx-1)
        

            Δnum = (Δx)
            @show X = [Δx * (i-1) for i ∈ range(1,Nx)]

      
            model=[Symbolics.value(substitute(β,Dict(x=>X[i]))) for i ∈ range(1,Nx)]
            models=[]
            models=push!(models, model)
            

            IneedExternalSources = true
            maskedRegionForSourcesInSpace = nothing

            force = [Symbolics.value(substitute(qₓ,Dict(x=>X[i]))) for i ∈ range(1,Nx)]
            sourceFull=reshape(force,Nx,1,1)

            #DrWatson configurations

            orderBtime=tmpOrderBtime
            orderBspace=tmpOrderBspace
            pointsInSpace=numPointsX[iPointsUsed]
            pointsInTime=0

            modelName = name*string(Nx)

            modelPoints = (Nx)
            maskedRegionForSourcesInSpace =nothing
            

            forceModels =((1.0)) # if your model does not have anything special material parameters then it's how it's written

            concreteModelParameters = @strdict famousEquationType Δnum orderBtime orderBspace pointsInSpace pointsInTime IneedExternalSources modelName models modelPoints forceModels maskedRegionForSourcesInSpace iExperiment


            opt,file=@produce_or_load(makeCompleteCostFunctions,concreteModelParameters,datadir("numOperators");filename = config -> savename("quasiNum",concreteModelParameters))

            Nt = 1
            
            
            syntheticData=timeMarchingScheme(opt, Nt, Δnum,modelName;videoMode=false,sourceType="Explicit",sourceFull=sourceFull,iExperiment=iExperiment)
            syntheticData=reduce(vcat,syntheticData)
            analyticalData = [Symbolics.value(substitute(u,Dict(x=>X[i]))) for i ∈ range(1,Nx)]

            fig =Figure()
            ax=Axis(fig[1,1]; title="Comparison for h=$Δx, model=$(cases[iCase].name), points=1+$pointsInSpace")
            lines!(ax,X,analyticalData,color=:blue,label="analytical")
            scatter!(ax,X,syntheticData,color=:red,marker=:circle,label="synthetic")        
            axislegend(ax)
            display(fig)
            
            #save("plot_$iH_$iCase_$iPointsUsed.png",fig)
            sleep(0.5)
            misfit[iH,iCase,iPointsUsed] = norm(syntheticData-analyticalData)

            

        end
    end
end

@save fileMisfit misfit
end
@load fileMisfit misfit

@show misfit

fig =Figure()
ax=Axis(fig[1,1]; title="Misfit")
N=length(cases)
colors = [get(Makie.colorschemes[:viridis], (i - 1) / (N - 1)) for i in 1:N]
for iCase in eachindex(cases)
    scatter!(ax,logsOfHinverse,log.(misfit[:,iCase,1]),color=colors[iCase],label=cases[iCase].name)
end

#O_1=log.(misfit[1,1,1]).-1.0*logsOfHinverse
O_2=log.(misfit[1,1,1]).-2.0*logsOfHinverse
O_4=log.(misfit[1,1,1]).-4.0*logsOfHinverse
#O_8=log.(misfit[1,1,1]).-1.0*logsOfHinverse
#lines!(ax,logsOfHinverse,O_1,color=:black,label="O1")
lines!(ax,logsOfHinverse,O_2,color=:black,label="O2")
lines!(ax,logsOfHinverse,O_4,color=:black,label="O4")
#lines!(ax,logsOfHinverse,O_8,color=:black,label="O8")
#axislegend(ax,position=:lb)
display(fig)