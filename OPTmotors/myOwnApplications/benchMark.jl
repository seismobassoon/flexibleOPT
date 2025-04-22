# Nobuaki Fuji (April 2025)
#
# Here in this project, I will try to benchmark OPT operators with different configurations for a Poisson 1D problem

using  Pkg, BenchmarkTools, Plots, Symbolics, UnPack

cd(@__DIR__)
Pkg.activate("../..")

include("../src/imageReader.jl") # read 2D images for models

include("../src/OPTwrappers.jl") 


famousEquationType="1DpoissonHetero"
exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂²=famousEquations(famousEquationType)

modelName="1D_for_Poisson"

logsOfHinverse = [0.2*i for i in 0:5]

cases=[]

# manufactured ExactSolutions (Cyrillic = analytical solution)
У(x) = cos(x)
Б(x) = sin(x) +2
cases = push!(cases,(У,Б))

#


for case in cases
    @unpack У,Б = case
    for tmpLogOfHinverse in logsOfHinverse
        
        Δx = 1.0/exp(tmpLogOfHinverse)
        Δnum = (Δx)

        

    end
end

