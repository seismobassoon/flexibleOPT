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

for numPointsX_used in numPointsX
    for case in cases
        @unpack name,u,β = case
        q = -mySimplify(β*∂[1](u))
        @show qₓ = mySimplify(∂[1](q))
        for tmpLogOfHinverse in logsOfHinverse
            
            Δx = 1.0/exp(tmpLogOfHinverse)
            Δnum = (Δx)



        end
    end
end
