# here I test the famousGT98

# Nobuaki Fuji (April 2025)
#
# Here in this project, I will try to benchmark OPT operators with different configurations for a 1D SH wave equation

using  Pkg, BenchmarkTools, Symbolics, UnPack,Symbolics,JLD2

cd(@__DIR__)
Pkg.activate("../..")


myInclude("../src/OPTwrappers.jl") 


# Choose backend depending on environment
if isdefined(Main, :IJulia) || get(ENV, "JULIA_EDITOR", "") == "code"
    # Use CairoMakie for inline if in notebook or VS Code
    using CairoMakie
else
    using GLMakie
end

famousEquationType="1DsismoTime"
exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂²=famousEquations(famousEquationType)

if length(coordinates)  !== 2
    @error "with this benchmark, we just look for 1D in space and 1D in time heterogeneous problems"
end

logsOfΔxinverse = [0.5*i for i in 0:6]
logsOfΔtinverse = [0.5*i for i in 0:6]

numPointsX = collect(2:3)
numPointsT = collect(2:3)

k₀ = 1.0 # wave number for the wavefield
ω₀ = 1.0 # wave frequency for the wavefield

kᵣ = 1.0 # wave number for the density
dᵣ = 0.0 # shift in space
kₘ = 1.0 # wave number for the rigidity
dₘ = 0.0 # shif in space

cases=[]
cases = push!(cases,(name=famousEquationType*"samewavelength",u=cos(k₀*x-ω₀*t),ρ=3.0+sin(kᵣ*x+dᵣ),μ=3.0+sin(kₘ*x+dₘ)))

L = 10.0*π # the length of the segment