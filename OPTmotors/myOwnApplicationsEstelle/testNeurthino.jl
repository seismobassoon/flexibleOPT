# we are trying to incorporate Neurthino inside flexOPT.jl since the original package is no more maintained
#
# (https://github.com/KM3NeT/Neurthino.jl) but still we have to acknowledge this repository when we use this

using  Pkg, BenchmarkTools, Graphs

cd(@__DIR__)
Pkg.activate("../..")

ParamFile = "../test/testparam.csv"


include("../src/DSM1D.jl")
include("../src_Neurthino/Neurthino.jl")

#using .Neurthino: OscillationParameters, setθ!, setδ!, setΔm²!
using .Neurthino
using .DSM1D


osc = OscillationParameters(3)

# mixing angles

setθ!(osc, 1=>2, 0.59)
setθ!(osc, 1=>3, 0.15)
setθ!(osc, 2=>3, 0.84)

setδ!(osc, 1=>3, 3.86)

setΔm²!(osc, 2=>3, -2.523e-3)
setΔm²!(osc, 1=>2, -7.39e-5)

p = Pνν(osc, 1, 10000)
@show p[Energy=1, Baseline=1, InitFlav=Muon, FinalFlav=Tau]

U = PMNSMatrix(osc)
H = Hamiltonian(osc)
@show Pνν(U, H, 1, 10000)


# here every zenith angle can have a different nₑ profile (note that nₑ= ρ * "Z/A" )

zenith = acos.(range(-1,stop=0,length=200))
@show paths = Neurthino.prempath(zenith, 2.5, samples=100, discrete_densities=0:0.1:14);