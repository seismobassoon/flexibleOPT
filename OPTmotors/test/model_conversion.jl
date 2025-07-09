# New version from March 2025 for OPT operators
# Nobuaki Fuji @ IPGP/UPC/IUF
using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")

ParamFile = "../test/testparam_Moon.csv"

include("../src/DSM1D.jl")
using .DSM1D



writeClassicDSM1DPSVmodel(DSM1D.my1DDSMmodel,"moon.inf")

arrayRadius, arrayParams=DSM1D.compute1DseismicParamtersFromPolynomialCoefficients(DSM1D.my1DDSMmodel,10)
f=Figure()
#lines(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=1,color=:red)
lines(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ,color=:red)
scatter!(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=3,color=:blue)
display(f)




display(DSM1D.input.ωᵣ)
@testset "First series of tests" verbose=true begin
    #@test DSM1D.greet() === nothing
    @test DSM1D.dsm1Dconfig.re === 0.01
    #@test DSM1D.dsm1Dconfig.omegai === 0.01
    #println("DSM1D.input.omegai: ", DSM1D.input.omegai)
end