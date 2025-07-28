

module Neurthino

using LinearAlgebra
using SparseArrays
using Graphs
#using StaticArrays
#using Polynomials
#using DocStringExtensions
#using LRUCache
#using LightGraphs
using AxisArrays

import Base

export oscprob, Pνν, OscillationParameters, PMNSMatrix, Hamiltonian, MatterOscillationMatrices
export masssquareddiff!, setΔm²!, cpphase!, setδ!, mixingangle!, setθ!
export cpphases, mixingangles

export NeutrinoFlavour, Electron, Muon, Tau

const N_A = 6.022e23 #[mol^-1]
const G_F = 8.961877245622253e-38 #[eV*cm^3]

# Julia 1.0 compatibility
isnothing(::Any) = false
isnothing(::Nothing) = true

include("Oscillation.jl")
include("Matter.jl")
#include("PREM.jl")

end # module
