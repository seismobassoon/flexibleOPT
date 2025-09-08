module DSM1D
# This module will include all the functions needed to run the 1D DSM
# Packages that should be installed before running the DSM:

# - SeisBase
# - Metal (or CUDA if we work for NVIDIA GPUs)
# - DelimitedFiles
# - Test

# external packages

#using Metal # This should be replaced by CUDA if we work for NVIDIA GPUs 
            # https://github.com/SciML/DiffEqGPU.jl
#using SeisBase





#greet() = print("Hello World!")


# preprocessing 
include("mainStructures.jl")
export compute1DseismicParamtersFromPolynomialCoefficients
include("inputParams.jl")
include("inputModels.jl")

# if geodynamic model -> use another function and get parameter conversion
modelFile=dirname(@__FILE__)*"/"*DSM1D.dsm1Dconfig.modelFolder*"/"*DSM1D.input.modelFile
@time my1DDSMmodel=read1DModel(modelFile)

# for kernel computations, it is still interesting to make a regular grid in r and θ



# Here we make vertical grids in radius for the DSM computation
include("DSM1Dmotor.jl")


end # module DSM1D
