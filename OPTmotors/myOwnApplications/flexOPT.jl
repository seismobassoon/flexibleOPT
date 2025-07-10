# New version from March 2025 for OPT operators
# Nobuaki Fuji @ IPGP/UPC/IUF
using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")

include("../src/imageReader.jl") # read 2D images for models

include("../src/OPTwrappers.jl") 




# important!!! You can call the coordinates as you like but if you want to make a timeMarching, then
# the coordinate should be "t" and it should be a the last coordinate




#region Physics to choose, method to define material variable models
#famousEquationType="2DsismoTimeIsoHetero" # you can write a new governing equation using x-y-z-t coordinates 
#famousEquationType="1DsismoFreqHomo"
#famousEquationType="2DpoissonHetero"
famousEquationType="2DacousticTime"
#famousEquationType="1DsismoTime"

modelName="marmousi"

modelDefinitionMethod="2DimageFile" # ToyModel or 2DimageFile (or 1DsphericalPlanet)
model =nothing

#endregion

if modelDefinitionMethod !== nothing
#region Model input - option i) Model domain definition

if modelDefinitionMethod === "ToyModel"
    DomainWindow=(DomainWindowT=1.0,DomainWindowX=1.0,DomainWindowY=1.0,DomainWindowZ=1.0)
    ModelSizeTXYZ=(ModelSizeT=101,ModelSizeX=101,ModelSizeY=101,ModelSizeZ=0)
end

#endregion

#region Model input - option ii) Read a file (2D or 3D) and define Δs

if modelDefinitionMethod === "2DimageFile"

    #imagefile="../data/model/random/colourful.jpg"
    #imagefile="../data/model/artemis/IMG_6098.jpeg"
    #imagefile="../data/model/random/tmp.png"
    imagefile = "../data/model/random/marmousi.png"
    colormap = "jet" #colormap can be RGB vector or predefined colormap

    #model=read2DimageModel(imagefile,colormap;Nwidth=10,Nheight=10,showRecoveredImage=false)
    model=read2DimageModel(imagefile,colormap;showRecoveredImage=false)
end
#endregion

#region Model input - option iii) Read a file (1D spherical planet models)

if modelDefinitionMethod ==="1DsphericalPlanet"
    # use some programmes that are developed during Xmas 2023
    # inputModels.jl
end

#endregion
end


#region numerical configuration


Δnum = (1.0,1.0,1.0) # this should be in the same order as coordinates 



IneedExternalSources = true
maskedRegionForSourcesInSpace = nothing

#DrWatson configurations

orderBtime=1
orderBspace=1
pointsInSpace=2
pointsInTime=2

WorderBspace=1
WorderBtime=1
supplementaryOrder=2

#endregion


#region model configuration

# here we need to give a numerical values 

#models = ((model.*0.5.+2), (1))

models=[] # you might need to make this empty tuple first, otherwise one-member tuple can be misinterpreted
models=push!(models, (model .* 0.2 .+ 0.4))
# if the dimension is degenerated, it is OK if the coordinate dependency is respected. The order will be taken based on the "coordinates" vector 


# put fake Nt here for quasi-numerical operator construction
exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = famousEquations(famousEquationType)

fakeNt = 1
timeMarching = any(a -> a === timeDimensionString, string.(coordinates)) 
if timeMarching
    fakeNt = pointsInTime+1
    modelPoints = (size(model)...,fakeNt) # Nx, Ny etc thing. Nt is also mentioned and it should be the last element!
else
    modelPoints = (size(model))
end


# if IneedExternalSources and if the source region is localised in space then
maskedRegionForSourcesInSpace  = Array{CartesianIndex,1}(undef,0) # it is important to decalre the type of this
maskedRegionForSourcesInSpace = push!(maskedRegionForSourcesInSpace, CartesianIndex(modelPoints[1:end-1].÷2))# in Ndimension (or Ndimension  - 1 if timeMarching)
# in this example, I put a point source at the centre of the model space

forceModels =((1.0)) # if your model does not have anything special material parameters then it's how it's written

concreteModelParameters = @strdict famousEquationType Δnum orderBtime orderBspace WorderBtime WorderBspace supplementaryOrder pointsInSpace pointsInTime IneedExternalSources modelName models modelPoints forceModels maskedRegionForSourcesInSpace

#endregion



#region Main programme

#region OPT symbolic derivation of objective functions to be minimised, first semi-symbolically then numerically

#opt,file=@produce_or_load(makeCompleteCostFunctions,concreteModelParameters,datadir("numOperators");filename = config -> savename("quasiNum",concreteModelParameters))

opt = myProduceOrLoad(makeCompleteCostFunctions,concreteModelParameters,"numOperators","quasiNum")

#endregion

#region use the quasi-numerical operators to start computing
Nt=1

# for a non time marching scheme, put Nt=1 and  Δnum as a fake value (just put as it is)
timeMarchingScheme(opt, Nt, Δnum,modelName)




#endregion

#region je râle, je râle
#
# if one wants to work with different Δvalues, OPTobj should be called each time
# this is a bit of frustration for me due to the incapability of Symbolics to deal with stupid simplification 
#
#endregion







#endregion