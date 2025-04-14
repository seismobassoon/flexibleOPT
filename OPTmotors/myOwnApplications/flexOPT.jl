# New version as of March 2025 for OPT operators
# Nobuaki Fuji @ IPGP/UPC/IUF
using  Pkg
Pkg.activate("../..")
#Pkg.add("DrWatson") 
#@quickactivate "flexibleDSM"


cd(Base.source_dir())       
              # active the project, with a  static environment
# Pkg.activate(; temp=true)    #  activate the project with a temporary environment
#Pkg.update()     

include("../src/imageReader.jl") # read 2D images for models
include("../src/batchNewSymbolics.jl")
include("../src/OPTnewEngines.jl") 
include("../src/famousEquations.jl")


# important!!! You can call the coordinates as you like but if you want to make a timeMarching, then
# the coordinate should be "t" and it should be a the last coordinate

using BenchmarkTools


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

    model=read2DimageModel(imagefile,colormap;Nwidth=21,Nheight=41,showRecoveredImage=false)

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


#DrWatson configurations

orderBtime=1
orderBspace=1
pointsInSpace=2
pointsInTime=2

#operatorConfigurations = @strdict famousEquationType Δnum orderBtime orderBspace pointsInSpace pointsInTime IneedExternalSources



#endregion



#region model configuration

concreteModelParameters = @strdict famousEquationType Δnum Nt orderBtime orderBspace pointsInSpace pointsInTime IneedExternalSources modelName models

# I know, I know, Nt is not very important to generate the numerical operators: I need to put Nt somewhere else ...

#endregion




#region Main programme

#region OPT symbolic derivation of objective functions to be minimised



f,file=produce_or_load(makeCompleteCostFunctions,concreteModelParameters,datadir("numOperators"))





#region je râle, je râle
#
# if one wants to work with different Δvalues, OPTobj should be called each time
# this is a bit of frustration for me due to the incapability of Symbolics to deal with stupid simplification 
#
#endregion




#constructingEquations(AjiννᶜU,)

#@show model

#@show AjiννᶜU


#endregion




#endregion