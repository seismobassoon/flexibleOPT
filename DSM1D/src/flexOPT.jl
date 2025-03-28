# New version as of March 2025 for OPT operators
# Nobuaki Fuji @ IPGP/UPC/IUF

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

modelDefinitionMethod="2DimageFile" # ToyModel or 2DimageFile (or 1DsphericalPlanet)
model2D =nothing
#endregion

#region Model input - option i) Model domain definition

if modelDefinitionMethod === "ToyModel"
    DomainWindow=(DomainWindowT=1.0,DomainWindowX=1.0,DomainWindowY=1.0,DomainWindowZ=1.0)
    ModelSizeTXYZ=(ModelSizeT=101,ModelSizeX=101,ModelSizeY=101,ModelSizeZ=0)
end

#endregion

#region Model input - option ii) Read a file (2D or 3D) and define Δs

if modelDefinitionMethod === "2DimageFile"

    #imagefile="DSM1D/data/model/random/colourful.jpg"
    #imagefile="DSM1D/data/model/artemis/IMG_6098.jpeg"
    #imagefile="DSM1D/data/model/random/tmp.png"
    imagefile = "DSM1D/data/model/random/marmousi.png"
    colormap = "jet" #colormap can be RGB vector or predefined colormap

    model=read2DimageModel(imagefile,colormap;Nwidth=201,Nheight=402,showRecoveredImage=false)

    # Either height or width should be decided 
    # (if the user precises the both, then they are respected)

    heightInMeter = 2000.0 # height in meter 
    widthInMeter = 3000.0 # width in meter 
    timeWindowInSecond = 200.0 # time window

    ΔheightInMeter = 20.0 # almost 
    ΔwidthInMeter = 20.0 # almost
    ΔtimeInSecond = 10.0 # almost

    #NF needs to make DomainWindow, ModelSizeTXYZ

    if heightInMeter !== nothing && ΔheightInMeter !== nothing
        ModelSizeX = trunc(Int, heightInMeter/ΔheightInMeter)+1
    else
        ModelSizeX = nothing
    end


    if widthInMeter !== nothing && ΔwidthInMeter !== nothing
        ModelSizeY = trunc(Int, widthInMeter/ΔwidthInMeter) +1
    else
        ModelSizeY = nothing
    end

    ModelSizeT = trunc(Int, timeWindowInSecond/ΔtimeInSecond)
    

end
#endregion

#region Model input - option iii) Read a file (1D spherical planet models)

if modelDefinitionMethod ==="1DsphericalPlanet"
    # use some programmes that are developed during Xmas 2023
    # inputModels.jl
end

#endregion



#region Main programme

#region OPT symbolic derivation of objective functions to be minimised

Δnum = (1.0,1.0,1.0) # this should be in the same order as coordinates 

exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = famousEquations(famousEquationType)
AjiννᶜU,dummyUtilities=OPTobj(exprs,fields,vars; coordinates=coordinates,CˡηSymbolicInversion=false,Δnum = Δnum)  
Γg,utilities=OPTobj(extexprs,extfields,extvars; coordinates=coordinates,CˡηSymbolicInversion=false,Δnum = Δnum)  

@show AjiννᶜU, Γg
#region je râle, je râle
#
# if one wants to work with different Δvalues, OPTobj should be called each time
# this is a bit of frustration for me due to the incapability of Symbolics to deal with stupid simplification 
#
#endregion

# here we need to give a numerical values 

#models = ((model.*0.5.+2), (1))
models=((model .* 0.5 .+ 2))
# if the dimension is degenerated, it is OK if the coordinate dependency is respected. The order will be taken based on the "coordinates" vector 


Nt= 120
modelPoints = (size(model)...,Nt) # Nx, Ny etc thing. Nt is also mentioned and it should be the last element!

constructingEquations(AjiννᶜU,Γg,coordinates,models,exprs,fields,vars,modelPoints,utilities;initialCondition=1.0)


#constructingEquations(AjiννᶜU,)

#@show model

#@show AjiννᶜU


#endregion




#endregion