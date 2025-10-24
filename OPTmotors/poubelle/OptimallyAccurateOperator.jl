myInclude("../src/OPTEngines.jl")
myInclude("../src/batchSymbolics.jl")
myInclude("../src/batchUseful.jl")

myInclude("../src/famousEquations.jl")
myInclude("../src/readOrMakeModels.jl")
myInclude("../src/imageReader.jl") # we can read 2D photos with a user-defined colorbar
myInclude("../src/plotTomography.jl") # myPcolor is inside (very useful for plot just a matrix)
using BenchmarkTools

# available equations: 1Dacceleration,1Dlaplacian,1DsismoFreq,2DsismoSHFreq,
#                      1DsismoTime,2DpoissonHomo,2DpoissonHetero,
#                      1DpoissonHetero,2DsismoTimeIso,highSchoolProblem
#                      3DsismoTimeIso is not yet working


#region Input parameters

#region 1. Physics to choose, method to define material variable models

famousEquationType="2DsismoTimeIsoHomo"
modelDefinitionMethod="ToyModel" # ToyModel or 2DimageFile (or 1DsphericalPlanet)

#endregion

#region 2.

#region 2 - option i) Model domain definition

if modelDefinitionMethod === "ToyModel"
    DomainWindow=(DomainWindowT=1.0,DomainWindowX=1.0,DomainWindowY=1.0,DomainWindowZ=1.0)
    ModelSizeTXYZ=(ModelSizeT=101,ModelSizeX=101,ModelSizeY=101,ModelSizeZ=0)
end

#endregion

#region 2 - option ii) Read a file (2D or 3D) and define Δs

if modelDefinitionMethod === "2DimageFile"

    showRecoveredImage=false # this is an option for read2DimageModel

    #imagefile="DSM1D/data/model/random/colourful.jpg"
    imagefile="DSM1D/data/model/artemis/IMG_6098.jpeg"
    #imagefile="DSM1D/data/model/random/tmp.png"
    #imagefile = "DSM1D/data/model/random/marmousi.png"
    colormap = "jet" #colormap can be RGB vector or predefined colormap

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

#region 2 - option iii) Read a file (1D spherical planet models)

if modelDefinitionMethod ==="1DsphericalPlanet"
    # use some programmes that are developed during Xmas 2023
    # inputModels.jl
end

#endregion

#endregion

#endregion




#region Main programme

#region 1. OPT symbolic derivation of Au-g
@show exprs,fields,vars,extexprs,extfields,extvars = famousEquations(famousEquationType)
@time Au,g,neighbourVars,pointsUsedLeftAndRight,Δs=CartesianOPTSymbolics(exprs,fields,vars,extexprs,extfields,extvars; iDoSTEP7=true)  
@show Au-g
#endregion

#region 2. option i) making discretised model with "1Dto3DcocentricSinePerturbation" option in makeHeterogeneousModels function

if modelDefinitionMethod === "ToyModel"
    # here we give grids for t, x, y, z
    Δnodes = regularGridConstruction(DomainWindow,ModelSizeTXYZ)  

    materials=readOrMakeDiscretisedModel(vars)
    modellingVariables=(materials=materials,Δnodes=Δnodes)
end

#endregion

#region 2. option ii) Read a file (2D or 3D) and define Δs

if modelDefinitionMethod === "2DimageFile"
    read2DimageModel(imagefile,colormap; showRecoveredImage=true) #colormap can be RGB vector or predefined colormap
end

#endregion

#region 2. option iii) Read a file (1D spherical planet models)

if modelDefinitionMethod ==="1DsphericalPlanet"
end

#endregion

#region 3. 

#region 4. make Au-g vectors numerically to be minimised

# here we concretise the Au-g numerically
makeNumericalAuMoinsG((Au-g),fields,vars,neighbourVars,pointsUsedLeftAndRight,ModelSizeTXYZ,modellingVariables,Δs)
#endregion

#region 5. fastdiff

#endregion

#endregion










#region some other ideas, old options that can help to improve the current workflow, etc.

# iDoSTEP7 should be false! "true" option is only to verify that our coefficients 
# are consistent with those from GT95/98 and others (NF Dec 2024)


# if we chose iDoSTEP7=true, then
# nouvelleExpr,neighbourVars,pointsUsedLeftAndRight,Δs =CartesianOPTSymbolics(exprs,fields,vars,extexprs,extfields,extvars)  

# here we read or make the model of material variables
# be careful, here each material can have <4 dimensions but it should be in the order of X-Y-Z then T (if ever T dependency is included in the material)


# Below is the old version and I think it is not relevant to keep

#@time A₀,δA,futureOnepoint,ModelSizeTXYZ,nodeVarsModel=ConcretisingAmatrix(Au,0,neighbourVars,pointsUsedLeftAndRight,fields,vars,ModelSizeTXYZ)

#@time makeNumericalAmatrixExplicit(A₀,δA,futureOnepoint,ModelSizeTXYZ,nodeVarsModel,modelingVariables,Δs)

#endregion


# to do: i) high-order B-spline; ii) external source; iii) boundary conditions; iv) irregular grids

