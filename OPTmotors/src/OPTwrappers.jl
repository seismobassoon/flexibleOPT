using Symbolics,DrWatson

include("../src/famousSourceFunctions.jl")

function makeCompleteCostFunctions(concreteModelParameters::Dict)
    # This is a kind of big wrapper to construct an explicit numerical Cost functions to be minimised during the simulation
    

    #operators=wload(datadir("semiSymbolics", savename(operatorConfigurations,"jld2")))
    

    @unpack famousEquationType, Δnum, orderBtime, orderBspace, pointsInSpace, pointsInTime, IneedExternalSources, modelName, models, modelPoints, maskedRegionForSourcesInSpace = concreteModelParameters
    exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = famousEquations(famousEquationType)
    global ∂,∂²
    
    # here we construct semi symbolic operators (with numerical Δnum)
    operatorConfigurations = @strdict famousEquationType Δnum orderBtime orderBspace pointsInSpace pointsInTime IneedExternalSources
    operators,file=produce_or_load(OPTobj, operatorConfigurations, datadir("semiSymbolics"))


    # constructing numerical operator (with still symbolic expression for time coordinates)

    costfunctions,fieldLHS,fieldRHS = quasiNumericalOperatorConstruction(operators,modelName,models,famousEquationType,modelPoints,IneedExternalSources;maskedRegionForSourcesInSpace=maskedRegionForSourcesInSpace) 
    

    # 

    numOperators=(costfunctions=costfunctions,fieldLHS=fieldLHS,fieldRHS=fieldRHS)
    
    return @strdict(numOperators)
end

function quasiNumericalOperatorConstruction(operators,modelName,models,famousEquationType,modelPoints,IneedExternalSources;maskedRegionForFieldInSpace = nothing,maskedRegionForSourcesInSpace=nothing)

    # this is a big wrapper that reads the semi symbolic expressions to give a set of numerical operators (with symbolic expression in time)
    # which will call wrappers of onstructingNumericalDiscretisedEquations(Masked)

    operators=operators["operators"]
    operatorPDE,operatorForce,eqInfo = operators
    exprs,fields,vars,extexprs,extfields,extvars,coordinates=eqInfo

    AjiννᶜU,utilities=operatorPDE
    Γg,utilitiesForce=operatorForce

    lhsConfigurations = @strdict semiSymbolicOpt=AjiννᶜU coordinates modelName models fields vars famousEquationType modelPoints utilities maskedRegion=maskedRegionForFieldInSpace 

    numOperators,file = @produce_or_load(constructingNumericalDiscretisedEquations,lhsConfigurations,datadir("numOperators",savename(concreteModelParameters));filename = config -> savename(concreteModelParameters; ignores=["vars", "fields"]))


    # left-hand side, which is far more recyclable than r.h.s.
    
    costfunctionsLHS,fieldLHS=numOperators["numOperators"]
    

    costfunctionsRHS = similar(costfunctionsLHS)
    costfunctionsRHS .= 0.
    fieldRHS = similar(fieldLHS)
    fieldRHS .= 0.
    

    if IneedExternalSources 
        rhsConfigurations = @strdict semiSymbolicOpt=Γg coordinates modelName models=((1.0)) fields=extfields vars=extvars famousEquationType modelPoints utilities=utilitiesForce maskedRegion=maskedRegionForSourcesInSpace 
        numOperators,file=produce_or_load(constructingNumericalDiscretisedEquations,rhsConfigurations,datadir("numOperators",savename(concreteModelParameters));filename = config -> savename("source",concreteModelParameters; ignores=["vars", "fields"]))
       
        costfunctionsRHS,fieldRHS=numOperators["numOperators"]

    end
  
    #@show size(costfunctionsLHS),size(costfunctionsRHS)
    #@show costfunctionsRHS[1,430],costfunctionsRHS[1,431],costfunctionsRHS[1,434]
    #costfunctions=0.#costfunctionsLHS[1,1]-costfunctionsLHS[1,1]
    costfunctions = costfunctionsLHS .- costfunctionsRHS
    #numOperators=(costfunctions=costfunctions)
    return costfunctions,fieldLHS,fieldRHS
end
