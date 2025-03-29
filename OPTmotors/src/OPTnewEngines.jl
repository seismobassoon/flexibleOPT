using Symbolics,UnPack,LinearAlgebra

include("../src/batchNewSymbolics.jl")
include("../src/batchUseful.jl")

# PDECoefFinder cannot detect the material partials × material partials for the moment!! 

# in the near future, xyzt coordinates shouold be more flexible! With CartesianIndices!!

function findCartesianDependency(expression,Ndimension)
    
    expressionDependency=ones(Int,Ndimension)
    for iDimension in 1:Ndimension
        if typeof(expand_derivatives(∂[iDimension](expression))==0) <: Bool 
            if expand_derivatives(∂[iDimension](expression))==0
                expressionDependency[iDimension] = 0
            end
        end
    end
    return expressionDependency
end

function makeMixPartials(orders,coordinates;field=identity)
    # this function will give a matrix of mixed partial derivative operators 
    # coordinates should be an array of Symbolics variables 
    # orders is a matrix 


    Ndimension=length(coordinates)
    if length(orders)!== Ndimension
        @error "the highest orders array has not the same dimension as that of the coordinates"
    end

    ∂ = []
    for iDim in 1:Ndimension
        ∂ = push!(∂,Differential(coordinates[iDim]))
    end

    ∇ = Array{Any,Ndimension}(undef, Tuple(orders))
    R=CartesianIndices(∇)
        
    ∇ .= field
    for I in R
        for iDim in 1:Ndimension
            ∇[I] = (∂[iDim]^(I[iDim]-1))(∇[I])
        end
    end
   
    return ∇
end

function PDECoefFinder(pointsUsed,coordinates,expr,field,vars)
    # PDECoefFinder cannot detect the material partials × material partials for the moment!! 
    # I know how to do it, but eq. 40 should be then more generalised (kind of the product of partials of different materials)

    # maxPolynomialOrderMaterial is also a chelou thing, that I need to work on more systematically
    # like the powers of partials should also be included but here search for Rm[1], yeah, that's what I am doing

    Ndimension = length(coordinates)
    alpha=[]
    
    maxPolynomialOrderMaterial = 2*(maximum(pointsUsed)-1)
    ∇=makeMixPartials(pointsUsed,coordinates;field=field)
    R=CartesianIndices(∇)
    expr=mySimplify(expr)

    varM=Array{Any,2}(undef,length(vars),length(R))
   
    for iVar in eachindex(vars)


        newstring=split(string(vars[iVar]),"(")[1]
     
        
        CartesianDependency=findCartesianDependency(vars[iVar],Ndimension)
       
        smallVarM=Symbolics.variables(Symbol(newstring),1:length(R))
        for j in R
            linearJ=LinearIndices(R)[j]
            realJ=(car2vec(j).-1).*CartesianDependency .+1 # if there is no dependence on a direction, it should get the same name
            linearRealJ=LinearIndices(R)[CartesianIndex(realJ...)]
            smallVarM[linearJ]=smallVarM[linearRealJ]
        end
        varM[iVar,:]=smallVarM
    end
    


    for i in R
        term_searched = ∇[i]

        tmpCoeff = myCoeff(expr,term_searched)
        if tmpCoeff !== 0
            isTmpCoeffAConstant=true
            for iVar in eachindex(vars)
                
                ∇m=makeMixPartials(pointsUsed,coordinates;field=vars[iVar]) # material partials
                Rm=CartesianIndices(∇m)
                for j in Rm
                    term_material_searched = ∇m[j]
                    tmpCoeffMaterial = myCoeff(tmpCoeff,term_material_searched)
                    
                    if tmpCoeffMaterial !==0
                        isTmpCoeffAConstant=false
                        isOKtoinclude =true
                        # This is to avoid partials of other material 
                        for jVar in eachindex(vars)
                            if jVar !== iVar
                                ∇n = makeMixPartials(pointsUsed,coordinates;field=vars[jVar]) 
                                for jj in Rm
                                    if jj !== Rm[1]
                                        term_material_searched_plus = ∇n[jj]
                                        differentialCoeff = myCoeff(tmpCoeffMaterial,term_material_searched_plus)
                                        if differentialCoeff !==0
                                            isOKtoinclude = false
                                        end
                                    end
                                end
                            end
                        end
                        if isOKtoinclude
                            #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                            specificMaterialTerm=tmpCoeffMaterial*vars[iVar]
                            tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=j,n=i) # the famous n prime and n in the equation 56 or 40
                            alpha = push!(alpha,tmpAlphaIJ)
                        end
                    end
                end
                for matPower in 2:maxPolynomialOrderMaterial
                    tmpCoeffMaterial = myCoeff(tmpCoeff,vars[iVar]^matPower)
                    #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                    if tmpCoeffMaterial !==0
                        specificMaterialTerm=tmpCoeffMaterial*vars[iVar]^matPower
                        tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=Rm[1],n=i) # the famous n prime and n in the equation 56 or 40
                        alpha = push!(alpha,tmpAlphaIJ)
                    end
                end
            end
            if isTmpCoeffAConstant
                specificMaterialTerm = tmpCoeff
                tmpAlphaIJ=(node=specificMaterialTerm,nᶜ=R[1],n=i)
                alpha = push!(alpha,tmpAlphaIJ)
            end
        end
    end
    alpha=unique(alpha)


    return alpha, varM # varM: iVar and linearised cartesian indices
end 

function illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=true,Δ=nothing)

    # here we propose a big Taylor expansion matrix with Δcoordinates, symbolically (when Δ=nothing or Δ as a symbolic array)
    #   or numerically otherwise

    Ndimension=length(coordinates)

    if Δ === nothing
       Δ = Symbolics.variables(:Δ,1:Ndimension)
    else
        if length(Δ) !== Ndimension
            @error "the numerical delta increment has not the same dimension!"
        end
    end

    numberOfEtas = length(multiPointsIndices)
    numberOfLs   = length(multiOrdersIndices)

    #TaylorExpansionCoeffs=Array{Any,Ndimension*2}(undef, multiOrdersIndices[end],multiPointsIndices[end])
    #TaylorExpansionCoeffs=Array{Any,Ndimension*2}(undef, Tuple(vcat(collect(Tuple(multiOrdersIndices[end])),collect(Tuple(multiPointsIndices[end])))))
    
    CˡηGlobal = Array{Any,3}(undef,numberOfEtas,numberOfLs,numberOfEtas)
    midLinearK = nothing # this is valid only for testOnlyCentre

    for k in multiPointsIndices
        linearK = LinearIndices(multiPointsIndices)[k]
        TaylorExpansionCoeffs = Array{Num,2}(undef,numberOfLs,numberOfEtas)
        if !testOnlyCentre || Tuple(k) === ((Tuple(multiPointsIndices[end] )).-1 ).÷2 .+1 
            midLinearK = linearK # this is valid only for testOnlyCentre
            for i in multiPointsIndices
                linearI = LinearIndices(multiPointsIndices)[i]
                η = car2vec(i-k)
                distances= η .* Δ
                for j in multiOrdersIndices
                    linearJ = LinearIndices(multiOrdersIndices)[j]
                    orders = car2vec(j).-1
                    numerator = prod(distances .^orders)
                    denominator=prod(factorial.(orders))
                    tmpTaylorCoeffs = numerator/denominator
                    TaylorExpansionCoeffs[linearJ,linearI]=tmpTaylorCoeffs 

                end
            end
            # here we do the famous inversion (ttttttt) even though this code is essentially a forward problem
            
            aa=transpose(TaylorExpansionCoeffs)*TaylorExpansionCoeffs
            invaa= myInv(aa)
            CˡηGlobal[:,:,linearK]=invaa*transpose(TaylorExpansionCoeffs)
        end
    end 

    if testOnlyCentre
        CˡηCentre = CˡηGlobal[:,:,midLinearK]
        CˡηGlobal = nothing
        return CˡηCentre,Δ,multiOrdersIndices
    else
        #@show CˡηGlobal
        return CˡηGlobal,Δ,multiOrdersIndices
    end
end

function integralBsplineTaylorKernels1D(BsplineOrder,Δ,l_n_variable,l_n_field)
    
    # this will compute \int dx Bspline K_{l-n} K_{lᶜ-nᶜ}
    middle_value = 0
    extreme_value = 0
    if BsplineOrder=== -1
        # this is for a delta function
        if l_n_variable === 0 && l_n_field === 0
            middle_value=1
            extreme_value=1
        else
            middle_value=0
            extreme_value=0
        end

    elseif BsplineOrder === 0
        #
    elseif BsplineOrder=== 1
        middle_value = (Δ^(l_n_variable+l_n_field+1)-(-Δ)^(l_n_variable+l_n_field+1))/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #extreme_value = (Δ^{l_n_variable+l_n_field+1})/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(l_n_variable)*factorial(l_n_field))
    end
    return middle_value,extreme_value
end

function OPTobj(exprs,fields,vars; coordinates=(x,y,z,t), trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, highestOrderPartial=2),CˡηSymbolicInversion=false,testOnlyCentre=true,Δnum = nothing)

    #region General introduction, some cautions

    # Nobuaki Fuji @ IPGP/UPC/IUF since 2024
    #
    # 
    # encouraged by Thibault Duretz @ U. Frankfurt Goethe, Kurama Okubo @ NIED
    #
    # with some inputs from Giacomo Aloisi @ ETH during Julia hackathon 
    #                                  in the black forest October 2024
    #
    # intermediate presentations: IPGP-CIA workshop October 2024; IPGP-ERI workshop November 2024; lighttalk @ systemI December 2024
    #

    # This function SHOULD offer an OPT operator for a given expression 
    # and the boundary conditions on the 8 sides of a local 4D Cartesian box 

     # B-spline order (-1: Dirac's delta function (finite difference); 0: step function (2 points); 1: triangle (3 points); 2+: beautiful curve) 
    
    # this should be corrected according to the orderPDEs that should be given by the gvn. eqn.
    # orderBtime and orderBspace should be given


    # This function should expand the compact PDE given by a textbook or by an imagination
    # in order to organise an OPT operators for 4D Cartesian coordinates. 

    # orderBtime and orderBspace are the order of B-spline (if -1: delta function) to be used as test functions
    # as the paper in prep. (as of August 2024) shows, the numbers of points used (pointsUsed) are not necessarily related to orderBtime and orderBspace

    # even for spherical harmonic version with only z-dependency, this should work (but the expansion should be done a posteriori)

    # Note that variableDependency defines the (maximum) dimension of physics to be solved. 



    # some notes on March 2025

    # highestOrderPartial should be somehow determined more naturally!!
    # (l-n) max is the same as pointsUsedForFields here!

    # if timeMarching === true then we consider that the last coordinate is time
    
    # CˡηSymbolicInversion is highly recommended to be false since it takes really a big effort for nothing
    # pointsUsed(ForFields) should be the highest order of PDE + 1 at least




    #endregion

    #region initialising, unpacking etc. 

    timeMarching = any(a -> a === t, coordinates)


    @unpack orderBtime, orderBspace, highestOrderPartial = trialFunctionsCharacteristics

    NtypeofExpr=length(exprs)   # number of governing equations
    NtypeofMaterialVariables=length(vars) # number of material coefficients
    NtypeofFields=length(fields) # number of unknown fields
    
    Ndimension = length(coordinates) # we do not change this for the moment, especially for the time-marching scheme
    pointsUsed = ones(Int, Ndimension).*(highestOrderPartial+1)


    if length(Δnum) !== Ndimension && !CˡηSymbolicInversion
        @error "the numerical delta increment has not the same dimension!"
    end

 


    #endregion

    #region investigation of all the fields and vars dependencies in terms of x-y-z-t

    variableDependency=ones(Int,Ndimension)
    fieldDependency=ones(Int,Ndimension)
    eachVariableDependency=ones(Int,Ndimension,NtypeofMaterialVariables) 
    eachFieldDependency=ones(Int,Ndimension,NtypeofFields)
  
    for iFields in 1:NtypeofFields
        eachFieldDependency[:,iFields]=findCartesianDependency(fields[iFields],Ndimension)
        fieldDependency = fieldDependency .* (ones(Int,Ndimension).-eachFieldDependency[:,iFields])
    end


    for iVars in 1:NtypeofMaterialVariables
        eachVariableDependency[:,iVars]=findCartesianDependency(vars[iVars],Ndimension)
        variableDependency = variableDependency .* (ones(Int,Ndimension).-eachVariableDependency[:,iVars])
    end

    

    fieldDependency = ones(Int,Ndimension).-fieldDependency
    variableDependency = ones(Int,Ndimension).-variableDependency

    # here we correct variableDependency with fieldDependency: if fieldDependency is zero then we do not take care of that dimension for the variables
    variableDependency = variableDependency .* fieldDependency

    #endregion

    #region definition of points in time and space to be used

    # heaviside(x) = x > 0 ? 1 : x == 0 ? 0 : -1

    # the orders of B-spline functions, depending on fields 

    orderBspline=zeros(Int,Ndimension)

    if timeMarching
        orderBspline[Ndimension]=orderBtime*fieldDependency[Ndimension]
        orderBspline[1:Ndimension-1]=orderBspace*fieldDependency[1:Ndimension-1]
    else
        orderBspline[1:Ndimension]=orderBspace*fieldDependency[1:Ndimension]
    end
    
    # the number of points used in the vicinity of the node, which is independent of the order of B-spline functions (see our paper)
    pointsUsedForFields=(pointsUsed.-1).*fieldDependency.+1
    
    # numbers of points to evaluate the integral for the governing equation filtered by the test functions
    
    # orderU is the maximal orders for the fields that we will use for OPT coefficients' exploration
    orderU = (pointsUsedForFields.-1).*2 .+1
    # we restore this orderU since we need to control this (we set this to the twice the size of the number of used points)

    #endregion

    #region analysis of expressions to obtain the α_{n'nji}

    bigα=Array{Any,2}(missing,NtypeofFields,NtypeofExpr)
    varM=nothing
    for iExpr in eachindex(exprs)
        for iField in eachindex(fields)
            
            tmpNonZeroAlphas,varM=PDECoefFinder(pointsUsedForFields,coordinates,exprs[iExpr],fields[iField],vars) 
            # we assume that the pointsUsedForFields represent the highest order of partials
            bigα[iField,iExpr]=unique(tmpNonZeroAlphas)
        end
    end
    @show bigα
    #endregion

    #region Preparation for Taylor expansion
    
    orderTaylors=Array{Any,Ndimension}(undef,Tuple(orderU))
    pointsInSpaceTime=Array{Any,Ndimension}(undef,Tuple(pointsUsedForFields))
    
    # Cartesian indices that can be interesting to use

    multiOrdersIndices=CartesianIndices(orderTaylors)
    multiPointsIndices=CartesianIndices(pointsInSpaceTime)

    NpointsMultiDim = length(multiPointsIndices)
    NordersMultiDim = length(multiOrdersIndices)

    #endregion

    #region obtaining Cˡη either symbolically either with Δcoordinates in a numerical way

    if CˡηSymbolicInversion # this seems super cool but it takes time
        Cˡη,Δ,multiLCar = illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=testOnlyCentre)
    else
        Cˡη,Δ,multiLCar = illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=testOnlyCentre,Δ=Δnum)
        # this clause can work only if the user gives Δcoordinates in advance!
    end

    #endregion

    #region making the (symbolic-numerical-hybrid) operator calling the factorial kernels and test functions
    
    l_minus_n = pointsUsedForFields # half the orderU 

    L_MINUS_N = CartesianIndices(Tuple(l_minus_n))
    L_MINUS_N = L_MINUS_N .-L_MINUS_N[1]

    AjiννᶜU = Array{Num,3}(undef,length(multiPointsIndices),NtypeofFields,NtypeofExpr)
    # yes indeed, (νᶜ,) ν, i, j are the oder here

    Ulocal = Array{Num,2}(undef,length(multiPointsIndices),NtypeofFields)
    for iField in eachindex(fields)
        newstring=split(string(fields[iField]),"(")[1]
        Ulocal[:,iField]=Symbolics.variables(Symbol(newstring),1:length(multiPointsIndices))
    end


    AjiννᶜU .= 0
    # the small dictionary map should be here (not inside the loop) but I am too tired that I let this go
    # why tired? since I need to prepare another set of theDiffNu that can run from minus to plus 

    middleLinearν = nothing # only valid for testOnlyCentre

    for iExpr in eachindex(exprs) # j in eq. 42
        for iField in eachindex(fields) # i in eq. 42
            α = bigα[iExpr,iField]
            
            for ν in multiPointsIndices # the relative centre of the local coordinates

                linearν = LinearIndices(multiPointsIndices)[ν]
                CoefU = 0
                
                if !testOnlyCentre || Tuple(ν) === ((Tuple(multiPointsIndices[end] )).-1 ).÷2 .+1 

                    middleLinearν = linearν # only valid for testOnlyCentre
                    tmpCˡη=nothing

                    if testOnlyCentre # Cⁿη size is not the same
                        tmpCˡη=Cˡη
                    else
                        tmpCˡη=Cˡη[:,:,linearν]
                    end

                    for νᶜ in multiPointsIndices 

                        linearνᶜ = LinearIndices(multiPointsIndices)[νᶜ]
                        #relativeDistanceνᶜ = Δ .* car2vec(νᶜ-ν)
                        #relativeDistanceνᶜ = car2vec(νᶜ-ν)
                        #localmapνᶜ = Dict(zip(coordinates, relativeDistanceνᶜ))
                        
                        #U_HERE = substitute(fields[iField],localmapνᶜ)
                        U_HERE = Ulocal[linearνᶜ,iField]
                        
                        for eachα in α
                            nodeValue=eachα.node
                            nᶜ = eachα.nᶜ
                            n = eachα.n

                            for ηᶜ in multiPointsIndices

                                linearηᶜ = LinearIndices(multiPointsIndices)[ηᶜ]
                                #relativeDistanceηᶜ = Δ .* car2vec(ηᶜ-ν)
                                #relativeDistanceηᶜ = car2vec(ηᶜ-ν)
                                #localmapηᶜ = Dict(zip(coordinates, relativeDistanceηᶜ))
                                localmapηᶜ=Dict()
                                for iVar in eachindex(vars)
                                    localmapηᶜ[vars[iVar]]=varM[iVar,linearηᶜ][]
                                end
                              
                                for l in n .+ L_MINUS_N
                                    linearl = LinearIndices(multiLCar)[l]
                                    for lᶜ in nᶜ.+L_MINUS_N
                                        linearlᶜ = LinearIndices(multiLCar)[lᶜ]
                                        kernelProducts = 1
                                        for iCoord in eachindex(coordinates)
                                            l_n_field = Tuple(l-n)[iCoord]
                                            l_n_variable = Tuple(lᶜ-nᶜ)[iCoord]
                                            # here I take only the middle_value
                                            kernelProducts*=integralBsplineTaylorKernels1D(1,Δ[iCoord],l_n_variable,l_n_field)[1]
                                        end
                                        
                                        #nodeValue=Symbol(nodeValue)
                                        #@show localExpression=substitute(nodeValue,localmap)
                                        #@show typeof(nodeValue)
                                        #newExpr = mySimplify.(map((e) -> substitute(e, Dict(localmap)), nodeValue))
                                        
                                        substitutedValue = substitute(nodeValue, localmapηᶜ)

                                        CoefU +=tmpCˡη[linearηᶜ,linearlᶜ]*tmpCˡη[linearνᶜ,linearl]*kernelProducts*substitutedValue*U_HERE
                                    end
                                end
                            end

                            

                        end
                    end
                end

                AjiννᶜU[linearν,iField,iExpr] += CoefU
            end
        end
    end

    #endregion

    #region outputs
    
    utilities=(middlepoint=middleLinearν,localPointsIndices=multiPointsIndices,localMaterials=varM,localFields=Ulocal)
    if testOnlyCentre
        smallAjiννᶜU = AjiννᶜU[middleLinearν,:,:]
        return smallAjiννᶜU,utilities
    else
        return AjiννᶜU,utilities
    end

    #endregion
    
end

function constructingEquations(AjiννᶜU,Γg,coordinates,models,exprs,fields,vars,modelPoints,utilities;absorbingBoundaries=nothing,initialCondition=0.0)

    #todo list
    #
    # I have to include some complex initial condition for 場
    #
    # external forces: do we need to construct a big matrix? 

    
    #region general introduction
    #
    #
    # after the construction of local (semi-)symbolic expressions with linearised operators#
    # here we will read the model parameters and construct the numerical operators
    #
    # Nobuaki Fuji @ IPGP/UPC/IUF since 2024
    #
    # 
    # encouraged by Thibault Duretz @ U. Frankfurt Goethe, Kurama Okubo @ NIED
    #
    # Julia hackathon October 2024, March 2025
    #
    #
    # intermediate presentations: IPGP-CIA workshop October 2024; IPGP-ERI workshop November 2024; lighttalk @ systemI December 2024
    #             EGU @ Vienna May 2025
    #    
    #     Fuji & Duretz in preparation
    #
    #
    #endregion

    #region unpacking, N-dimensionalising all the models 

    timeMarching = any(a -> a === t, coordinates) # important to know if we need to construct a time marching scheme

    @unpack middlepoint,localPointsIndices,localMaterials,localFields = utilities
    Ndimension=length(coordinates)
   
    modelPoints=collect(modelPoints)

    Models=[]

    NtypeofMaterialVariables = length(vars)
    NtypeofFields = length(fields)
    NtypeofExpr = length(exprs)

    Models=Array{Any,1}(undef,NtypeofMaterialVariables)

    if length(models) !== NtypeofMaterialVariables 
        @error "Each material has to have its own model"
    end
    
    
    for iVar in eachindex(vars)
        CartesianDependency=findCartesianDependency(vars[iVar],length(coordinates))
        if ndims(models[iVar]) === CartesianDependency
            @error "model parameter dimension is not what you declared in the equation!"
        end
        if sum(CartesianDependency) === 0 # when it is a constant
            tmpModel=Array{Any,Ndimension}(undef,(ones(Int, Ndimension)...)...)
            tmpModel[CartesianIndex(Tuple(ones(Int, Ndimension))...)] = models[iVar]
            Models[iVar]=tmpModel
        else
            newCoords=expandVectors(size(models[iVar]),CartesianDependency)
            @show newCoords
            
            tmpModel=reshape(models[iVar],newCoords...)
            Models[iVar]=tmpModel

            for iCoord in eachindex(newCoords)
                if newCoords[iCoord]!== modelPoints[iCoord] || newCoords[iCoord] !== 1
                    @error "the model should have the same dimension! (or constant)"
                end
            end
        end
     
    end

    #endregion
    
    #region construction of the fields

    wholeRegionPoints = nothing

    if absorbingBoundaries === nothing
        wholeRegionPoints=modelPoints
    else
        # absorbingBoundaries should be two column array 
        if size(absorbingBoundaries)[1] !== 2
            @error "you have to give us the left and right values for absorbing boundaries"
        elseif size(absorbingBoundaries)[2] !== size(modelPoints)[1] && !timeMarching
            @error "you have to give us the values for each direction for absorbing boundaries"
        elseif size(absorbingBoundaries)[2] === size(modelPoints)[1]-1 && timeMarching
            absorbingBoundaries=[absorbingBoundaries; 0 0]
        end
        wholeRegionPoints=modelPoints.+sum(absorbingBoundaries,1)
    end

    # below should be parallelised at some points

    @show size(localPointsIndices)

    timeSteps=1

    if timeMarching
        @show pointsForOneStep=car2vec(localPointsIndices)[end]
        timeSteps=wholeRegionPoints[end]-pointsForOneStep+1
        @show wholeRegionPoints[end]=pointsForOneStep
    end

    場=Array{Any,1}(undef,NtypeofFields)
    for iField in eachindex(fields)
        newstring=split(string(fields[iField]),"(")[1]*"_mod"
        場[iField]=string_as_varname(newstring, Array{Any,Ndimension}(undef,Tuple(wholeRegionPoints)))
    end
    場 .= initialCondition
    
    #endregion

    #region we construct the numerical operators for each test functions that is related to the points




    #endregion

end