using Symbolics,UnPack,LinearAlgebra,DrWatson

include("../src/batchNewSymbolics.jl")
include("../src/batchUseful.jl")
include("../src/batchDrWatson.jl")
include("../src/CerjanBoundary.jl")
include("../src/IntegrateBsplineAndPolynomials.jl")

# PDECoefFinder cannot detect the material partials × material partials for the moment!! 

timeDimensionString="t" 
# if the user does not want to use "t" for the time marching scheme, it should be changed
# and this "t" should be the last element in coordinates

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

function varMmaker(maxPointsUsed,coordinates,vars) 
    # this will make an array of material coeffs for with a local Cartesian grid (max points used for a node)
    Ndimension = length(coordinates)

    R = CartesianIndices(Tuple(maxPointsUsed))

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
    
    return varM
end

function PDECoefFinder(orders,coordinates,expr,field,vars)
    # PDECoefFinder cannot detect the material partials × material partials for the moment!! 
    # I know how to do it, but eq. 40 should be then more generalised (kind of the product of partials of different materials)

    # maxPolynomialOrderMaterial is also a chelou thing, that I need to work on more systematically
    # like the powers of partials should also be included but here search for Rm[1], yeah, that's what I am doing


    Ndimension = length(coordinates)
    alpha=[]
    
    maxPolynomialOrderMaterial = 2*(maximum(orders)-1)
    ∇=makeMixPartials(orders,coordinates;field=field)
    R=CartesianIndices(∇)
    expr=mySimplify(expr)


    for i in R
        term_searched = ∇[i]

        tmpCoeff = myCoeff(expr,term_searched)
        if tmpCoeff !== 0
            isTmpCoeffAConstant=true
            for iVar in eachindex(vars)
                
                ∇m=makeMixPartials(orders,coordinates;field=vars[iVar]) # material partials
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
                                ∇n = makeMixPartials(orders,coordinates;field=vars[jVar]) 
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


    return alpha # varM: iVar and linearised cartesian indices
end 

function TaylorCoefInversion(coefInversionDict::Dict)

    # the user might want to have a look on illposedTaylorCoefficientsInversion_legend, which is deprecated as of 10/06/2025.


    # based on the equation 27 (of the version 10/06/2025 FD2025 : \psi_{;\mu,\nu}^{(l)}[\nu+\mu]=\sum_\eta C_{\mu+\eta;\mu,\nu}^{(l)} \psi[\nu+\mu+\eta]), we need to perform this inversion anyways for all the point \mu inside L(\nu) (the concerned points for \nu)


    # be careful that pointsIndices is now a 1D array of integer vectors!!

    @unpack coordinates,multiOrdersIndices,pointsIndices, Δ = coefInversionDict

    Ndimension=length(coordinates)

    if length(Δ) !== Ndimension
        @error "the numerical delta increment has not the same dimension!"
    end


    numberOfEtas = length(pointsIndices)
    numberOfLs   = length(multiOrdersIndices)

    CˡηGlobal = Array{Any,3}(undef,numberOfEtas,numberOfLs,numberOfEtas)

    # this is the C^{(l)}_{\mu+\eta; \mu, \nu}
   
    for k in eachindex(pointsIndices)
        CˡηGlobal[:,:,k]=TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,Δ,k)
    end 

    return CˡηGlobal

end

function TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,Δ,k)
    # the old version is : illposedTaylorCoefficientsInversionSingleCentre
    TaylorExpansionCoeffs = Array{Num,2}(undef,numberOfLs,numberOfEtas)
    for i in eachindex(pointsIndices)
        η = pointsIndices[i]-pointsIndices[k]
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
    Cˡηlocal=invaa*transpose(TaylorExpansionCoeffs)
    return Cˡηlocal
end



function illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=true,Δ=nothing,timeMarching=false)

    # this function is deprecated, no more testOnlyCentre/timeMarching options are allowed (nor Symbolics Δ)

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
   

    tmpVecForMiddlePoint = (car2vec(multiPointsIndices[end]).-1 ).÷2 .+1 # only valid for testOnlyCentre
    midTimeCoord=nothing
    if timeMarching
        midTimeCoord=car2vec(multiPointsIndices[end])[end]-1
        tmpVecForMiddlePoint[end]=midTimeCoord
    end
    midK=vec2car(tmpVecForMiddlePoint)

    for k in multiPointsIndices
        linearK = LinearIndices(multiPointsIndices)[k]
        
        if !testOnlyCentre || k === midK || (timeMarching && car2vec(k)[end] === midTimeCoord && !testOnlyCentre) # because we cannot predict more than one futures
            CˡηGlobal[:,:,linearK]=illposedTaylorCoefficientsInversionSingleCentre(numberOfLs,numberOfEtas,multiOrdersIndices,multiPointsIndices,Δ,k)
        end
    end 

    if testOnlyCentre
        CˡηCentre = CˡηGlobal[:,:,LinearIndices(multiPointsIndices)[midK]]
        CˡηGlobal = nothing
        return CˡηCentre,Δ,multiOrdersIndices
    else
        #@show CˡηGlobal
        return CˡηGlobal,Δ,multiOrdersIndices
    end
end

function illposedTaylorCoefficientsInversionSingleCentre(numberOfLs,numberOfEtas,multiOrdersIndices,multiPointsIndices,Δ,k)

    # this function is deprecated as of 10/06/2025

    # this should be completely numerical
    TaylorExpansionCoeffs = Array{Num,2}(undef,numberOfLs,numberOfEtas)
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
    Cˡηlocal=invaa*transpose(TaylorExpansionCoeffs)
    return Cˡηlocal
end

function integralBsplineTaylorKernels1DWithWindow1D(BsplineOrder,WBsplineOrder,μᶜ,μ,ν,L,Δ,l_n_variable,l_n_field)
    # this computes the analytical value of the 1D integral between B-spline fns and weighted Taylor kernels
    # \int dx Bspline Y_μᶜ Y_μ  K_{lᶜ-nᶜ}(y-y_μᶜ) K_{l-n}(y-y_μ)

    # unlike the previous integralBsplineTaylorKernels1D, it computes for a specific ν
    # Cˡη;μ are computed for a specific geometry, so even though the boundary condition reduce
    # the number of available points, each Taylor expansion for K_{l-n}(y-y_μ) should be Okubo
    
    # however, the 'forgotten' μ (due to the whole) should be treated carefully 
    # (which I do not yet implemented here): I think Y_ignoredμ should be added to the Y_availableneighbouringμ

    # or maybe the 'forgotten' μ is anyways not available (and thus very probably not continuous)
    # so we just let this be forgotten 

    kernelValue=0.0
   

    if BsplineOrder=== -1
        # this is for a delta function
        if l_n_variable === 0 && l_n_field === 0
            kernelValue=1.0
        else
            kernelValue=0
        end
    else
        maximumOrder = (BsplineOrder,WBsplineOrder)
        
        params=@strdict maximumOrder numberNodes = L

        output,_=@produce_or_load(BsplineTimesPolynomialsIntegrated,params,datadir("BsplineInt");filename = config -> savename("Bspline",params))
        nodeIndices,nodesSymbolic,b_deriv,integral_b,Δx,extFns,x =output["BsplineIntegraters"]

        # here we make a function Y_μ' Y_μ K_μ' K_μ (details ommitted)
        # note that ν is somewhere middle or at extremeties and 'ν+' expression is ommitted 

        Y_μᶜ=[:,μᶜ,1,WBsplineOrder]
        Y_μ =[:,μ ,1,WBsplineOrder]
        K_μᶜ=(x-nodesSymbolic[μᶜ])^l_n_variable/BigInt(factorial(l_n_variable))
        K_μ =(x-nodesSymbolic[μ])^l_n_field/BigInt(factorial(l_n_field))

        

        dictionaryForSubstitute = Dict()
    
    
        for i in 1:1:maximumOrder-1
            taylorNum *= N+i
            dictionaryForSubstitute[gvec[i]]=x^(N+i)/taylorNum
        
        end


    end


end



function integralBsplineTaylorKernels1D(BsplineOrder,Δ,l_n_variable,l_n_field)

    #call_new_fn(fn, args...) = Base.invokelatest(fn, args...)
    
    # this will compute \int dx Bspline K_{l-n} K_{lᶜ-nᶜ}
    middle_value = 0
    extreme_value = 0
    midPoint = BsplineOrder+2
    maxPoint = (BsplineOrder+1)*2 + 3
    nearboundaries_values=Array{Any,1}(undef,maxPoint)

    if BsplineOrder=== -1
        # this is for a delta function
        if l_n_variable === 0 && l_n_field === 0
            middle_value=1
            nearboundaries_values=ones(3)
        else
            middle_value=0
            nearboundaries_values=zeros(3)
        end

    elseif BsplineOrder >= 0
        maximumOrder = BsplineOrder
        params=@strdict maximumOrder
        output,_=@produce_or_load(BsplineTimesPolynomialsIntegrated,params,datadir("BsplineInt");filename = config -> savename("Bspline",params))
        numberNodes,integral_b_polys,N,Δx=output["BsplineIntegraters"]
        #fns=eval.(build_function.(integral_b_polys,N,Δx))
        middleNode = numberNodes ÷ 2
        #@show call_new_fn(fns[middleNode, BsplineOrder+1], l_n_variable + l_n_field + 1, Δ)
        #@show middle_value = call_new_fn(fns[middleNode,BsplineOrder+1](l_n_variable+l_n_field+1,Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #@show integral_b_polys[middleNode, BsplineOrder+1]
     
        middle_value = Symbolics.substitute(integral_b_polys[middleNode, BsplineOrder+1],Dict(N=>l_n_variable+l_n_field,Δx => Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #@show (Δ^(l_n_variable+l_n_field+1)-(-Δ)^(l_n_variable+l_n_field+1))/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        for iNode in 1:midPoint+2
            nearboundaries_values[iNode] = Symbolics.substitute(integral_b_polys[iNode, BsplineOrder+1],Dict(N=>l_n_variable+l_n_field,Δx => Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
            nearboundaries_values[maxPoint-iNode+1] = Symbolics.substitute(integral_b_polys[numberNodes-iNode+1, BsplineOrder+1],Dict(N=>l_n_variable+l_n_field,Δx => Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        end
    end

#=
    elseif BsplineOrder === 0
        #
    elseif BsplineOrder=== 1
        middle_value = (Δ^(l_n_variable+l_n_field+1)-(-Δ)^(l_n_variable+l_n_field+1))/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #extreme_value = (Δ^{l_n_variable+l_n_field+1})/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(l_n_variable)*factorial(l_n_field))
    end
    =#


    return middle_value,nearboundaries_values
end

function spaceCoordinatesConversionfunctions(absorbingBoundaries, NdimensionMinusTime)
    offset_model = vec2car(absorbingBoundaries[1, 1:NdimensionMinusTime])
    #offset_empty = vec2car(spacePointsUsed)

    model2whole(a::CartesianIndex) = a + offset_model
    whole2model(a::CartesianIndex) = a - offset_model
    #whole2empty(a::CartesianIndex) = a + offset_empty
    #empty2whole(a::CartesianIndex) = a - offset_empty
    #model2empty(a::CartesianIndex) = whole2empty(model2whole(a))
    #empty2model(a::CartesianIndex) = whole2model(empty2whole(a))
    return(; model2whole, whole2model)
    #return (; model2whole, whole2model, whole2empty, empty2whole, model2empty, empty2model)
end

function BouncingCoordinates(a::CartesianIndex,PointsUsed)
    #
    # this will bounce the boundary inside the PointsUsed vector
    #
    # i.e. get the nearby coordinates inside the domain to fake the continuity

    if length(a) !== length(PointsUsed)
        @error "cannot bound this CartesianIndex due to the dimension mismatch"
    end
    avector=car2vec(a)
    for iCoord in eachindex(avector)
        if avector[iCoord] < 1
            avector[iCoord] = 1
        elseif avector[iCoord] > PointsUsed[iCoord]
            avector[iCoord] = PointsUsed[iCoord]
        end
    end
    a=vec2car(avector)
    return a
end

function OPTobj(operatorConfigurations::Dict)
    # this is just a wrapper for the OPTobj function below, for DrWatson package
    @unpack famousEquationType, Δnum, orderBtime, orderBspace, pointsInSpace, pointsInTime,IneedExternalSources, iExperiment= operatorConfigurations
    exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = famousEquations(famousEquationType)
  
    trialFunctionsCharacteristics=(orderBtime=orderBtime,orderBspace=orderBspace,pointsInSpace=pointsInSpace,pointsInTime=pointsInTime)
    @time operatorData=OPTobj(exprs,fields,vars; coordinates=coordinates,trialFunctionsCharacteristics=trialFunctionsCharacteristics,Δnum = Δnum,iExperiment=iExperiment)
    #AjiννᶜU=operatorData[1]
    #utilities=operatorData[2]   

    operatorForceData=nothing
    # if you do not want to apply external forces, it is possible to skip below
    if IneedExternalSources 
        @time operatorForceData=OPTobj(extexprs,extfields,extvars; coordinates=coordinates,trialFunctionsCharacteristics=trialFunctionsCharacteristics,Δnum = Δnum,iExperiment=iExperiment)  
        #@show Γg = operatorForceData[1]
        #utilitiesForce = operatorForceData[2]
    end
    eqInfo=(exprs=exprs,fields=fields,vars=vars,extexprs=extexprs,extfields=extfields,extvars=extvars,coordinates=coordinates)
    operators=(operatorPDE=operatorData, operatorForce=operatorForceData,eqInfo=eqInfo)
    return @strdict(operators)
end

function OPTobj(exprs,fields,vars; coordinates=(x,y,z,t), TaylorOptions=(WorderBtime=1,WorderBspace=1,supplementaryOrder=2), trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, pointsInSpace=2,pointsInTime=2),CˡηSymbolicInversion=false,testOnlyCentre=true,Δnum = nothing,iExperiment =nothing)

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



    timeMarching = any(a -> a === timeDimensionString, string.(coordinates))


    @unpack orderBtime, orderBspace, pointsInSpace, pointsInTime = trialFunctionsCharacteristics
    @unpack WorderBtime, WorderBspace,supplementaryOrder = TaylorOptions

    NtypeofExpr=length(exprs)   # number of governing equations
    NtypeofMaterialVariables=length(vars) # number of material coefficients
    NtypeofFields=length(fields) # number of unknown fields
    
    Ndimension = length(coordinates) # we do not change this for the moment, especially for the time-marching scheme
    pointsUsed = ones(Int, Ndimension).*(pointsInSpace+1)
    if timeMarching
        pointsUsed[end]=pointsInTime+1
    end


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
    WorderBspline=zeros(Int,Ndimension)

    if timeMarching
        orderBspline[Ndimension]=orderBtime*fieldDependency[Ndimension]
        orderBspline[1:Ndimension-1]=orderBspace*fieldDependency[1:Ndimension-1]
        WorderBspline[Ndimension]=WorderBtime*fieldDependency[Ndimension]
        WorderBspline[1:Ndimension-1]=WorderBspace*fieldDependency[1:Ndimension-1]
    else
        orderBspline[1:Ndimension]=orderBspace*fieldDependency[1:Ndimension]
        WorderBspline[1:Ndimension]=WorderBspace*fieldDependency[1:Ndimension]
    end
    
    # the maximum number of points used in the vicinity of the node, which is independent of the order of B-spline functions (see our paper)
    pointsUsedForFields=(pointsUsed.-1).*fieldDependency.+1

    # orderExpressions is the maximal orders of partials that we could expect in the expressions
    orderExpressions=pointsUsedForFields
    
    # numbers of points to evaluate the integral for the governing equation filtered by the test functions
    
    # orderU is the maximum orders for the fields that we will use for OPT coefficients' exploration
    @show orderU = (orderExpressions .-1) .+ (supplementaryOrder .*fieldDependency).+1 
    # we restore this orderU since we need to control this 

    #endregion

    #region analysis of expressions to obtain the α_{n'nji}

    bigα=Array{Any,2}(missing,NtypeofFields,NtypeofExpr)
    varM=nothing
    for iExpr in eachindex(exprs)
        for iField in eachindex(fields)
            
            tmpNonZeroAlphas=PDECoefFinder(orderExpressions,coordinates,exprs[iExpr],fields[iField],vars) 
            # we assume that the pointsUsedForFields represent the highest order of partials
            bigα[iField,iExpr]=unique(tmpNonZeroAlphas)
        end
    end
    varM=varMmaker(pointsUsedForFields,coordinates,vars)
    @show bigα,varM
    #endregion

    #region Preparation for Taylor expansion
    
    orderTaylors=Array{Any,Ndimension}(undef,Tuple(orderU))
    pointsInSpaceTime=Array{Any,Ndimension}(undef,Tuple(pointsUsedForFields))
    
    
    multiOrdersIndices=CartesianIndices(orderTaylors)

    availablePointsConfigurations = []
    centrePointConfigurations=[]

    #endregion

    #region Cartesian indices that can be available to use (normally: iGeometry=1)

    multiPointsIndices=CartesianIndices(pointsInSpaceTime)
    # this is the whole local Cartesian grids (without any lacking points)

    tmpVecForMiddlePoint = (car2vec(multiPointsIndices[end]).-1 ).÷2 .+1 # only valid for testOnlyCentre
    midTimeCoord = nothing
    if timeMarching
        midTimeCoord=car2vec(multiPointsIndices[end])[end]-1
        tmpVecForMiddlePoint[end]=midTimeCoord
        #AjiννᶜU = Array{Num,2}(undef,length(multiPointsIndices)÷(midTimeCoord+1),NtypeofExpr)
    end
    middleν=vec2car(tmpVecForMiddlePoint)


    availablePointsConfigurations=push!(availablePointsConfigurations,car2vec.(multiPointsIndices))
    centrePointConfigurations=push!(centrePointConfigurations,LinearIndices(multiPointsIndices)[middleν])

    #endregion


    #region obtaining the semi-symbolic expression of cost function based on eqns. 52 and 53.

    # before calling AuSymbolic we can manipulate pointsIndices for various boundary configurations


    if CˡηSymbolicInversion # this seems super cool but it takes time
        #Cˡη,Δ,multiLCar = illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=testOnlyCentre,timeMarching=timeMarching)
        Δ = Symbolics.variables(:Δ,1:Ndimension)
    else
        Δ = Δnum
    end

    AjiννᶜU=[]
    Ulocal=[]

    for iConfigGeometry in eachindex(availablePointsConfigurations) 
        pointsIndices=availablePointsConfigurations[iConfigGeometry]
        middleLinearν=centrePointConfigurations[iConfigGeometry]
        #varM is given above for the max number of points used 
        tmpAjiννᶜU,tmpUlocal=AuSymbolic(coordinates,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δ,varM,bigα,orderBspline,WorderBspline)
        AjiννᶜU=push!(AjiννᶜU,tmpAjiννᶜU)
        Ulocal=push!(Ulocal,tmpUlocal)
    end


    #endregion



    #region outputs
    
    utilities=(middlepoint=middleν,middlepointLinear=middleLinearν,localPointsIndices=multiPointsIndices,localMaterials=varM,localFields=Ulocal)
    if testOnlyCentre
        smallAjiννᶜU = Array{Num,2}(undef,1,NtypeofExpr) # shrinking but the dimension is still the same
        smallAjiννᶜU[1,:] = AjiννᶜU[middleLinearν,:]
        return smallAjiννᶜU,utilities
    else
        return AjiννᶜU,utilities
    end

    #endregion
    
end


function AuSymbolic(coordinates,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δ,varM,bigα,orderBspline,WorderBspline)

    # the contents of OPTobj which is now renamed as AuSymbolic since we compute Au for different pointsIndices

    #region obtaining Cˡη either symbolically either with Δcoordinates in a numerical way


    coefInversionDict = @strdict coordinates multiOrdersIndices pointsIndices Δ
    Cˡη, _ = produce_or_load(TaylorCoefInversion,coefInversionDict,datadir("taylorCoefInv");filename = config -> savename("TaylorInv",coefInversionDict))
   

    #endregion

    #region making the (symbolic-numerical-hybrid) operator calling the factorial kernels and test functions

    L_MINUS_N = multiOrdersIndices
    L_MINUS_N = L_MINUS_N .-L_MINUS_N[1]

    
    # the order is: (νᶜ,) ν, i, j  here

    Ulocal = Array{Num,2}(undef,length(pointsIndices),NtypeofFields)
    for iField in eachindex(fields)
        newstring=split(string(fields[iField]),"(")[1]
        Ulocal[:,iField]=Symbolics.variables(Symbol(newstring),1:length(pointsIndices))
    end

    AjiννᶜU = Array{Num,1}(undef,NtypeofExpr)
    
    # this is the cost function for ν point so the number of elements is just the number of expressions (governing equations)
    

    AjiννᶜU .= 0

    for iExpr in eachindex(exprs) # j in eq. 52
        for iField in eachindex(fields) # i in eq. 52
            α = bigα[iExpr,iField]

            # this is the ν point in the relative space-time domain

            linearν = middleLinearν  

            CoefU = 0

            for linearμᶜ in eachindex(pointsIndices)
                    
                tmpCˡημᶜ=Cˡη[:,:,linearμᶜ] # C^{(l')}_{μ'+η';μ',ν}

                for linearμ in eachindex(pointsIndices)

                    tmpCˡημ=Cˡη[:,:,linearμ] # C^{(l)}_{μ+η;μ,ν}

                    for linearμ_plus_η in eachindex(pointsIndices) # relative position νᶜ-ν

                        U_HERE = Ulocal[linearμ_plus_η,iField]
                        
                        for eachα in α
                            
                            nodeValue=eachα.node
                            nᶜ = eachα.nᶜ
                            n = eachα.n

                            for linearμᶜ_plus_ηᶜ in eachindex(pointsIndices)
                                
                                linearμᶜ_plus_ηᶜ_in_the_whole = LinearIndices(multiOrdersIndices)[vec2car(pointsIndices[linearμᶜ_plus_ηᶜ])]

                                localmapηᶜ=Dict()

                                for iVar in eachindex(vars)
                                    localmapηᶜ[vars[iVar]]=varM[iVar,linearμᶜ_plus_ηᶜ_in_the_whole][]
                                end
                                
                                for l in n .+ L_MINUS_N
                                    if l ∈ L_MINUS_N
                                        linearl = LinearIndices(multiOrdersIndices)[l]
                                        for lᶜ in nᶜ.+ L_MINUS_N
                                            if l ∈ L_MINUS_N
                                                linearlᶜ = LinearIndices(multiOrdersIndices)[lᶜ]
                                                kernelProducts = 1
                                                for iCoord in eachindex(coordinates)
                                                    l_n_field = Tuple(l-n)[iCoord]
                                                    l_n_variable = Tuple(lᶜ-nᶜ)[iCoord]
                                                    # here I take only the middle_value
                                                    #kernelProducts*=integralBsplineTaylorKernels1D(orderBspline[iCoord],Δ[iCoord],l_n_variable,l_n_field)[1]

                                                    kernelProducts*=integralBsplineTaylorKernelsWithWindow1D(orderBspline[iCoord],WorderBspline[iCoord],pointsIndices[linearμᶜ][iCoord],pointsIndices[linearμ][iCoord],pointsIndices[linearν][iCoord],multiOrdersIndices[end][iCoord], Δ[iCoord],l_n_variable,l_n_field)[1]
                                                end
                                                
                                                #nodeValue=Symbol(nodeValue)
                                                #@show localExpression=substitute(nodeValue,localmap)
                                                #@show typeof(nodeValue)
                                                #newExpr = mySimplify.(map((e) -> substitute(e, Dict(localmap)), nodeValue))
                                                
                                                substitutedValue = substitute(nodeValue, localmapηᶜ)

                                                CoefU +=tmpCˡημᶜ[linearμᶜ_plus_ηᶜ,linearlᶜ]*tmpCˡημ[linearμ_plus_η,linearl]*kernelProducts*substitutedValue*U_HERE
                                            end
                                        end
                                    end
                                end
                                
                            end
                            
                        end
              
                    end
                end
            end
       

            AjiννᶜU[iExpr] += CoefU
            
        end
    end



    #endregion
    return AjiννᶜU,Ulocal


end 



function constructingNumericalDiscretisedEquations(config::Dict)
    # just a wrapper
    @unpack semiSymbolicOpt,coordinates,modelName,models,fields,vars,famousEquationType,modelPoints,utilities, maskedRegion, NpointsUsed = config

    costfunctions,場,champsLimité=constructingNumericalDiscretisedEquations(semiSymbolicOpt,coordinates,models,fields,vars,modelPoints,utilities, maskedRegion;initialCondition=0.0)
    numOperators=(costfunctions=costfunctions,場=場,champsLimité=champsLimité)

    #@show costfunctions
    return @strdict(numOperators)
end

function constructingNumericalDiscretisedEquations(semiSymbolicsOperators,coordinates,models,fields,vars,modelPoints,utilities,maskedRegionInSpace;absorbingBoundaries=nothing,initialCondition=0.0)

    #todo list
    # 
    # this function is tooooooo complicated! I think I can simplify very much this!
    #
    #
    #  need to work on the bc, same like the masked thing (limited region of source)
    #
    # absorbing boundaries : I think we can already put the bc inside the numerical operators but be careful with the time marching: search for weightingCerjan
    # 
    # need extend to 4 points with the same test functions (3 points) -> staggered grid
    #  
    # I have to include some complex initial condition for 場
    #
    # have to write:
    #  function illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices,midPoint,Δ)

    
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

    # coordinates: Model: the real model domain; Whole: computation domain with absorbing boundaries; 
    #              Empty: Whole + some more points to avoid missing reference to the field and material (they should be just zeros)
    #

    # intermediate presentations: IPGP-CIA workshop October 2024; IPGP-ERI workshop November 2024; lighttalk @ systemI December 2024
    #             EGU @ Vienna May 2025
    #    
    #     Fuji & Duretz in preparation
    #
    #
    #
    #
    #endregion

    #region unpacking, N-dimensionalising all the models 

    testOnlyCentre = false
 
    if size(semiSymbolicsOperators)[1] === 1
        testOnlyCentre = true
    elseif ndims(semiSymbolicsOperators) !==2
        @error "the semi symbolic operators are not computed correctly!"
    end

    timeMarching = any(a -> a === timeDimensionString, string.(coordinates)) # important to know if we need to construct a time marching scheme

    @unpack middlepoint,middlepointLinear,localPointsIndices,localMaterials,localFields = utilities

    # the last coordinate should be cosidered as time

    if !timeMarching
        localPointsIndices=CartesianIndices(Tuple([car2vec(localPointsIndices[end]);1]))
        middlepoint=CartesianIndex([car2vec(middlepoint);1]...)
        tmpT=Symbolics.variable(timeDimensionString)
        coordinates = (coordinates...,tmpT)
        modelPoints = (modelPoints...,1)
        tmp_del = Symbolics.variable("∂"*timeDimensionString)
        tmp_del = Differential(tmpT)
        ∂ .= push!(∂,tmp_del)
        tmp_del_2 = Symbolics.variable("∂"*timeDimensionString*"²")
        tmp_del_2 = Differential(tmpT)*Differential(tmpT)
        ∂² .= push!(∂²,tmp_del_2)
    end

    #@show coordinates,∂,∂²

    localPointsSpaceIndices=CartesianIndices(Tuple(car2vec(localPointsIndices[end])[1:end-1]))
    Ndimension=length(coordinates)
    
   
    modelPoints=collect(modelPoints)

    Models=[]

    NtypeofMaterialVariables = length(vars)
    NtypeofFields = length(fields)
    NtypeofExpr = size(semiSymbolicsOperators)[end]

    Models=Array{Any,1}(undef,NtypeofMaterialVariables)
    ModelPoints=Array{Int,2}(undef,Ndimension,NtypeofMaterialVariables)

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
            ModelPoints[:,iVar] = ones(Int, Ndimension)
            tmpModel[vec2car(ones(Int, Ndimension))] = models[iVar]
            Models[iVar]=tmpModel
        else
            #@show models[iVar],iVar,CartesianDependency, vars[iVar]
            newCoords=expandVectors(size(models[iVar]),CartesianDependency)
            ModelPoints[:,iVar] = newCoords
 
            tmpModel=reshape(models[iVar],newCoords...)
            Models[iVar]=tmpModel

            for iCoord in eachindex(newCoords)
                if newCoords[iCoord]!== modelPoints[iCoord] && newCoords[iCoord] !== 1
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
        absorbingBoundaries = zeros(Int,2, Ndimension)
    elseif absorbingBoundaries === "CerjanBoundary"
        wholeRegionPoints=modelPoints
        absorbingBoundaries = ones(Int,2, Ndimension-1)*CerjanGridPoints
        absorbingBoundaries=[absorbingBoundaries; 0 0]
        wholeRegionPoints=modelPoints.+sum(absorbingBoundaries,1) 
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


    # some useful stuff


    spacePointsUsed=car2vec(localPointsIndices[end])[1:end-1]
    timePointsUsedForOneStep=car2vec(localPointsIndices[end])[end]
    wholeRegionPointsSpace=wholeRegionPoints[1:end-1]
   
 
    # we need to put the left and right regions in order that centre ν configuration can pass

    #emptyRegionPointsSpace=wholeRegionPointsSpace.+ 2 .* spacePointsUsed

    #場dummy=Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
    場 = Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
    

    for it in 1:timePointsUsedForOneStep
        for iField in eachindex(fields)
            newstring=split(string(fields[iField]),"(")[1]*"_mod"*"_t="*string(it)
            場[iField,it]=Symbolics.variables(Symbol(newstring),Base.OneTo.(Tuple(wholeRegionPointsSpace))...)
        end
    end

    #since everything is super clumsy, here we make several useful functions to change one coordinate to another
    
    conv=spaceCoordinatesConversionfunctions(absorbingBoundaries[:,1:end-1], Ndimension-1)

    #endregion 

    #region making a maskingField (for limited source areas, boundary conditions, etc.)

    maskingField=Array{Any,Ndimension-1}(undef,Tuple(wholeRegionPointsSpace)) # maskingField is defined only for whole domain
    champsLimité = nothing
    if maskedRegionInSpace === nothing
        maskingField .= 1.0
    elseif typeof(maskedRegionInSpace) === Array{CartesianIndex,1}
        champsLimité = Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
        for it in 1:timePointsUsedForOneStep
            for iField in eachindex(fields)
                newstring=split(string(fields[iField]),"(")[1]*"_mod_limited"*"_t="*string(it)
                champsLimité[iField,it] = Array{Any,1}(undef,length(maskedRegionInSpace))
            end
        end
        maskingField .= 0.0
        tmpIndex=1
        for iSpace in maskedRegionInSpace
            jSpace = conv.model2whole(iSpace)
            maskingField[jSpace] =1.0
            for it in 1:timePointsUsedForOneStep
                for iField in eachindex(fields)
                    
                    #tmpChampsLimitéContents= (jSpace,場[iField,it][jSpace])
                    champsLimité[iField,it][tmpIndex]=場[iField,it][jSpace]
                end
            end
            tmpIndex += 1
        end
    else
        @error "maskedRegionInSpace should be a 1D array of CartesianIndex (if it is CartesianIndices, you need to collect(Tuple()))"
    end


    #endregion

    #region relative ν to be considered, especially around the boundaries, useful for the following sections

    PointsSpace=CartesianIndices(Tuple(wholeRegionPointsSpace))
    NpointsSpace=length(PointsSpace) # number of points in space

    NtestfunctionsInSpace=NpointsSpace # this assumption is valid only for test functions related to grid points

    νWhole=Array{Any,1}(undef,NtestfunctionsInSpace) # the coordinate in wholeRegionPointsSpace: for the moment mapping from testfunction to ν is bijective

    # below is only for the bijective projection between test functions and ν

    for iPoint in eachindex(νWhole)

        νWhole[iPoint] = PointsSpace[iPoint] # this should be not true for higher B-spline test functions

    end


    νRelative=Array{Any,1}(undef,NtestfunctionsInSpace) # the relative coordinate to take (the coordinate used for the semi-symbolic operator derivation)
    νRelative.=middlepoint


    
    #endregion

    #region we construct the numerical operators for each test function that is related to its corresponding point

    # first we compute the νRelative more seriously if testOnlyCentre we might do nothing at all

    if testOnlyCentre
        # If we compute only the operators without boundaries, we use kind of 'truncated' crazy operators 
        # derived at the centre point and we do not talk about it, just believe the absorbing boundaries
        # like, tant pis, il n'y a pas de points donc j'ignore juste !

        # maybe we do not do anything!

    else

        # we use this clause only if we are interested in a serious boundary conditions, 
        # i.e. not the artificial cartesian box boundaries (we can just apply some stupid absorbing boundaries in that case)
        # Hence this clause should be more generalised even it kills the performance
        #  like, we give an array of free surface or discontinuities in CartesianIndex arrays 
        # and we look all the points concerned

        # we need to explore everywhere in the wholeRegionPoints! Free surface etc. should be very much affected 
        #
      
        boundaryPointsSpace=[]
        for iDimSpace in 1:Ndimension-1 # we take care of the boundaries of the Cartesian box (it should be the same for internal/external topography)
            # points concerned
            leftstart=1
            leftend=spacePointsUsed[iDimSpace]÷2
            rightstart=NpointsSpace[iDimSpace]-spacePointsUsed[iDimSpace]÷2+1
            rightend=NpointsSpace[iDimSpace]
            # suppose that the domain is sufficiently big (maybe it is not the case for some crazy topography ...)
            for iCoord in range(leftstart,leftend)
                


            end

            for iCoord in range(rightstart,rightend)

            end
        end
    end

    # here the number of test functions should not be necessarily the number of points but I will work later

    costFunctions=Array{Any,2}(undef,NtypeofExpr,NtestfunctionsInSpace)

    #@show semiSymbolicsOperators
    #@show localMaterials[1,15],localFields,size(localFields)
    #@show Models[1][10,15,1]

    for iTestFunctions in eachindex(νWhole)
        # here each test function is connected to one ν point 
        # We need to be careful that this can be no more true for different basis functions other than linear B-spline
        iPoint = iTestFunctions 
        νtmpWhole=νWhole[iPoint]
        # be careful with the two lines above, they are based on the assumption that each test function is linked to one collocated point


        #νtmpModel=conv.whole2model(νtmpWhole)
        νᶜtmpWhole = localPointsSpaceIndices .+ (νtmpWhole - carDropDim(νRelative[iPoint])) # this is the shift vector
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)

        # examine νᶜtmpWhole if it is out of the range 



        for iExpr in eachindex(semiSymbolicsOperators[1,:])

            tmpMapping=Dict()


            for iT in 1:timePointsUsedForOneStep
                
                for iVar in eachindex(vars)
                    
                    spaceModelBouncedPoints=ModelPoints[1:end-1,iVar]

                    if ModelPoints[end,iVar] > 1
                        iiT=iT
                    else
                        iiT = 1
                    end


                    # model parameters should be bounced at the whole region limits
                    νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                    for jPoint in νᶜtmpWhole
                        jPointLocal = jPoint - νtmpWhole + carDropDim(νRelative[iPoint])
                        jPointTLocal = carAddDim(jPointLocal,iT)
                        linearjPointTLocal=LinearIndices(localPointsIndices)[jPointTLocal]
              
                        tmpMapping[localMaterials[iVar,linearjPointTLocal]] = Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal],iiT)]
                        
                    end

                end

                for jPoint in νᶜtmpWhole
                    #@show iPoint, jPoint
                    jPointLocal = jPoint - νtmpWhole + carDropDim(νRelative[iPoint])
                    jPointTLocal = carAddDim(jPointLocal,iT)
                    linearjPointTLocal=LinearIndices(localPointsIndices)[jPointTLocal]
                    #jPointT=carAddDim(jPoint,iT)
                    #linearjPointT=LinearIndices(localPointsIndices)[jPointT]
                    for iField in eachindex(fields)
                        if is_all_less_than_or_equal(CartesianIndex(ones(Int, Ndimension-1)...),conv.whole2model(jPoint)) && is_all_less_than_or_equal(conv.whole2model(jPoint),vec2car(ModelPoints[1:end-1]))
                            # when it is inside the model domain box
                            tmpMapping[localFields[linearjPointTLocal,iField]] = 場[iField,iT][jPoint]*maskingField[jPoint]

                        elseif is_all_less_than_or_equal(νWhole[1],jPoint) && is_all_less_than_or_equal(jPoint,νWhole[end])
                            # if it is in the absorbing boundary zones we apply a simple Cerjan
                            if iT === timePointsUsedForOneStep # the last one (the future) will be using un-weighted operators
                                tmpMapping[localFields[linearjPointTLocal,iField]] = 場[iField,iT][jPoint]*maskingField[jPoint]
                            else
                                distance2 = distance2_point_to_box(conv.whole2model(jPoint),CartesianIndex(ones(Int, Ndimension-1)...), vec2car(ModelPoints[1:end-1]))
                                tmpMapping[localFields[linearjPointTLocal,iField]] = 場[iField,iT][jPoint]*maskingField[jPoint]*CerjanBoundaryCondition(distance2)
                            end
                        else
                            tmpMapping[localFields[linearjPointTLocal,iField]]=0.
                            #jPoint, νWhole[1],νWhole[end]
                        end
                    end
                end

                        

                        #場[iField,it]=string_as_varname(newstring, Array{Any,Ndimension-1}(undef,Tuple(wholeRegionPointsSpace)))
                        #
                        #νᶜtmpWholeMissing = ReplacerHorsLimiteParMissing(νᶜtmpWhole,PointsSpace[end])
                        # field values are defined only at the whole region and not at the Empty
                        #replace!(x -> x>0.2 ? missing : x, Array{Union{Float64, Missing}}(A) )
                        #replace!(x -> )
                        
                        
                
               
            end
            tmpAddress = nothing
            if testOnlyCentre
                tmpAddress = 1
            else
                tmpAddress=carDropDim(νRelative[iPoint])
            end
            costFunctions[iExpr,iTestFunctions]=substitute(semiSymbolicsOperators[tmpAddress,iExpr],tmpMapping)
            # be careful that semiSymbolicsOperators could be 2D
        end
    end

    #@show costFunctions

    #endregion

    return costFunctions,場,champsLimité

end

