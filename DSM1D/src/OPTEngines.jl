using Symbolics,UnPack

function findTXYZDependency(expression)
    # this subroutine gives the t-x-y-z dependency: e.g. expression(t,x,y) will gives [1,1,1,0]
    expressionDependency=ones(Int,4)
    for iDimension in 1:4
        if typeof(expand_derivatives(∂[iDimension](expression))==0) <: Bool 
            if expand_derivatives(∂[iDimension](expression))==0
                expressionDependency[iDimension] = 0
            end
        end
    end
    return expressionDependency
end

function CartesianOPTSymbolics(exprs,fields,vars, exts,extfields,extvars=nothing ; trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, pointsUsed=[3,3,3,3]),iDoSTEP7=false)
    # this is an outer function that demands the r.h.s. of the governing equation (external forces)
    # the user can use this function as : 
    #       i) CartesianOPTSymbolics(exprs,fields,vars,0,0)  for no external source
    #       ii) CartesianOPTSymbolics(exprs,fields,vars,exts,extfields) if the force is a simple function 
    #    or iii) CartesianOPTSymbolics(exprs,fields,vars, exts,extfields,extvars) for a force described with a differential expressions
    @variables C_dummy_dummy force_dummy
    if exts === 0
        exts =zeros(length(exprs))
        extfields = force_dummy
    end
    @show length(exprs)
    @show length(exts)
    if length(exprs) != length(exts)
        @error "the number of expressions is not equal to the number of external forces"
    end
    if extvars===nothing
        extvars = C_dummy_dummy
    end
    # Below if iDoSTEP7 is true, A will be a matrix expression but normally we do not need it (just to check with Geller & Takeuchi papers)
    A,neighbourVars,pointsUsedLeftAndRight,Δs=CartesianOPTSymbolics(exprs,fields,vars; trialFunctionsCharacteristics=trialFunctionsCharacteristics,iDoSTEP7=iDoSTEP7)
    g,neighbourVars_ext,pointsUsedLeftAndRight_ext,Δs_ext=CartesianOPTSymbolics(exts,extfields,extvars; trialFunctionsCharacteristics=trialFunctionsCharacteristics,iDoSTEP7=iDoSTEP7)
    return A,g,neighbourVars,pointsUsedLeftAndRight,Δs
end

function CartesianOPTSymbolics(exprs,fields,vars ; trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, pointsUsed=[3,3,3,3]), iDoSTEP7=false)

    #region General introduction

    # Nobuaki Fuji @ IPGP/UPC/IUF since 2024
    #
    # 
    # encouraged by Thibault Duretz @ U. Frankfurt Goethe, Kurama Okubo @ NIED
    #
    # with some inputs from Giacomo Aloisi @ ETH during Julia hackathon in the black forest October 2024
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
    #endregion

    #region STEP 0: initialising, unpacking etc. 

 

    @unpack orderBtime, orderBspace, pointsUsed = trialFunctionsCharacteristics

    NtypeofExpr=length(exprs)
    NtypeofMaterialVariables=length(vars)
    NtypeofFields=length(fields)
    

    #endregion

    #region STEP 1: investigation of all the fields and vars dependencies in terms of t-x-y-z

    variableDependency=ones(Int,4)
    fieldDependency=ones(Int,4)
    forceDependency=ones(Int,4)
    eachVariableDependency=ones(Int,4,NtypeofMaterialVariables) 
    eachFieldDependency=ones(Int,4,NtypeofFields)
  
    for iFields in 1:NtypeofFields
        eachFieldDependency[:,iFields]=findTXYZDependency(fields[iFields])
        fieldDependency = fieldDependency .* (ones(Int,4).-eachFieldDependency[:,iFields])
    end


    for iVars in 1:NtypeofMaterialVariables
        eachVariableDependency[:,iVars]=findTXYZDependency(vars[iVars])
        variableDependency = variableDependency .* (ones(Int,4).-eachVariableDependency[:,iVars])
    end
    fieldDependency = ones(Int,4).-fieldDependency
    variableDependency = ones(Int,4).-variableDependency
    forceDependency = ones(Int,4).-forceDependency
    # here we correct variableDependency and forceDependency with fieldDependency: if fieldDependency is zero then we do not take care of that dimension for the variables
    variableDependency = variableDependency .* fieldDependency
    forceDependency = forceDependency .* fieldDependency

    #endregion

    #region STEP 2: definition of points in time and space to be used

    # heaviside(x) = x > 0 ? 1 : x == 0 ? 0 : -1
    ### NF needs to work for a pure FD code too (orderBspline = -1)
    #
    #
    #
    #
    #
    ### NF needs to do this to complete the work!!


    # the orders of B-spline functions, depending on fields 
    orderBspline=zeros(Int,4)
    orderBspline[1]=orderBtime*fieldDependency[1]
    orderBspline[2:4]=orderBspace*fieldDependency[2:4]

    
    # the number of points used in the vicinity of the node, which is independent of the order of B-spline functions (see our paper)
    pointsUsedForFields=(pointsUsed.-1).*fieldDependency.+1
    
    # numbers of points to evaluate the integral for the governing equation filtered by the test functions
    
    # orderU is the maximal orders for the fields that we will use for OPT coefficients' exploration
    # orderU = (pointsUsedForFields.-1).*2 .+1
    # We don't do this anymore since we just need to go till the number of points used
    # orderU should be the same number as the number of orders considered
    orderU=pointsUsedForFields

    # The number of points to be used for time will need always one point after the present (for the future)

    Ltₗ = 0 
    Ltᵣ= 0
    if pointsUsedForFields[1]>1
        Ltₗ  = pointsUsedForFields[1]-2 # the leftmost point for t direction (past)
        Ltᵣ = 1 # the rightmost point (future)
    end
    # NF
    # The number of points to be used for space will depend on the boundary 
    # Below shoud be more flexible 

    Lxₗ  = (pointsUsedForFields[2])÷2 # the leftmost point for x direction
    Lxᵣ = (pointsUsedForFields[2]-1)÷2  # the rightmost point
    Lyₗ = (pointsUsedForFields[3])÷2 # the leftmost point for y direction
    Lyᵣ = (pointsUsedForFields[3]-1)÷2 # the rightmost point
    Lzₗ  = (pointsUsedForFields[4])÷2 # the leftmost point for z direction
    Lzᵣ = (pointsUsedForFields[4]-1)÷2 # the rightmost point    


    pointsUsedLeftAndRight=(Ltₗ=Ltₗ,Ltᵣ=Ltᵣ,Lxₗ=Lxₗ,Lxᵣ=Lxᵣ,Lyₗ=Lyₗ,Lyᵣ=Lyᵣ,Lzₗ=Lzₗ,Lzᵣ=Lzᵣ)

    # the number of points used to compute the deirvatives of material variables in the vicinity of the node
    pointsUsedForVariables=(pointsUsed.-1).*variableDependency.+1

    # orderC is the maximal orders for the structure variables that we use for inverting material variables for their derivatives
    orderC = pointsUsedForVariables

    # The number of points to be used for time will need always one point after the present (for the future)
    Stₗ = 0
    Stᵣ=0
    if pointsUsedForVariables[1]>1
        Stₗ  = pointsUsedForVariables[1]-2 # the leftmost point for t direction (past)
        Stᵣ = 1 # the rightmost point (future)
    end
    # NF
    # The number of points to be used for space will depend on the boundary 
    # Below shoud be more flexible 

    Sxₗ  = (pointsUsedForVariables[2])÷2 # the leftmost point for x direction
    Sxᵣ = (pointsUsedForVariables[2]-1)÷2  # the rightmost point
    Syₗ = (pointsUsedForVariables[3])÷2 # the leftmost point for y direction
    Syᵣ = (pointsUsedForVariables[3]-1)÷2 # the rightmost point
    Szₗ  = (pointsUsedForVariables[4])÷2 # the leftmost point for z direction
    Szᵣ = (pointsUsedForVariables[4]-1)÷2 # the rightmost point    

    




    # A generalised model of discretised (given) material variables
    C_given_generalised=Symbolics.variables(:C_given_generalised,1:Stₗ+Stᵣ+1,1:Sxₗ+Sxᵣ+1,1:Syₗ+Syᵣ+1,1:Szₗ+Szᵣ+1)
    
    #endregion

    #region STEP 3: definition of grids, inversion for material variable derivatives for general case

    print("STEP 3 started\n")

    # heterogeneous grids (not mesh!) should be adapted but for the moment it takes so much time 
    # NF -> maybe a subject for a hackathon

    Δt = Symbolics.variables(:Δt,1:Ltₗ+Ltᵣ)
    Δx = Symbolics.variables(:Δx,1:Lxₗ+Lxᵣ)
    Δy = Symbolics.variables(:Δy,1:Lyₗ+Lyᵣ)
    Δz = Symbolics.variables(:Δz,1:Lzₗ+Lzᵣ)

    # homogeneous grids

   @variables Δt₀ Δx₀ Δy₀ Δz₀

    Δt .= Δt₀
    Δx .= Δx₀
    Δy .= Δy₀
    Δz .= Δz₀

    #derivatives at the node
    @variables dt dx dy dz


   
    @variables C_deriv_generalised[1:orderC[1]*orderC[2]*orderC[3]*orderC[4]]

   
    # Generalised functions for C (material variables)

    C_Taylor_generalised(dt,dx,dy,dz)=sum(sum(sum(sum(C_deriv_generalised[it+(ix-1)*orderC[1]+(iy-1)*orderC[1]*orderC[2]+(iz-1)*orderC[1]*orderC[2]*orderC[3]]
    *dt^(it-1)*dx^(ix-1)*dy^(iy-1)*dz^(iz-1)
    // factorial(BigInt(it-1)) // factorial(BigInt(ix-1)) // factorial(BigInt(iy-1)) // factorial(BigInt(iz-1)) 
    for it in 1:orderC[1]; init=0) for ix in 1:orderC[2]; init=0) for iy in 1:orderC[3]; init=0) for iz in 1:orderC[4]; init=0)

    # Here we compute the derivatives of material variables with the maximum dimension.
    # The expression will then be concretised for each variable with different dependency

    eqnsC=Vector{Equation}(undef, 0)

    # NF exclusively tested this with left=right=1 or 0 for the moment 

    nodeZ = sum(Δz[iz] for iz in 1:Szₗ; init=0)
    for iz in 1:Szₗ+Szᵣ+1
        ddz = sum(Δz[jz] for jz in 1:iz-1; init=0)-nodeZ

        nodeY = sum(Δy[iy] for iy in 1:Lyₗ; init=0)
        for iy in 1:Syₗ+Syᵣ+1
            ddy = sum(Δy[jy] for jy in 1:iy-1; init=0)-nodeY

            nodeX = sum(Δx[ix] for ix in 1:Lxₗ; init=0)
            for ix in 1:Sxₗ+Sxᵣ+1
                ddx = sum(Δx[jx] for jx in 1:ix-1; init=0)-nodeX

                nodeT = sum(Δt[it] for it in 1:Ltₗ; init=0)
                for it in 1:Stₗ+Stᵣ+1
                    ddt = sum(Δt[jt] for jt in 1:it-1; init=0)-nodeT
                            
                    eqnsC = push!(eqnsC, C_given_generalised[it,ix,iy,iz]~C_Taylor_generalised(ddt,ddx,ddy,ddz))
                end
            end
        end
    end
    
    @time tmpvec=mySolvefor(eqnsC,C_deriv_generalised)
    newC_Coefs_deriv_generalised = reshape(tmpvec,(orderC[1],orderC[2],orderC[3],orderC[4])) # solve_for will give you a vector
    
    # Jacky's idea (and mine several in the spring 2024) is not to put this newC_Taylor_generalised expression during integral

    #endregion

    #region STEP 4: Taylor expansion for material variables and prepare mapping from partials of material variables to nodeVars
    Vars=[]
    neighbourVars=[]
    specificMaterialPartialMapping=Dict()

    for iVars in 1:NtypeofMaterialVariables
        #@show vars[iVars]
        specificMaterialMapping=Dict()
        newstring=split(string(vars[iVars]),"(")[1]
        nodeVars =Symbolics.variables(Symbol(newstring),1:(Stₗ+Stᵣ)*eachVariableDependency[1,iVars]+1,1:(Sxₗ+Sxᵣ)*eachVariableDependency[2,iVars]+1,1:(Syₗ+Syᵣ)*eachVariableDependency[3,iVars]+1,1:(Szₗ+Szᵣ)*eachVariableDependency[4,iVars]+1)

        newstring_for_partials="∂"*split(string(vars[iVars]),"(")[1]
        partialVars = Symbolics.variables(Symbol(newstring_for_partials),1:(orderC[1]-1)*eachVariableDependency[1,iVars]+1,1:(orderC[2]-1)*eachVariableDependency[2,iVars]+1,1:(orderC[3]-1)*eachVariableDependency[3,iVars]+1,1:(orderC[4]-1)*eachVariableDependency[4,iVars]+1)
        

        tmpVarExpression=sum(sum(sum(sum(partialVars[it,ix,iy,iz]*t^(it-1)*x^(ix-1)*y^(iy-1)*z^(iz-1)
    // factorial(BigInt(it-1)) // factorial(BigInt(ix-1)) // factorial(BigInt(iy-1)) // factorial(BigInt(iz-1)) 
    for it in 1:(orderC[1]-1)*eachVariableDependency[1,iVars]+1; init=0) for ix in 1:(orderC[2]-1)*eachVariableDependency[2,iVars]+1; init=0) for iy in 1:(orderC[3]-1)*eachVariableDependency[3,iVars]+1; init=0) for iz in 1:(orderC[4]-1)*eachVariableDependency[4,iVars]+1; init=0)

       
        for iz in 1:(Szₗ+Szᵣ)*eachVariableDependency[4,iVars]+1
            for iy in 1:(Syₗ+Syᵣ)*eachVariableDependency[3,iVars]+1
                for ix in 1:(Sxₗ+Sxᵣ)*eachVariableDependency[2,iVars]+1
                    for it in 1:(Stₗ+Stᵣ)*eachVariableDependency[1,iVars]+1
                        itt=it*eachVariableDependency[1,iVars]+(Stₗ+1)*(1-eachVariableDependency[1,iVars])
                        ixx=ix*eachVariableDependency[2,iVars]+(Sxₗ+1)*(1-eachVariableDependency[2,iVars])
                        iyy=iy*eachVariableDependency[3,iVars]+(Syₗ+1)*(1-eachVariableDependency[3,iVars])
                        izz=iz*eachVariableDependency[4,iVars]+(Szₗ+1)*(1-eachVariableDependency[4,iVars])
                        specificMaterialMapping[C_given_generalised[itt,ixx,iyy,izz]]=nodeVars[it,ix,iy,iz]
                    end
                end
             
            end
        end
        for iz in 1:(orderC[4]-1)*eachVariableDependency[4,iVars]+1
            for iy in 1:(orderC[3]-1)*eachVariableDependency[3,iVars]+1 
                for ix in 1:(orderC[2]-1)*eachVariableDependency[2,iVars]+1
                    for it in 1:(orderC[1]-1)*eachVariableDependency[1,iVars]+1
                        specificMaterialPartialMapping[partialVars[it,ix,iy,iz]]  = substitute(newC_Coefs_deriv_generalised[it,ix,iy,iz],specificMaterialMapping)
                    end
                end
            end
        end
        Vars =push!(Vars,tmpVarExpression)
        neighbourVars=push!(neighbourVars,nodeVars)
    end 

    
    #endregion

    #region STEP 5: translate the expressions with explicit variable names

    @variables U_deriv[1:orderU[1],1:orderU[2],1:orderU[3],1:orderU[4],1:NtypeofFields] # fields that we care of

    mapping=Dict()
    for iVar in 1:NtypeofMaterialVariables
        mapping[vars[iVar]]=Vars[iVar]
    end

    for iField in 1:NtypeofFields
        mapping[fields[iField]]=sum(sum(sum(sum(U_deriv[it,ix,iy,iz,iField]*t^(it-1)*x^(ix-1)*y^(iy-1)*z^(iz-1)
        // factorial(BigInt(it-1)) // factorial(BigInt(ix-1)) // factorial(BigInt(iy-1)) // factorial(BigInt(iz-1)) 
        for it in 1:orderU[1]; init=0) for ix in 1:orderU[2]; init=0) 
        for iy in 1:orderU[3]; init=0) for iz in 1:orderU[4]; init=0) 
    end

    #newExpr = expand_derivatives.(map((e) -> substitute(e, Dict(mapping)), exprs))
    newExpr = mySimplify.(map((e) -> substitute(e, Dict(mapping)), exprs))
   

    #endregion

    #region STEP 6: integrate the expressions over t, x, y, z

    print("STEP6 started\n")

    # orderBspline defines the shape of test functions: 0-> Dirac's delta (finite difference); 1-> triangle; 2+ -> B-spline curves

   

    nouvelleExpr=[]

    for iExpr in 1:NtypeofExpr
        variableArray=[t x y z]
        ΔVariableArray=[Δt₀ Δx₀ Δy₀ Δz₀]
       # @show deltaVariableArray=[Δt Δx Δy Δz] # we need think how to implement irregular grids
        newExpression = newExpr[iExpr]
        for iComponent in 1:4        
           
            if orderBspline[iComponent] == 0
                # orderBspline = 0 then it takes the integral is equivalent to the point 0
                variableArray[iComponent]
                newExpression=substitute(newExpression,Dict(variableArray[iComponent]=>0))
                #@show substitute(newExpr[iExpr], Dict(variableArray[iComponent]=>0))
            elseif orderBspline[iComponent] == 1
                # orderBspline = 1 then this integrates over the famous triangle
                # first we need to find the second-order antiderivative
                @time f1=integrateTaylorPolynomials(newExpression,variableArray[iComponent])
                @time f2=integrateTaylorPolynomials(f1,variableArray[iComponent])

                f2_right=substitute(f2,Dict([variableArray[iComponent]=>ΔVariableArray[iComponent]]))
                f2_left=substitute(f2,Dict([variableArray[iComponent]=>-ΔVariableArray[iComponent]]))

                ϕ_deriv_right=-1//ΔVariableArray[iComponent]
                ϕ_deriv_left = 1//ΔVariableArray[iComponent]
                desired_value=-(ϕ_deriv_right*f2_right-ϕ_deriv_left*f2_left)
                
                newExpression=mySimplify(desired_value)
              

            else
                error("high-order B-spline is not yet implemented, sorry.")
            end
            print("Integration over $(variableArray[iComponent]) is proceeded\n")
        end
        push!(nouvelleExpr,newExpression)
       
    end
    # here nouvelleExpr is the list of desired value (symbolic integration over the 'element')

    #endregion

    #region STEP 7: solve for equations (desired value ~ inner product) from lower derivatives

    if iDoSTEP7

        print("STEP 7 started\n")


        # The explicit optimally accuate operators
        A= Symbolics.variables(:A,1:Ltₗ+Ltᵣ+1,1:Lxₗ+Lxᵣ+1,1:Lyₗ+Lyᵣ+1,1:Lzₗ+Lzᵣ+1,1:NtypeofFields,1:NtypeofExpr)
        #
        innerProducts=[]

        for iExpr in 1:NtypeofExpr
            tmpInnerProduct=0
            for iField in 1:NtypeofFields
                for iz in 1:Lzₗ+Lzᵣ+1
                    dz = Δz₀*(iz-Lzₗ-1)
                    for iy in 1:Lyₗ+Lyᵣ+1
                        dy = Δy₀*(iy-Lyₗ-1)
                        for ix in 1:Lxₗ+Lxᵣ+1
                            dx = Δx₀*(ix-Lxₗ-1)
                            for it in 1:Ltₗ+Ltᵣ+1
                                dt = Δt₀*(it-Ltₗ-1)
                                tmpInnerProduct+=A[it,ix,iy,iz,iField,iExpr]*sum(sum(sum(sum(U_deriv[it,ix,iy,iz,iField]*dt^(it-1)*dx^(ix-1)*dy^(iy-1)*dz^(iz-1)
                                // factorial(BigInt(it-1)) // factorial(BigInt(ix-1)) // factorial(BigInt(iy-1)) // factorial(BigInt(iz-1)) 
                                for it in 1:orderU[1]; init=0) for ix in 1:orderU[2]; init=0) 
                                for iy in 1:orderU[3]; init=0) for iz in 1:orderU[4]; init=0) 
                            end
                        end
                    end
                end
            end
            push!(innerProducts,tmpInnerProduct)
        end

        hierarchical_equations=Vector{Equation}(undef, 0)

        for iExpr in 1:NtypeofExpr
            difference = mySimplify(innerProducts[iExpr]-nouvelleExpr[iExpr])

            for iField in 1:NtypeofFields
                for iz in 1:Lzₗ+Lzᵣ+1
                    for iy in 1:Lyₗ+Lyᵣ+1
                        for ix in 1:Lxₗ+Lxᵣ+1
                            for it in 1:Ltₗ+Ltᵣ+1
                                
                                tmpCoeffDiff = Symbolics.coeff(difference,U_deriv[it,ix,iy,iz,iField])
                                
                                push!(hierarchical_equations,tmpCoeffDiff~0)
                                
                            end
                        end
                    end
                end
            end
        end

        @time tmpvec=mySolvefor(hierarchical_equations,A)
        newA=reshape(tmpvec,(Ltₗ+Ltᵣ+1,Lxₗ+Lxᵣ+1,Lyₗ+Lyᵣ+1,Lzₗ+Lzᵣ+1,NtypeofFields,NtypeofExpr))
        newA=mySimplify(newA)
        #for iz in 1:Lzₗ+Lzᵣ+1
        #   for iy in 1:Lyₗ+Lyᵣ+1
        #        for ix in 1:Lxₗ+Lxᵣ+1
        #            for it in 1:Ltₗ+Ltᵣ+1
        #                @show it, ix, iy, iz,newA[it,ix,iy,iz,1,1]
        #            end
        #        end
        #    end
        #end
    end
   
    #endregion

    #region STEP 7bis: if we do not explicitly obtain the A matrix (which is cheaper and we can get what we want at the end)
    if !iDoSTEP7
        print("STEP 7bis started\n")
        finalExpr=Array{Num,1}(undef,NtypeofExpr)
        

        # we need to find the expressions for U_deriv in this case

        #Ufield= Symbolics.variables(:Ufield,1:Ltₗ+Ltᵣ+1,1:Lxₗ+Lxᵣ+1,1:Lyₗ+Lyᵣ+1,1:Lzₗ+Lzᵣ+1,1:NtypeofFields)
        Utmp=Symbolics.variables(:Utmp,1:Ltₗ+Ltᵣ+1,1:Lxₗ+Lxᵣ+1,1:Lyₗ+Lyᵣ+1,1:Lzₗ+Lzᵣ+1)
        Utmp_deriv= Symbolics.variables(:Utmp_deriv,1:Ltₗ+Ltᵣ+1,1:Lxₗ+Lxᵣ+1,1:Lyₗ+Lyᵣ+1,1:Lzₗ+Lzᵣ+1)
        eqnsC=Vector{Equation}(undef, 0)

    
        for iz in 1:Lzₗ+Lzᵣ+1
            dz = Δz₀*(iz-Lzₗ-1)
            for iy in 1:Lyₗ+Lyᵣ+1
                dy = Δy₀*(iy-Lyₗ-1)
                for ix in 1:Lxₗ+Lxᵣ+1
                    dx = Δx₀*(ix-Lxₗ-1)
                    for it in 1:Ltₗ+Ltᵣ+1
                        dt = Δt₀*(it-Ltₗ-1)
                        tmpValue=Utmp[it,ix,iy,iz]-sum(sum(sum(sum(Utmp_deriv[itt,ixx,iyy,izz]*dt^(itt-1)*dx^(ixx-1)*dy^(iyy-1)*dz^(izz-1)
                        // factorial(BigInt(itt-1)) // factorial(BigInt(ixx-1)) // factorial(BigInt(iyy-1)) // factorial(BigInt(izz-1)) 
                        for itt in 1:Ltₗ+Ltᵣ+1; init=0) for ixx in 1:Lxₗ+Lxᵣ+1; init=0) 
                        for iyy in 1:Lyₗ+Lyᵣ+1; init=0) for izz in 1:Lzₗ+Lzᵣ+1; init=0) 
                            eqnsC = push!(eqnsC,tmpValue~0)
                    end
                end
            end
        end

        tmpvec=mySolvefor(eqnsC,Utmp_deriv)
        newU_deriv = reshape(tmpvec,(Ltₗ+Ltᵣ+1,Lxₗ+Lxᵣ+1,Lyₗ+Lyᵣ+1,Lzₗ+Lzᵣ+1))

        specificFieldPartialMapping=Dict()

      
        for iField in 1:NtypeofFields

            tmpDictionary=Dict()
            newstring=split(string(fields[iField]),"(")[1]
            Ufield =Symbolics.variables(Symbol(newstring),1:Ltₗ+Ltᵣ+1,1:Lxₗ+Lxᵣ+1,1:Lyₗ+Lyᵣ+1,1:Lzₗ+Lzᵣ+1)

            for iz in 1:Lzₗ+Lzᵣ+1
                for iy in 1:Lyₗ+Lyᵣ+1
                    for ix in 1:Lxₗ+Lxᵣ+1
                        for it in 1:Ltₗ+Ltᵣ+1
                            tmpDictionary[Utmp[it,ix,iy,iz]]=Ufield[it,ix,iy,iz]
                        end
                    end
                end
            end


            for iz in 1:Lzₗ+Lzᵣ+1
                for iy in 1:Lyₗ+Lyᵣ+1
                    for ix in 1:Lxₗ+Lxᵣ+1
                        for it in 1:Ltₗ+Ltᵣ+1
                            specificFieldPartialMapping[U_deriv[it,ix,iy,iz,iField]]=substitute(newU_deriv[it,ix,iy,iz],tmpDictionary)
                        end
                    end
                end
            end
        end



        for iExpr in 1:NtypeofExpr
            tmpExpr=nouvelleExpr[iExpr]
            tmpExpr=substitute(tmpExpr,specificMaterialPartialMapping)
            tmpExpr=substitute(tmpExpr,specificFieldPartialMapping)
            finalExpr[iExpr]=mySimplify(tmpExpr)
        end

        Δs=(Δt₀=Δt₀,Δx₀=Δx₀,Δy₀=Δy₀,Δz₀=Δz₀)
        return finalExpr,neighbourVars,pointsUsedLeftAndRight,Δs 
    end
    #endregion
   
    #region STEP 8: translate the matrix with node material variables

    if iDoSTEP7
        print("STEP 8 started\n")
        #newA=substitute.(newA,specificMaterialPartialMapping)
        finalA=Array{Num,6}(undef,Ltₗ+Ltᵣ+1,Lxₗ+Lxᵣ+1,Lyₗ+Lyᵣ+1,Lzₗ+Lzᵣ+1,NtypeofFields,NtypeofExpr)
        for iExpr in 1:NtypeofExpr
            for iField in 1:NtypeofFields
                for iz in 1:Lzₗ+Lzᵣ+1
                    for iy in 1:Lyₗ+Lyᵣ+1
                        for ix in 1:Lxₗ+Lxᵣ+1
                            for it in 1:Ltₗ+Ltᵣ+1
                                tmpExpr=newA[it,ix,iy,iz,iField,iExpr]
                                tmpExpr=substitute(tmpExpr,specificMaterialPartialMapping)
                                finalA[it,ix,iy,iz,iField,iExpr]=mySimplify(tmpExpr)
                            end
                        end
                    end
                end
            end
        end
        Δs=(Δt₀=Δt₀,Δx₀=Δx₀,Δy₀=Δy₀,Δz₀=Δz₀)
        return finalA, neighbourVars, pointsUsedLeftAndRight,Δs
    end
    #endregion

end

function makeNumericalAuMoinsG(AuMoinsG,neighbourVars,fields,vars,pointsUsedLeftAndRight,ModelSizeTXYZ,modelingVariables,Δs)
    #region General introduction
    # This function is making the numerical vector with unknown fields to be 0
    # OPT does not need to explicitly derive A matrix so Au is given beforehand.
    #endregion

    #region STEP 1: analyse Au and dimensions
    NtypeofMaterialVariables=length(vars)
    NtypeofFields=length(fields)
    @unpack ModelSizeT, ModelSizeX, ModelSizeY, ModelSizeZ = ModelSizeTXYZ
    @unpack Ltₗ,Ltᵣ,Lxₗ,Lxᵣ,Lyₗ,Lyᵣ,Lzₗ,Lzᵣ=pointsUsedLeftAndRight
    @show size(AuMoinsG)
    #endregion

    # NF in fact, STEP 2 can already try to read the model parameters since there is no
    # reason to play again with the variable names, normally

    # here is the rule of the road: 
    # i) first define the 3D times time grids for fields
    # ii) then read the real values with input files
    # iii) make Au-g expression and give it to the fastdiff

    #region STEP 2: creating field and var vector (attention with the time marching part!)


        TXYZDependencies=[]
        #for iField in 1:NtypeofFields
        #    newstring=split(string(fields[iField]),"(")[1]
        #    tmpFields=Symbolics.variables(Symbol(newstring),1:ModelSizeX,1:ModelSizeY,1:ModelSizeZ)
        #    nodeFields=push!(nodeFields,tmpFields)
        #end
        nodeVarsModel=[]
        nodeVarsDepedencyConsideredModel=[]
        isVariablesTimeDependent=0
        for iVars in 1:NtypeofMaterialVariables
            newstring=split(string(vars[iVars]),"(")[1]*"M"
            TXYZDependency=findTXYZDependency(vars[iVars])
            TXYZDependencies=push!(TXYZDependencies,TXYZDependency)
            isVariablesTimeDependent+=TXYZDependency[1]
            # tmpVars is the parameters that we will get from input model or something (if constant, we give a scalar and not an array)
            tmpVars=Symbolics.variables(Symbol(newstring),1:(NpointsT-1)*TXYZDependency[1]+1,1:(ModelSizeX-1)*TXYZDependency[2]+1,1:(ModelSizeY-1)*TXYZDependency[3]+1,1:(ModelSizeZ-1)*TXYZDependency[4]+1)

            # then we create another matrix that should be distributed over the model nodes
            # Be careful of rearrangement of indices (time is the last one now)
            # Maybe this should be like that from the start 
            # but my mathematical philosophy stacked me for t-x-y-z order
            # But I agree that it should have been x-y-z-t order 
            tmpDistributedVars=Array{Num,4}(undef,ModelSizeX,ModelSizeY,ModelSizeZ,NpointsT)
            for it in 1:NpointsT
                for iz in 1:ModelSizeZ
                    for iy in 1:ModelSizeY
                        for ix in 1:ModelSizeX
                            itt=(it-1)*TXYZDependency[1]+1
                            ixx=(ix-1)*TXYZDependency[2]+1   
                            iyy=(iy-1)*TXYZDependency[3]+1
                            izz=(iz-1)*TXYZDependency[4]+1                    
                            tmpDistributedVars[ix,iy,iz,it]=tmpVars[itt,ixx,iyy,izz]
                        end
                    end
                end
            end

            nodeVarsModel=push!(nodeVarsModel,tmpDistributedVars)
            nodeVarsDepedencyConsideredModel=push!(nodeVarsDepedencyConsideredModel,tmpVars)
            # Be careful that time dependent variables are treated symbolically during the whole simulation!
        end


        #endregion

    #region STEP 3: read the model parameters
     #check the size of material variable size
     NtypeofMaterialVariables=length(nodeVarsModel)
     if length(materials)!=NtypeofMaterialVariables
         @error "ooops the input model does not have the same number of material variable types"
     end
end



# Below are the dinasaur functions (late 2024) that I think is not relevant any more

function ConcretisingAmatrix(Au,g, neighbourVars,pointsUsedLeftAndRight,fields,vars,ModelSizeTXYZ; PCscheme=true)
    #region General introduction
    # this function is extracting the OPT operators
    # before making a numerical matrix, this function still works only symbolically.
    # PCscheme = true will prepare A and δA in order to prepare predictor-corrector scheme in time marching
    #endregion

    #region STEP 1: analyse Au, correction of dimensions

    NtypeofMaterialVariables=length(vars)
    NpointsT, NpointsX, NpointsY, NpointsZ, NtypeofFields, NtypeofExpr = size(A)
    @unpack ModelSizeT, ModelSizeX, ModelSizeY, ModelSizeZ = ModelSizeTXYZ
    @unpack Ltₗ,Ltᵣ,Lxₗ,Lxᵣ,Lyₗ,Lyᵣ,Lzₗ,Lzᵣ=pointsUsedLeftAndRight
    if NpointsT == 1
        PCscheme=false
        ModelSizeT=1
    end
    if NpointsX == 1
        ModelSizeX=1
    end
    if NpointsY == 1
        ModelSizeY=1
    end
    if NpointsZ == 1
        ModelSizeZ=1
    end
    ModelSizeTXYZ=(ModelSizeT=ModelSizeT, ModelSizeX=ModelSizeX, ModelSizeY=ModelSizeY, ModelSizeZ=ModelSizeZ)
    #endregion
    
    #region STEP 2: creating field and var vector (attention with the time marching part!)
   
    
    TXYZDependencies=[]
    #for iField in 1:NtypeofFields
    #    newstring=split(string(fields[iField]),"(")[1]
    #    tmpFields=Symbolics.variables(Symbol(newstring),1:ModelSizeX,1:ModelSizeY,1:ModelSizeZ)
    #    nodeFields=push!(nodeFields,tmpFields)
    #end
    nodeVarsModel=[]
    nodeVarsDepedencyConsideredModel=[]
    isVariablesTimeDependent=0
    for iVars in 1:NtypeofMaterialVariables
        newstring=split(string(vars[iVars]),"(")[1]*"M"
        TXYZDependency=findTXYZDependency(vars[iVars])
        TXYZDependencies=push!(TXYZDependencies,TXYZDependency)
        isVariablesTimeDependent+=TXYZDependency[1]
        # tmpVars is the parameters that we will get from input model or something (if constant, we give a scalar and not an array)
        tmpVars=Symbolics.variables(Symbol(newstring),1:(NpointsT-1)*TXYZDependency[1]+1,1:(ModelSizeX-1)*TXYZDependency[2]+1,1:(ModelSizeY-1)*TXYZDependency[3]+1,1:(ModelSizeZ-1)*TXYZDependency[4]+1)

        # then we create another matrix that should be distributed over the model nodes
        # Be careful of rearrangement of indices (time is the last one now)
        # Maybe this should be like that from the start 
        # but my mathematical philosophy stacked me for t-x-y-z order
        # But I agree that it should have been x-y-z-t order 
        tmpDistributedVars=Array{Num,4}(undef,ModelSizeX,ModelSizeY,ModelSizeZ,NpointsT)
        for it in 1:NpointsT
            for iz in 1:ModelSizeZ
                for iy in 1:ModelSizeY
                    for ix in 1:ModelSizeX
                        itt=(it-1)*TXYZDependency[1]+1
                        ixx=(ix-1)*TXYZDependency[2]+1   
                        iyy=(iy-1)*TXYZDependency[3]+1
                        izz=(iz-1)*TXYZDependency[4]+1                    
                        tmpDistributedVars[ix,iy,iz,it]=tmpVars[itt,ixx,iyy,izz]
                    end
                end
            end
        end

        nodeVarsModel=push!(nodeVarsModel,tmpDistributedVars)
        nodeVarsDepedencyConsideredModel=push!(nodeVarsDepedencyConsideredModel,tmpVars)
        # Be careful that time dependent variables are treated symbolically during the whole simulation!
    end
  
   
    #endregion

    #region STEP 3: making A matrix for the past, the present and the future

    # Be careful! We are obliged to reorganise the A

    BigA=Array{Num,9}(undef, NpointsX, NpointsY, NpointsZ, NtypeofFields, NtypeofExpr, ModelSizeX, ModelSizeY, ModelSizeZ, NpointsT)
    

    for itModel in 1:NpointsT
        for izModel in 1:ModelSizeZ
            for iyModel in 1:ModelSizeY
                for ixModel in 1:ModelSizeX
                    # making map of variables
                    tmpDic=Dict()
                    
                    for iVars in 1:NtypeofMaterialVariables
                        for iNeighbourZ in 1:(NpointsZ-1)*TXYZDependencies[iVars][4]+1
                            for iNeighbourY in 1:(NpointsY-1)*TXYZDependencies[iVars][3]+1
                                for iNeighbourX in 1:(NpointsX-1)*TXYZDependencies[iVars][2]+1
                                    for iNeighbourT in 1:(NpointsT-1)*TXYZDependencies[iVars][1]+1
                                        generalisedNeighbourParam=neighbourVars[iVars][iNeighbourT,iNeighbourX,iNeighbourY,iNeighbourZ]
                                        #@show itModel, ixModel, iyModel, izModel, iNeighbourT,iNeighbourX,iNeighbourY,iNeighbourY
                                        #@show Ltₗ, Lxₗ, Lyₗ, Lzₗ
                                        # find the real coordinates
                                        itt=itModel+iNeighbourT-Ltₗ-1
                                        ixx=ixModel+iNeighbourX-Lxₗ-1
                                        iyy=iyModel+iNeighbourY-Lyₗ-1
                                        izz=izModel+iNeighbourZ-Lzₗ-1
                                        
                                        # check if it is not out of the boundary
                                        # Attention!!! Normally here we need to use another set of boundary condition operators
                                        #if 1<=itt<=NpointsT && 1<=ixx<=NpointsX && 1<=iyy<=NpointsY && 1<=izz<=NpointsZ
                                        #specifiedNeighbourParam=nodeVarsModel[iVars][ixx,iyy,izz,itt]
                                       
                                        if itt < 1
                                            itt = 1
                                        end
                                        if ixx < 1
                                            ixx = 1
                                        end
                                        if iyy < 1
                                            iyy = 1
                                        end
                                        if izz < 1
                                            izz = 1
                                        end
                                        if itt > NpointsT
                                            itt = NpointsT
                                        end
                                        if ixx > ModelSizeX
                                            ixx = ModelSizeX
                                        end
                                        if iyy > ModelSizeY
                                            iyy = ModelSizeY
                                        end
                                        if izz > ModelSizeZ
                                            izz = ModelSizeZ
                                        end
                                        #@show ixx,iyy,izz,itt
                                        specifiedNeighbourParam=nodeVarsModel[iVars][ixx,iyy,izz,itt]
                                
                                        tmpDic[generalisedNeighbourParam]=specifiedNeighbourParam
                                    end
                                end
                            end
                        end
                    end
                
                    for iExpr in 1:NtypeofExpr
                        for iField in 1:NtypeofFields
                            for iNeighbourZ in 1:NpointsZ
                                for iNeighbourY in 1:NpointsY
                                    for iNeighbourX in 1:NpointsX
                                        tmpExpr=A[itModel,iNeighbourX,iNeighbourY,iNeighbourZ,iField,iExpr]
                                        
                                        BigA[iNeighbourX, iNeighbourY, iNeighbourZ, iField, iExpr, ixModel, iyModel, izModel, itModel]=substitute(tmpExpr,tmpDic)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    #endregion

    #region STEP 4: PC scheme
    if PCscheme
        A₀=Array{Num,9}(undef, NpointsX, NpointsY, NpointsZ, NtypeofFields, NtypeofExpr, ModelSizeX, ModelSizeY, ModelSizeZ, NpointsT-1)
        δA=Array{Num,8}(undef, NpointsX, NpointsY, NpointsZ, NtypeofFields, NtypeofExpr, ModelSizeX, ModelSizeY, ModelSizeZ)
        futureOnepoint=Array{Num,5}(undef, NtypeofFields, NtypeofExpr, ModelSizeX, ModelSizeY, ModelSizeZ)
        #A₀[:,:,:,:,:,:,:,:,1:NpointsT-1]=BigA[:,:,:,:,:,:,:,:,1:NpointsT-1]
        #δA[:,:,:,:,:,:,:,:]=BigA[:,:,:,:,:,:,:,:,NpointsT]
        #δA[Lxₗ+1,Lyₗ+1,Lzₗ+1,:,:,:,:,:]=0
        #futureOnepoint=BigA[Lxₗ+1,Lyₗ+1,Lzₗ+1,:,:,:,:,:,NpointsT]
        A₀=BigA
        δA=nothing
        futureOnepoint=nothing
    else
        A₀=BigA
        δA=nothing
        futureOnepoint=nothing
    end
    #endregion
    return A₀,δA,futureOnepoint,ModelSizeTXYZ,nodeVarsDepedencyConsideredModel
end

function makeNumericalAmatrixExplicit(A₀,δA,futureOnepoint,ModelSizeTXYZ,nodeVarsModel,modelingVariables,Δs)
    
    #region General introduction
    # Finally this function will substitute all the symbols (for the time-dependent material variables, we need to call this every time step)
    # 
    # for the T-dependent variables, we need to be careful because here materials, nodeVarsModel are in the order of XYZT (not TXYZ)

    #endregion
      
    #region STEP 1: Check the configuration, check the dimensions and copy material variables

    @unpack ModelSizeT, ModelSizeX, ModelSizeY, ModelSizeZ = ModelSizeTXYZ
    @unpack materials,Δnodes = modelingVariables
    @unpack ΔnodeX,ΔnodeY,ΔnodeZ,ΔnodeT=Δnodes
    @unpack Δt₀,Δx₀,Δy₀,Δz₀=Δs

    NpointsT=size(A₀)[9]
    NtypeofFields=size(A₀)[4]
    NtypeofExpr=size(A₀)[5]
    
    PCscheme=true
    if futureOnepoint === nothing
        PCscheme=false
    end

    #check the size of material variable size
    NtypeofMaterialVariables=length(nodeVarsModel)
    if length(materials)!=NtypeofMaterialVariables
        @error "ooops the input model does not have the same number of material variable types"
    end
    
    disretisedMaterials=copy(nodeVarsModel) # this should have 4D 

    bigDic=Dict()
    for iVar in 1:NtypeofMaterialVariables
      
        if length(nodeVarsModel[iVar]) != length(materials[iVar])
            @error "ooops the input model seems incompatible with the proposed operators"
        end
        disretisedMaterials[iVar] = reshape(materials[iVar],size(nodeVarsModel[iVar]))


        nx = size(nodeVarsModel[iVar])[1]
        ny = size(nodeVarsModel[iVar])[2]
        nz = size(nodeVarsModel[iVar])[3]
        nt = size(nodeVarsModel[iVar])[4]

        for it in 1:nt
            for iz in 1:nz
                for iy in 1:ny
                    for ix in 1:nx
                        bigDic[nodeVarsModel[iVar][ix,iy,iz,it]]=disretisedMaterials[iVar][ix,iy,iz,it]
                    end
                end
            end
        end


    end


    #endregion
    
    #region STEP 2: substitute Δnodes
    
    bigDic[Δx₀] = ΔnodeX
    bigDic[Δy₀] = ΔnodeY
    bigDic[Δz₀] = ΔnodeZ
    bigDic[Δt₀] = ΔnodeT
  
    @show substitute(A₀,bigDic)
    #endregion

    #region STEP 3: make a cost function (in fact the next generation of this code should pass without the construction of A)

    #BigA=Array{Num,9}(undef, NpointsX, NpointsY, NpointsZ, NtypeofFields, NtypeofExpr, ModelSizeX, ModelSizeY, ModelSizeZ, NpointsT)
    
    nodeFields=Array{Any,5}(undef, ModelSizeX, ModelSizeY, ModelSizeZ, NtypeofFields, NpointsT)

    # nodeFields[:,:,:,:,1] should be the future and nodeFields[:,:,:,:,2:NpointsT] are the past

    #endregion


end