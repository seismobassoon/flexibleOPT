
include("../src/batchSymbolics.jl")
include("../src/batchUseful.jl")
using BenchmarkTools

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

function CartesianOPTSymbolics(exprs,fields,vars, exts = nothing ; trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, pointsUsed=[3,3,3,3]))

    #region General introduction

    # Nobuaki Fuji @ IPGP
    # with some inputs from Giacomo Aloisi @ ETH during Julia hackathon in the black forest

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

    if exts === nothing
        exts = 0
    end

    @unpack orderBtime, orderBspace, pointsUsed = trialFunctionsCharacteristics

    NtypeofExpr=length(exprs)
    NtypeofMaterialVariables=length(vars)
    NtypeofFields=length(fields)
    NtypeofForces=length(exts)

    #endregion

    #region STEP 1: investigation of all the fields and vars dependencies in terms of t-x-y-z

    variableDependency=ones(Int,4)
    fieldDependency=ones(Int,4)
    forceDependency=ones(Int,4)
    eachVariableDependency=ones(Int,4,NtypeofMaterialVariables) 
    eachFieldDependency=ones(Int,4,NtypeofFields)
    eachForceDependency=ones(Int,4,NtypeofForces)

    for iFields in 1:NtypeofFields
        eachFieldDependency[:,iFields]=findTXYZDependency(fields[iFields])
        fieldDependency = fieldDependency .* (ones(Int,4).-eachFieldDependency[:,iFields])
    end

    for iForces in 1:NtypeofForces
        eachForceDependency[:,iForces]=findTXYZDependency(exts[iForces])
        forceDependency=forceDependency .* (ones(Int,4).-eachForceDependency[:,iForces])
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
    orderU = (pointsUsedForFields.-1).*2 .+1

    # The number of points to be used for time will need always one point after the present (for the future)

    Ltₗ  = pointsUsedForFields[1]-1-pointsUsedForFields[1]%2 # the leftmost point for t direction (past)
    Ltᵣ = pointsUsedForFields[1]%2 # the rightmost point (future)

    # NF
    # The number of points to be used for space will depend on the boundary 
    # Below shoud be more flexible 

    Lxₗ  = (pointsUsedForFields[2])÷2 # the leftmost point for x direction
    Lxᵣ = (pointsUsedForFields[2]-1)÷2  # the rightmost point
    Lyₗ = (pointsUsedForFields[3])÷2 # the leftmost point for y direction
    Lyᵣ = (pointsUsedForFields[3]-1)÷2 # the rightmost point
    Lzₗ  = (pointsUsedForFields[4])÷2 # the leftmost point for z direction
    Lzᵣ = (pointsUsedForFields[4]-1)÷2 # the rightmost point    


    # the number of points used to compute the deirvatives of material variables in the vicinity of the node
    pointsUsedForVariables=(pointsUsed.-1).*variableDependency.+1

    # orderC is the maximal orders for the structure variables that we use for inverting material variables for their derivatives
    orderC = pointsUsedForVariables

    # The number of points to be used for time will need always one point after the present (for the future)

    Stₗ  = pointsUsedForVariables[1]-1-pointsUsedForVariables[1]÷2  # the leftmost point for t direction (past)
    Stᵣ = pointsUsedForVariables[1]÷2  # the rightmost point (future)

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
    / factorial(BigInt(it-1)) / factorial(BigInt(ix-1)) / factorial(BigInt(iy-1)) / factorial(BigInt(iz-1)) 
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
    specificMaterialPartialMapping=Dict()

    for iVars in 1:NtypeofMaterialVariables
    
        specificMaterialMapping=Dict()
        newstring=split(string(vars[iVars]),"(")[1]
        nodeVars =Symbolics.variables(Symbol(newstring),1:(Stₗ+Stᵣ)*eachVariableDependency[1,iVars]+1,1:(Sxₗ+Sxᵣ)*eachVariableDependency[2,iVars]+1,1:(Syₗ+Syᵣ)*eachVariableDependency[3,iVars]+1,1:(Szₗ+Szᵣ)*eachVariableDependency[4,iVars]+1)

        newstring_for_partials="∂"*split(string(vars[iVars]),"(")[1]
        partialVars = Symbolics.variables(Symbol(newstring_for_partials),1:(orderC[1]-1)*eachVariableDependency[1,iVars]+1,1:(orderC[2]-1)*eachVariableDependency[2,iVars]+1,1:(orderC[3]-1)*eachVariableDependency[3,iVars]+1,1:(orderC[4]-1)*eachVariableDependency[4,iVars]+1)
        

        tmpVarExpression=sum(sum(sum(sum(partialVars[it,ix,iy,iz]*t^(it-1)*x^(ix-1)*y^(iy-1)*z^(iz-1)
    / factorial(BigInt(it-1)) / factorial(BigInt(ix-1)) / factorial(BigInt(iy-1)) / factorial(BigInt(iz-1)) 
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
        / factorial(BigInt(it-1)) / factorial(BigInt(ix-1)) / factorial(BigInt(iy-1)) / factorial(BigInt(iz-1)) 
        for it in 1:orderU[1]; init=0) for ix in 1:orderU[2]; init=0) 
        for iy in 1:orderU[3]; init=0) for iz in 1:orderU[4]; init=0) 
    end

    newExpr = expand_derivatives.(map((e) -> substitute(e, Dict(mapping)), exprs))

   

    #endregion

    #region STEP 6: integrate the expressions over t, x, y, z

    print("STEP6 started\n")

    # orderBspline defines the shape of test functions: 0-> Dirac's delta (finite difference); 1-> triangle; 2+ -> B-spline curves

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
                            / factorial(BigInt(it-1)) / factorial(BigInt(ix-1)) / factorial(BigInt(iy-1)) / factorial(BigInt(iz-1)) 
                            for it in 1:orderU[1]; init=0) for ix in 1:orderU[2]; init=0) 
                            for iy in 1:orderU[3]; init=0) for iz in 1:orderU[4]; init=0) 
                        end
                    end
                end
            end
        end
        push!(innerProducts,tmpInnerProduct)
    end

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

    print("STEP 7 started\n")

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

   
    #endregion

    #region STEP 8 : translate the matrix with node material variables
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

    return finalA
    #endregion
end



# here we define the governing equation (the dependency should be explicitly declared with parentheses!!)
#@variables ρ(x) μ(x,t) u(x,y,t)

iExecute = 1

if iExecute == 0
    @variables u(x)
    expr = ∂x²(u)

    exprs = expr
    vars = 1
    exts = 0
    fields = u

    CartesianOPTSymbolics(exprs,fields,vars)
elseif iExecute == 1 

    # 1D SH in freq. problem
    @variables ρ μ  ω u(x) f(x)
    expr = ρ*ω^2*u + μ*∂x²(u)

    # declaration of physics
    exprs = expr
    vars = ρ, μ, ω
    exts = 0
    fields = u # for the moment we need to put like ux uy uz to distinguish different components



elseif iExecute == 2

    # 1D SH in time problem
    @variables ρ(x) μ(x) u(t,x) f(t,x) 
    expr = ρ*∂t²(u)- ∂x(μ*∂x(u)) # scalar PDE(s) to be solved

    # declaration of physics
    exprs = expr
    vars = ρ, μ
    exts = 0
    fields = u # for the moment we need to put like ux uy uz to distinguish different components



elseif iExecute == 3

    # 2D Poisson homo

    @variables κ T(x,y)
    expr = ∂x²(κ*T)+∂y²(κ*T)

    exprs = expr
    vars = κ
    exts=0

    fields = T
 

elseif iExecute == 4

    # 2D Poisson hetero

    @variables κ(x,y) T(x,y)
    #expr = ∂x²(κ*T)+∂y²(κ*T)
    expr = ∂x(κ*∂x(T))+∂y(κ*∂y(T))
    exprs = mySimplify(expr)
    vars = κ
    exts=0
    fields = T
  

elseif iExecute == 5

    # 1D Poisson hetero
    @variables κ(x) T(x)
    expr = ∂x²(κ*T)

    exprs = mySimplify(expr)
    vars = κ
    exts = 0
    fields = T

    
end

A=CartesianOPTSymbolics(exprs,fields,vars)  # or CartesianOPTSymbolics(exprs,fields,vars,exts)

