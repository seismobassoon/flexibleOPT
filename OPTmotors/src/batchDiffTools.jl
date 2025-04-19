using SparseDiffTools,SparseArrays,Symbolics


function buildNumericalFunctions(costfunctions,symbUnknownField,symbKnownField,symbKnownForce)
    # this function will translate the costfunctions fully numerically
    knownInputs = vcat(reduce(vcat,reduce(vcat,symbKnownField) ),reduce(vcat,symbKnownForce))
    unknownInputs = reduce(vcat,symbUnknownField)
    all_inputs = vcat(unknownInputs,knownInputs)
    residual_func=Array{Any,1}(undef,length(costfunctions))
    for i in eachindex(costfunctions)
        residual_func_expr = build_function(costfunctions[i], all_inputs; expression = Val{false})
        residual_func[i] = eval(residual_func_expr)
    end
    return residual_func

    
end


function Residual_OLD!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)

    mapping = Dict()
    
    for k in eachindex(knownField)
        for j in eachindex(knownField[k])
            #mapping[symbKnownField[k][j]] = knownField[k][j]
        end
    end

    for j in eachindex(unknownField)
        #mapping[symbUnknownField[j]] = unknownField[j]
    end

    for j in eachindex(knownForce)
        #mapping[symbKnownForce[j]] = knownForce[j]
    end

    for i in eachindex(F)
        F[i] = substitute(costfunctions[i],mapping)
    end
    return
end

function makeInputsForNumericalFunctions(unknownField,knownField,knownForce)
    knownInputs = vcat(reduce(vcat,reduce(vcat,knownField) ),reduce(vcat,knownForce))
    unknownInputs = reduce(vcat,unknownField)
    all_inputs = vcat(unknownInputs,knownInputs)
    return all_inputs
end

function Residual!(F,f,unknownField,knownField,knownForce)
    all_inputs=makeInputsForNumericalFunctions(unknownField,knownField,knownForce)
    for i in eachindex(f)
        @show F[i]=f[i](all_inputs...)
    end
end


function sparseColouring(f,unknownField,knownField,knownForce)
    nCostfunctions = length(f)
    nUnknownField = length(unknownField)
    #input = Vector{Float64}(undef,nUnknownField)
    input = rand(nUnknownField)
    #output = Vector{Float64}(undef,nCostfunctions)
    output = rand(nCostfunctions)
    F=zeros(nCostfunctions)
    Res_closed! = (F,unknownField) -> Residual!(F,f,unknownField,knownField,knownForce)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!,output, input)
    J           = Float64.(sparse(sparsity))
    colors      = matrix_colors(J)
    return J,colors
end


function timeStepOptimisation!(F, f,unknownField,knownField,knownForce,J,colors;nIteration=10,smallNumber =1.e-8)

    nEq = length(f)
    # normalisation by the number of equations
    normalisation = 1.0/nEq
    r1 = 1.0
    #unknownField .= 0.0
    for iter in 1:nIteration
        
        #Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
        Res_closed! = (F,unknownField) -> Residual!(F,f,unknownField,knownField,knownForce)
        Res_closed!(F,unknownField)
  
        

        @show r = norm(F)#*normalisation
        
        if iter==1 r1 = r; end
        if r === 0.0 break end
        if r/r1 < smallNumber break end
      
            
        # Jacobian assembly
        forwarddiff_color_jacobian!(J, Res_closed!, unknownField, colorvec = colors)

        # Solve
        δunknownField   .= .-J\F

        # update
        unknownField    .+= δunknownField
    end

    return unknownField
end