using SparseDiffTools,SparseArrays,Symbolics

function Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)

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

function sparseColouring(costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
    nCostfunctions = length(costfunctions)
    nUnknownField = length(unknownField)
    input = Vector{Float64}(undef,nUnknownField)
    output = Vector{Float64}(undef,nCostfunctions)
    F=zeros(nCostfunctions)
    Res_closed! = (F,unknownField) -> Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!,output, input)
    J           = Float64.(sparse(sparsity))
    colors      = matrix_colors(J)
    return J,colors
end


function timeStepOptimisation!(F, costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce,J,colors;nIteration=10,smallNumber =1.e-8)
    @show knownForce
    nEq = length(costfunctions)
    # normalisation by the number of equations
    normalisation = 1.0/nEq
    r1 = 1.0
    #unknownField .= 0.0
    for iter in 1:nIteration
        unknownField.=0.0
        Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
        Res_closed! = (F,unknownField) -> Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
    
  
        r = norm(F)*normalisation
        
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