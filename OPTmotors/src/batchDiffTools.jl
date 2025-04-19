using SparseDiffTools,SparseArrays,Symbolics

function Residual!(F,1dcostfunctions,unknownField,knownField,knownForce)
    for i in eachindex(F)
        mapping = Dict()
        for j in eachindex(knownField)
            mapping[symbKnownField[j]] = knownField[j]
        end
        for j in eachindex(knownForce)
            mapping[symbKnownForce[j]] = knownForce[j]
        end
        F[i] = substitute(1dcostfunctions[i],mapping)
    end
    return
end

function sparseColouring(1dcostfunctions,unknownField,knownField,knownForce)
    nCostfunctions = length(1dcostfunctions)
    nUnknownField = length(unknownField)
    input = rand(nUnknownField)
    output = rand(nCostfunctions)
    F=zeros(nCostfunctions)
    Res_closed! = (F,U) -> Residual!(F,1dcostfunctions,unknownField,knownField,knownForce)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))
    colors      = matrix_colors(J)
    return J, colors
end


function myForwarddiff_color_jacobian!(F, U,nEq)


    # Sparsity pattern
    input       = rand(ncx)
    output      = similar(input)
    Res_closed! = (F, U) -> Residual!(F, U, U0, U00, f, G, Gc, ρ, Δx, Δt, x, t)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed!, output, input)
    J           = Float64.(sparse(sparsity))

    # Makes coloring
    colors      = matrix_colors(J)

    # normalisation by the number of equations
    normalisation = 1.0/nEq
    r1 = 1.0
    for iter in 1:10
        Res_closed! = (F, U) -> Residual!(F, U)
        r = norm(F)*normalisation
        
        if iter==1 r1 = r; end
       
        if r/r1 < 1e-8 break end
            
        # Jacobian assembly
        forwarddiff_color_jacobian!(J, Res_closed!, U, colorvec = colors)

        # Solve
        δU   .= .-J\F

        # update
        U    .+= δU
    end

    return 
end