using SparseDiffTools,SparseArrays

function sparseColouring()
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