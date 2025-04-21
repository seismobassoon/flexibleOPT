


function timeStepOptimisation!(f,unknownField,knownField,knownForce,J,cache,pointsFieldSpace;nIteration=10,smallNumber =1.e-8)


    nEq = length(f)
   
    nUnknownField = length(unknownField)
    #input = Vector{Float64}(undef,nUnknownField)
    input = rand(nUnknownField)
    #output = Vector{Float64}(undef,nCostfunctions)
    output = zeros(nEq)
    F=zeros(nEq)

    U,knownInputs = makeInputsForNumericalFunctions(input,knownField,knownForce)

    Res_closed_look! = (F,U) -> Residual!(F,f,U,knownInputs)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed_look!,output,U)
    rows, cols, _ = findnz(sparse(sparsity))

    J = spzeros(length(F), length(U))
   

    # start the simu

    # normalisation by the number of equations
    normalisation = 1.0/nEq
    r1 = 1.0
    unknownField .= 0.0
    U,knownInputs = makeInputsForNumericalFunctions(unknownField,knownField,knownForce)
    δU = U
    @show typeof(f)
    
    #f_specific = temporaryConstantResidualFunction(f,U,knownInputs)
    #@show f_specific = make_static_residual_generator(f,knownInputs)
    #residual_wrapper(F::Vector{Float64}, U::Vector{Float64}) = f_specific(F, U)
    f_tuple = tuple(f...)  # convert vector to tuple
    f_specific = StaticResidual(f_tuple, knownInputs)
    #f_specific = StaticResidual(f, knownInputs)

    for iter in 1:nIteration
        #f_specific = (F, U) -> Residual!(F, f, U, knownInputs)
        #f_specific = StaticResidual(f, knownInputs)
        #Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
        #Res_closed! = (F,U) -> Residual!(F,f,U,knownInputs)
        #Res_closed!(F,U)
        f_specific(F,U)
        @show r = norm(F)*normalisation
        
        if iter==1 r1 = r; end
        if r === 0.0 break end
        if r/r1 < smallNumber break end
      

        compute_jacobian!(J, f_specific, F, U, rows, cols)

        # coloring
        #V=rand(length(U))
        #cache=ForwardColorJacCache(Res_closed!, V)
        #@show V, cache
        # Jacobian assembly
        
        # Solve
        @time factor = lu(J)  # Or try `ldlt`, `cholesky`, or `qr` depending on J's properties
        @time δU .= - (factor \ F)
        #@time δU   .= .-J\F

        # update
        α = 0.5
        U    .+= α .* δU
    end

    unknownField = reshape(U,pointsFieldSpace)
    return unknownField
end


function compute_jacobian!(J, f_specific, F, U, rows, cols)
    # Define the residual wrapper function to give to Enzyme
    residual_wrapper(F::Vector{Float64}, U::Vector{Float64}) = f_specific(F, U)

    # perturbation 
    ΔU = maximum((1.e-8,maximum(U)))

    # Loop through the non-zero entries
    for k in eachindex(rows)
        i = rows[k]
        j = cols[k]

        # Perturbation vector
        dU = zeros(length(U))
        dU[j] = ΔU

        # dF: output for directional derivative
        dF = zeros(length(F))
        
        # Call Enzyme's forward-mode autodiff on the wrapper function
        Enzyme.autodiff(Enzyme.Forward,residual_wrapper,Const, Duplicated(F, dF), Duplicated(U, dU))


        # Store entry in the sparse Jacobian
        J[i, j] = dF[i]/dU[j]
    end
end

struct StaticResidual{F, K}
    f::F  # Tuple of concrete function types
    known::K
end

StaticResidual(f::Tuple, known::Vector{Float64}) = StaticResidual{typeof(f), typeof(known)}(f, known)
function (sr::StaticResidual)(F::Vector{Float64}, U::Vector{Float64})
    all_inputs = vcat(U, sr.known)
    for i in eachindex(sr.f)
        F[i] = sr.f[i](all_inputs)
    end
    return nothing
end
