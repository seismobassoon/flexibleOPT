using SparseDiffTools,SparseArrays,Symbolics,Enzyme

include("../src/batchUseful.jl")

function buildNumericalFunctions(costfunctions, symbUnknownField, symbKnownField, symbKnownForce)
    # Collect all known inputs (constants) and unknown inputs (variables)
    knownInputs = vcat(reduce(vcat, reduce(vcat, symbKnownField)), reduce(vcat, symbKnownForce))
    unknownInputs = reduce(vcat, symbUnknownField)
    all_inputs = vcat(unknownInputs, knownInputs)
    
    residual_func = Array{Function, 1}(undef, length(costfunctions))
    
    for i in eachindex(costfunctions)
        # Create the symbolic function
        residual_func_expr = build_function(costfunctions[i], all_inputs; expression = Val{false})
        
        # Evaluate the symbolic function and store it as a numerical function
        residual_func[i] = eval(residual_func_expr)
    end
    #@show residual_func[120]
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
    knownInputs = vcat(reduce(vcat,reduce(vcat,knownField) ),reduce(vcat,reduce(vcat,knownForce)))
    unknownInputs = reduce(vcat,unknownField)
    #all_inputs = vcat(unknownInputs,knownInputs)
    #@show length(knownInputs),length(unknownInputs)
    return unknownInputs,knownInputs
end

function Residual!(F,f,unknownInputs,knownInputs)
    all_inputs=vcat(unknownInputs,knownInputs)
    for i in eachindex(f)
        F[i]=f[i]((all_inputs))
    end
end
function temporaryConstantResidualFunction(f::Vector{F}, U::AbstractVector, knownInputs::AbstractVector) where {F <: Function}
    N_U = length(U)
    N_K = length(knownInputs)
    function f_specific!(Fout::AbstractVector{T}, unknownInputs::AbstractVector{T}) where {T}
        all_inputs = Vector{T}(undef, N_U + N_K)
        @inbounds begin
            all_inputs[1:N_U] .= unknownInputs
            for i in 1:N_K
                all_inputs[N_U + i] = T(knownInputs[i])
            end
            for i in eachindex(f)
                Fout[i] = f[i](all_inputs)
            end
        end
        return nothing
    end
    return f_specific!
end


function sparseColouring(f,unknownField,knownField,knownForce)
    # this function is colouring the matrix
    nCostfunctions = length(f)
    nUnknownField = length(unknownField)
    #input = Vector{Float64}(undef,nUnknownField)
    input = rand(nUnknownField)
    #output = Vector{Float64}(undef,nCostfunctions)
    output = zeros(nCostfunctions)
    F=zeros(nCostfunctions)

    U,knownInputs = makeInputsForNumericalFunctions(input,knownField,knownForce)

    Res_closed_look! = (F,U) -> Residual!(F,f,U,knownInputs)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed_look!,output,U)
    J           = Float64.(sparse(sparsity))
    V=rand(length(U))
    cache=ForwardColorJacCache(Res_closed_look!, V)
    #J = (sparse(sparsity))
    #colors      = matrix_colors(J)
    #coloring_cache = ForwardColorJacCache(Res_closed_look!, J, U; sparsity = J)
    return J,cache
end


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


    
    # normalisation by the number of equations
    normalisation = 1.0/nEq
    r1 = 1.0
    unknownField .= 0.0
    U,knownInputs = makeInputsForNumericalFunctions(unknownField,knownField,knownForce)
    δU = U
    F=zeros(nEq)
    f_specific! = temporaryConstantResidualFunction(f,U,knownInputs)
    for iter in 1:nIteration
        
        #Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
        #Res_closed! = (F,U) -> Residual!(F,f,U,knownInputs)
        #Res_closed!(F,U)
        f_specific!(F,U)
        @show r = norm(F)
        
        if iter==1 r1 = r; end
        if r === 0.0 break end
        if r/r1 < smallNumber break end
      
            

        # coloring
        #V=rand(length(U))
        #cache=ForwardColorJacCache(Res_closed!, V)
        #@show V, cache
        # Jacobian assembly
    
        @time forwarddiff_color_jacobian!(J, f_specific!, U, cache)
        #@time handMadeJacobianComputation!(J, f_specific!, F, U, rows, cols)

        # Solve
        #@time factor = lu(J)  # Or try `ldlt`, `cholesky`, or `qr` depending on J's properties
        #invJac=inv(factor)
        invJac = myInv(J)

        @time δU = - invJac * F
        #@time δU   .= .-J\F

        # update
        α = 1.0
        U    .+= α .* δU
    end

    unknownField = reshape(U,pointsFieldSpace)
    return unknownField
end


function handMadeJacobianComputation!(J, f_specific!, F, U, rows, cols)
    # Define the residual wrapper function to give to Enzyme
    #residual_wrapper(F::Vector{Float64}, U::Vector{Float64}) = f_specific(F, U)

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
        #dF = zeros(length(F))
        
        # Call Enzyme's forward-mode autodiff on the wrapper function
        #Enzyme.autodiff(Enzyme.Forward,residual_wrapper,Const, Duplicated(F, dF), Duplicated(U, dU))
        newU = U .+ dU
        newF = similar(F)
        f_specific!(newF, newU)
        dF_i = newF[i]-F[i]

        # Store entry in the sparse Jacobian
        J[i, j] = dF_i/dU[j]
    end
end