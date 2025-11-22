using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")

include("../src/imageReader.jl") # read 2D images for models

include("../src/OPTwrappers.jl") 

using Enzyme
using SparseArrays


function detect_backend()
    _has_cuda = false
    _has_metal = false

    try
        if Sys.isapple()
            @eval using Metal
            devs = Metal.devices()
            _has_metal = !isempty(devs)
        else
            @eval using CUDA
            _has_cuda = CUDA.has_cuda()
        end
    catch e
        @warn "GPU backend not available: $e"
    end

    return (_has_cuda, _has_metal)
end

_has_cuda, _has_metal = detect_backend()


"""
    enzyme_column_jacobian!(J, f!, u; backend=:cpu)

Computes the sparse Jacobian J of f!(res, u) by Enzyme autodiff, column-by-column.

# Arguments
- `J`: preallocated SparseMatrixCSC (m × n)
- `f!`: in-place residual function f!(res, u)
- `u`: input vector (Array, CuArray, MtlArray)
- `backend`: `:cpu`, `:cuda`, or `:metal`

# Notes
- The result is stored in `J`
- `f!(res, u)` must support the selected backend array type
"""
function enzyme_column_jacobian!(J::SparseMatrixCSC{T}, f!::Function, u::AbstractVector{T}; backend::Symbol = :cpu) where T
    n = length(u)
    tmp_u = copy(u)
    
    res = similar(u)  # for shape
    f!(res, u)
    m = length(res)

    tmp_res = similar(res)
    df = similar(res)
    du = similar(u)

    for j in 1:n
        fill!(du, zero(T))
        du[j] = one(T)
        fill!(df, zero(T))

        Enzyme.autodiff(
            Enzyme.Reverse,
            f!,
            Duplicated(df, tmp_res),
            Duplicated(du, tmp_u)
        )

        for i in 1:m
            if df[i] != 0
                J[i, j] = df[i]
            end
        end
    end

    return J
end


function f!(resid::AbstractVector{T}, u::AbstractVector{T}) where T
    resid .= zero(T)
    resid[1] = u[1]^2 + abs2(u[2])
    resid[2] = sin(real(u[1])) + imag(u[2])
    return nothing
end



# Unified wrapper
function allocate_array(T, N)
    if _has_cuda
        return CUDA.fill(zero(T), N)
    else
        return fill(zero(T), N)
    end
end

# Driver function
function run()
    T = ComplexF64
    N = 2
    u = allocate_array(T, N)
    res = allocate_array(T, N)

    u .= 1.0 + 1.0im

    Enzyme.autodiff(Enzyme.Reverse, f!, Duplicated(res, res), Duplicated(u, u))

    println("res = ", res)
end

#run()


if _has_cuda
    println("Running on CUDA GPU")
    # use CUDA arrays
elseif _has_metal
    println("Running on Metal GPU")
    # use Metal arrays (MetalArray or MtlArray)
else
    println("Running on CPU")
    # use standard Array
end
