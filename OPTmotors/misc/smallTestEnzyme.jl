using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")

include("../src/imageReader.jl") # read 2D images for models

include("../src/OPTwrappers.jl") 

using Enzyme
using SparseArrays
using Enzyme

_has_cuda = try
    @eval using CUDA
    CUDA.has_cuda()  # ✅ call the function, don't shadow or reassign
catch
    false
end


function enzyme_column_jacobian!(J::SparseMatrixCSC, f!, u::CuArray, res::CuArray)
    n = length(u)
    du = similar(u)
    df = similar(res)

    for j in 1:n
        CUDA.fill!(du, zero(eltype(u)))
        du[j] = one(eltype(u))
        CUDA.fill!(df, zero(eltype(res)))

        Enzyme.autodiff(Enzyme.Reverse, f!, Duplicated(df, res), Duplicated(du, u))

        # Now df contains the j-th column of J
        for i in 1:length(res)
            if df[i] != 0
                J[i,j] = df[i]
            end
        end
    end
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

run()


