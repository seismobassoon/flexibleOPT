using Enzyme
using Requires

function __init__()
    if Sys.isapple()
        @require Metal="dde4c033-4e86-420c-a63e-0dd931031962" begin
            println("Metal backend available")
        end
    else
        @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" begin
            println("CUDA backend available")
        end
    end
end


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
