using Pkg

cd(@__DIR__)
Pkg.activate("../..")
using Metal


@kernel function mul2_kernel!(A)
    i = @index(Global)
    A[i]=2 * A[i]
end


@kernel function mul2_kernel2!(A)
    i = @index(Global,1)
    j = @index(Global,2)
    A[i,j]=2.0* A[i,j]
end



dev = CPU()
A=ones(1024*8,1024*8)
@time ev = mul2_kernel(dev,64)(A,ndrange=size(A))


dev =MetalBackend()
A=ones(1024*8,1024*8)

A_gpu= Adapt.adapt(MetalBackend(), Float32.(A))

#ev = mul2_kernel(dev,64)(A_gpu,ndrange=size(A_gpu))

@time ev=mul2_kernel!(dev,64)(A_gpu; ndrange=size(A_gpu))
#KernelAbstractions.wait(ev)
@print ev
KernelAbstractions.synchronize(dev)


B_gpu= Adapt.adapt(dev,ones(Float32,20,10))
@time ev=mul2_kernel2!(dev,64)(B_gpu;ndrange=([20,10]))




# -----------------------------
# Dummy dimensions
# -----------------------------
P = 4          # nPoints
L = 3          # nLs
nDim = 3       # nCoordinates
nalpha = 2     # nTotalSmallα

# -----------------------------
# Dummy GPU arrays
# -----------------------------
output_gpu       = Adapt.adapt(backendTmp, zeros(Float32, P, P, nalpha))
C_gpu            = Adapt.adapt(backendTmp,  rand(Float32, P, L, P))
int_gpu          = Adapt.adapt(backendTmp,  rand(Float32, nDim, L, L, P, P))
tableForLoop_gpu = Adapt.adapt(backendTmp,  zeros(Int32, 2+nDim*2, L*L, nalpha))
tableForPoints_gpu = Adapt.adapt(backendTmp,  zeros(Int32, nDim, P))



dummy = 1.0f0

# -----------------------------
# Kernel definition
# -----------------------------
@kernel function windowContraction2!(
    output::Float32, C::Float32
)
    xC = @index(Global,1)
    x  = @index(Global,2)
    α  = @index(Global,3)

    # Example: write some dummy data to output
    if xC <= size(output,1) && x <= size(output,2) && α <= size(output,3)
        output[xC, x, α] = 1.0 # simple test write
    end
end

# -----------------------------
# Launch kernel on Metal
# -----------------------------
ndrange=([P,P,nalpha])
windowContraction2!(dev,64)(output_gpu,C_gpu;ndrange=ndrange)

