using WriteVTK
#using FortranFiles
using Printf
using LinearAlgebra


function write_vtk_structured(filename::String, x, y, z, data_dict::Dict{String,Array{Float32,3}})
    vtk_grid(filename, x, y, z) do vtk
        for (name, field) in data_dict
            vtk[name] = field
        end
    end
end

# Fonction pour lire un enregistrement Fortran séquentiel
function read_fortran_record(f, count::Int, dtype::Type)
    rec_start = read(f, Int32)

    # Promote dtype to a preferred Julia-native type
    T = dtype === Int32 ? Int :
        dtype === Float32 ? Float64 :
        dtype  # keep original type if no match

    data = [convert(T, read(f, dtype)) for _ in 1:count]

    rec_end = read(f, Int32)
    @assert rec_start == rec_end "Fortran record length mismatch"

    return data
end

function binRead(io,dataTypeTmp::DataType; byte_reverse_in=false)

    value=read(io,dataTypeTmp)
    if byte_reverse_in
        value=bswap(value)
    end
        
    return value
end

# Fonction pour lire binaire direct (fichier ouvert normalement)
function read_direct_binary(f::IO, count::Int, dtype::Type)
    data = [binRead(f,dtype) for _ in 1:count]
    return data
end

# Chemin vers le fichier grid
fname_grid = "eq.Explosion.Z.Z.alphaV.grid"

# Ouvre et lit le fichier grid une fois (informations constantes)
#fgrid = FortranFile(fname_grid, "r")
fgrid = open(fname_grid,"r")
dims = read_fortran_record(fgrid, 4, Int32)
@show (nr, nphi, ntheta, nktype) = dims
nktype += 1  # Ajustement selon votre remarque


@show radii  = read_fortran_record(fgrid, nr, Float32)
phis   = read_fortran_record(fgrid, nphi * ntheta, Float32)
thetas = read_fortran_record(fgrid, nphi * ntheta, Float32)
close(fgrid)




# Reshape des angles en 2D (ntheta × nphi)
#phis = reshape(phis, ntheta, nphi)
#thetas = reshape(thetas, ntheta, nphi)

phis = reshape(phis, nphi, ntheta)
thetas = reshape(thetas, nphi, ntheta)

@info "Lecture fichier grid terminée."
println("Dimensions: nr=$nr, nphi=$nphi, ntheta=$ntheta")
println("Types kernel: nktype=$nktype")
println("Rayons: $(radii[1]) → $(radii[end])")
println("Phis: $(phis[1]) → $(phis[end])")
println("Thetas: $(thetas[1]) → $(thetas[end])")

# Boucle sur les snapshots
for itime in 1:135
    num_snap = @sprintf("%07d", itime)
    fname_kernel = "tmpvideo/eq.Explosion.Z.Z.100s2s.$num_snap.video"
    fname_vtk = "newVTK/ker$num_snap.vts"

    # Ouvre fichier kernel en binaire direct
    fker = open(fname_kernel, "r")
    npoints = nr * nphi * ntheta
    kernel_raw = read_direct_binary(fker, npoints * nktype, Float32)
    close(fker)

    # Reshape kernel : (nktype, ntheta, nphi, nr)
    #kernel = reshape(kernel_raw, (nktype, ntheta, nphi, nr))
    kernel = reshape(kernel_raw, (nr, nphi, ntheta, nktype))
  
    # Calcul des coordonnées cartésiennes
    sinθ = sin.(deg2rad.(thetas))
    cosθ = cos.(deg2rad.(thetas))
    cosϕ = cos.(deg2rad.(phis))
    sinϕ = sin.(deg2rad.(phis))

    #===
    xgrid = [sinθ[i,j] * cosϕ[i,j] * r for i in 1:ntheta, j in 1:nphi, r in radii]
    ygrid = [sinθ[i,j] * sinϕ[i,j] * r for i in 1:ntheta, j in 1:nphi, r in radii]
    zgrid = [cosθ[i,j] * r for i in 1:ntheta, j in 1:nphi, r in radii]

    xgrid = reshape(xgrid, (ntheta, nphi, nr))
    ygrid = reshape(ygrid, (ntheta, nphi, nr))
    zgrid = reshape(zgrid, (ntheta, nphi, nr))
    ===#


    xgrid = [sinθ[j,i] * cosϕ[j,i] * r for r in radii, j in 1:nphi, i in 1:ntheta]
    ygrid = [sinθ[j,i] * sinϕ[j,i] * r for r in radii, j in 1:nphi, i in 1:ntheta]
    zgrid = [cosθ[j,i] * r for r in radii, j in 1:nphi, i in 1:ntheta]

    xgrid = reshape(xgrid, (nr,nphi,ntheta))
    ygrid = reshape(ygrid, (nr,nphi,ntheta))
    zgrid = reshape(zgrid, (nr,nphi,ntheta))


    # Crée un fichier VTK (structured grid)
    #vtk = vtk_grid(fname_vtk, xgrid, ygrid, zgrid)
    vtk=Dict()
    # Ajout des données du kernel comme point data
    for i in 1:nktype
        vtk["kernel_$(i-1)"] = kernel[:, :, :,i]    
        @show maximum(kernel[:,:,:,i])
        
    end

    vtk = Dict{String, Array{Float32,3}}(vtk)
    xgrid = Float32.(xgrid)
    ygrid = Float32.(ygrid)
    zgrid = Float32.(zgrid)
    write_vtk_structured(fname_vtk, xgrid, ygrid, zgrid,vtk)
    

    @info "Snapshot $num_snap écrit dans $fname_vtk"
end
