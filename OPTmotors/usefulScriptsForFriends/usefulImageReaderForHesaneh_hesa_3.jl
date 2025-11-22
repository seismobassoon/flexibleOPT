using Pkg
using Printf        

cd(@__DIR__)
Pkg.activate("../..")
Pkg.update()
Pkg.update("Makie")
include("../src/batchRevise.jl")

include("../src/imageReader_hesa.jl")


imagefile = "/Users/hessiemohammadi/Documents/FUJI/events/rock_chamber_new.jpg"

#  Define geology color - velocity mapping 
geologyColors = [
    RGB(1.0, 1.0, 1.0),    # Air (white)
    RGB(1.0, 0.0, 0.0),    # Magma (red)
    RGB(1.0, 0.65, 0.0),   # Mush (orange)
    RGB(0.55, 0.27, 0.07)  # Rock (brown)
]

targetVp = [
    1500.0,   # Air (m/s)
    2000.0,  # Magma
    4000.0,  # Mush
    6000.0   # Rock
]


# Brocher (2005, BSSA 95:2081–2092)
# Valid for 1.5–8.5 km/s; here Vp in m/s, ρ in kg/m³
function rho_from_vp_brocher(vp)
    vp64 = Float64(vp)  
    if vp64 ≤ 343.0
        return 100.0
    else
        vp_km = vp64 / 1000.0
        ρ_gcm3 = 1.6612*vp_km - 0.4721*vp_km^2 + 0.0671*vp_km^3 -
                  0.0043*vp_km^4 + 0.000106*vp_km^5
        return ρ_gcm3 * 1000.0   # kg/m³
    end
end

# Shear-wave velocity from Vp (approximate crustal ratio)
vs_from_vp(vp) = vp ≤ 343.0 ? 0.0 : vp / 1.7

# Read categorical model (Vp-based) 
floatMatrix_vp = read2DimageModel(
    imagefile;
    colorbar = geologyColors,
    values   = targetVp,
    showRecoveredImage = true
)

println("\n=== Model loaded ===")
println("Matrix size: ", size(floatMatrix_vp))
println(@sprintf("Vp range: min = %.2f m/s, max = %.2f m/s",
    minimum(floatMatrix_vp), maximum(floatMatrix_vp)))

#  Compute ρ and Vs from Vp 
rhoMatrix = rho_from_vp_brocher.(floatMatrix_vp)
vsMatrix  = vs_from_vp.(floatMatrix_vp)


# Save to binary files 
open("rock_chamber_new.vp", "w") do io
    write(io, Float32.(vec(floatMatrix_vp)))
end

open("rock_chamber_new.rho", "w") do io
    write(io, Float32.(vec(rhoMatrix)))
end

open("rock_chamber_new.vs", "w") do io
    write(io, Float32.(vec(vsMatrix)))
end


println(@sprintf("Vp range : %.2f – %.2f m/s", minimum(floatMatrix_vp), maximum(floatMatrix_vp)))
println(@sprintf("ρ range  : %.2f – %.2f kg/m³", minimum(rhoMatrix), maximum(rhoMatrix)))
println(@sprintf("Vs range : %.2f – %.2f m/s", minimum(vsMatrix), maximum(vsMatrix)))
