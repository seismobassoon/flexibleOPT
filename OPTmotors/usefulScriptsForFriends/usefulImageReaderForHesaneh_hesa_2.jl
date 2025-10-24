using  Pkg
#cd(Base.source_dir())       
cd(@__DIR__)
Pkg.activate("../..")           # active the project, with a  static environment
# Pkg.activate(; temp=true)    #  activate the project with a temporary environment
Pkg.update()      
Pkg.update("Makie")


include("../src/imageReader_hesa.jl")

#imagefile="../data/model/artemis/IMG_6098.jpeg"
#imagefile="../data/model/random/tmp.png"
#imagefile = "../data/model/moi/ground_canyon.png"
imagefile = "/Users/hessiemohammadi/Documents/FUJI/events/rock_potato_topo_air.jpg"
#colormap = "hot" #colormap can be RGB vector or predefined colormap

# Define geology color - velocity mapping
geologyColors = [
    
    RGB(1.0, 1.0, 1.0),     # Air (white)
    RGB(1.0, 0.0, 0.0),    # Magma (red)
    RGB(1.0, 0.65, 0.0),   # Mush (orange)
    RGB(0.55, 0.27, 0.07)  # Rock (brown)
]

targetVp = [
   
    343.0,    # air : very low Vp
    2000.0,  # magma : lowest Vp
    4000.0,  # mush  : intermediate Vp
    6000.0   # rock  : highest Vp   
]

#rhoValues  = targetVp ./ 1.5    # density values corresponding to Vp

# Adjust rhoValues to avoid zero density for air
#rhoValues = [ vp==343.0 ? 100 : 1000 + 0.32 * vp for vp in targetVp ]
rhoValues = [ vp==343.0 ? 100 : vp ./ 1.5 for vp in targetVp ]
 # small nonzero rho for air

# Read categorical model  
floatMatrix = read2DimageModel(imagefile;
    colorbar = geologyColors,
    values   = rhoValues,
    showRecoveredImage = true
)
#floatMatrix=read2DimageModel(imagefile,colormap;min=1000,max=3300, showRecoveredImage=true) 
@show size(floatMatrix)

#size(floatMatrix) = (nz, nx)


# Brocher (2005)
"""
function vp_from_rho(rho::Float64)
    return rho*1.5
end

function vs_from_vp(vp::Float64)
    return vp/1.7
end
"""

#vp_from_rho(rho) = (rho - 1000) / 0.32
vp_from_rho(rho) = rho * 1.5
vs_from_vp(vp)   = vp / 1.7

# this is only for Lyon concours use!

# Density
data=floatMatrix
@show maximum(data),minimum(data)
  # Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("rock_potato_topo_air.rho", "w") do io
    write(io, data_vec)
end


# Velocity
data_vp = vp_from_rho.(data)
data_vp[abs.(data .- 100) .< 1e-3] .= 343.0
 
#data_vp[floatMatrix .== rhoValues[1]] .= 343.0   # restore 0 Vp where air
@show maximum(data_vp), minimum(data_vp)

open("rock_potato_topo_air.vp", "w") do io
    write(io, Float32.(vec(data_vp)))
end


# Shear Velocity
data_vs = vs_from_vp.(data_vp)
data_vs[data_vp .== 343.0] .= 0               # ensure Vs=0 where air
@show maximum(data_vs), minimum(data_vs)

open("rock_potato_topo_air.vs", "w") do io
    write(io, Float32.(vec(data_vs)))
end