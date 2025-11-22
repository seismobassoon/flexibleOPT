using  Pkg
#cd(Base.source_dir())       
cd(@__DIR__)
Pkg.activate("../..")           # active the project, with a  static environment
# Pkg.activate(; temp=true)    #  activate the project with a temporary environment
Pkg.update()      
Pkg.update("Makie")
Pkg.update("Polynomials")

include("../src/imageReader_hesa.jl")

#imagefile="../data/model/artemis/IMG_6098.jpeg"
#imagefile="../data/model/random/tmp.png"
#imagefile = "../data/model/moi/ground_canyon.png"
imagefile = "/Users/hessiemohammadi/Documents/FUJI/events/rock_no_topo_air.jpg"
#colormap = "hot" #colormap can be RGB vector or predefined colormap

# Define geology color - velocity mapping
geologyColors = [
    
    RGB(1.0, 1.0, 1.0),     # Air (white)
    RGB(1.0, 0.0, 0.0),    # Magma (red)
    RGB(1.0, 0.65, 0.0),   # Mush (orange)
    RGB(0.55, 0.27, 0.07)  # Rock (brown)
]

targetVp = [
   
    300.0,    # air : very low Vp
    2000.0,  # magma : lowest Vp
    4000.0,  # mush  : intermediate Vp
    6000.0   # rock  : highest Vp   
]

rhoValues  = targetVp ./ 1.5    # density values corresponding to Vp

# Read categorical model  
floatMatrix = read2DimageModel(imagefile;
    colorbar = geologyColors,
    values   = rhoValues,
    showRecoveredImage = true
)
#floatMatrix=read2DimageModel(imagefile,colormap;min=1000,max=3300, showRecoveredImage=true) 
#@show size(floatMatrix)

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

vp_from_rho(rho) = rho * 1.5
vs_from_vp(vp)   = vp / 1.7

# this is only for Lyon concours use!

data=floatMatrix
@show maximum(data),minimum(data)
  # Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("rock_no_topo_air2.rho", "w") do io
    write(io, data_vec)
end



data=vp_from_rho.(data)
@show maximum(data),minimum(data)
# Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("rock_no_topo_air2.vp", "w") do io
    write(io, data_vec)
end


data=vs_from_vp.(data)
@show maximum(data),minimum(data)
# Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("rock_no_topo_air2.vs", "w") do io
    write(io, data_vec)
end