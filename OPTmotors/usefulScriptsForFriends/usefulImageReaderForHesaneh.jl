using  Pkg
#cd(Base.source_dir())       
cd(@__DIR__)
Pkg.activate("../..")           # active the project, with a  static environment
# Pkg.activate(; temp=true)    #  activate the project with a temporary environment
Pkg.update()      

include("../src/imageReader.jl")

#imagefile="DSM1D/data/model/artemis/IMG_6098.jpeg"
#imagefile="DSM1D/data/model/random/tmp.png"
imagefile = "../data/model/moi/ground_canyon.png"
colormap = "hot" #colormap can be RGB vector or predefined colormap

floatMatrix=read2DimageModel(imagefile,colormap;min=1.0,max=3.3, showRecoveredImage=true) 


# Brocher (2005)
function vp_from_rho(rho::Float64)
    return -24.629 + 116.00*rho - 182.34*rho^2 + 121.90*rho^3 - 25.478*rho^4
end

function vs_from_vp(vp::Float64)
    return -1.361 + 1.483*vp - 0.822*vp^2 + 0.260*vp^3 - 0.032*vp^4
end

# this is only for Lyon concours use!

data=floatMatrix
  # Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("GCNF.rho", "w") do io
    write(io, data_vec)
end



data=vp_from_rho.(data)

# Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("GCNF.vp", "w") do io
    write(io, data_vec)
end


data=vs_from_vp.(data)
# Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("GCNF.vs", "w") do io
    write(io, data_vec)
end