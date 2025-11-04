using  Pkg
#cd(Base.source_dir())       
cd(@__DIR__)
Pkg.activate("../..")           # active the project, with a  static environment
# Pkg.activate(; temp=true)    #  activate the project with a temporary environment
Pkg.update()      

include("../src/imageReader.jl")

#imagefile="../data/model/artemis/IMG_6098.jpeg"
#imagefile="../data/model/random/tmp.png"
imagefile = "../data/model/moi/ground_canyon.png"
colormap = "hot" #colormap can be RGB vector or predefined colormap

floatMatrix=read2DimageModel(imagefile,colormap;min=1000,max=3300, showRecoveredImage=true) 
@show size(floatMatrix)



# Brocher (2005)
function vp_from_rho(rho::Float64)
    return rho*1.5
end

function vs_from_vp(vp::Float64)
    return vp/1.7
end

# this is only for Lyon concours use!

data=floatMatrix
@show maximum(data),minimum(data)
  # Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("GCNF.rho", "w") do io
    write(io, data_vec)
end



data=vp_from_rho.(data)
@show maximum(data),minimum(data)
# Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("GCNF.vp", "w") do io
    write(io, data_vec)
end


data=vs_from_vp.(data)
@show maximum(data),minimum(data)
# Convert to Float32 (single precision)
data_single = Float32.(data)

# Flatten in column-major order (Fortran-style)
data_vec = vec(data_single)

# Write to binary file
open("GCNF.vs", "w") do io
    write(io, data_vec)
end