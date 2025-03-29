using  Pkg
cd(Base.source_dir())       
Pkg.activate("../")                  # active the project, with a  static environment
# Pkg.activate(; temp=true)    #  activate the project with a temporary environment
Pkg.update()      

include("../src/imageReader.jl")

#imagefile="DSM1D/data/model/artemis/IMG_6098.jpeg"
#imagefile="DSM1D/data/model/random/tmp.png"
imagefile = "DSM1D/data/model/random/marmousi.png"
colormap = "jet" #colormap can be RGB vector or predefined colormap

floatMatrix=read2DimageModel(imagefile,colormap; showRecoveredImage=true) 