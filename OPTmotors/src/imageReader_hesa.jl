include("../src/batchImages.jl")

using FileIO,CairoMakie, LinearAlgebra

# Convert each pixel color to a float value based on the nearest color in the colorbar
function color2Float(image, colorbar, values)
    nz, nx = size(image)[1:2]
    float_array = zeros(Float32, nz, nx)
    for i in 1:nz, j in 1:nx
        c = RGB(image[i,j])
        # find nearest reference color in RGB space
        dists = [norm([red(c)-red(r), green(c)-green(r), blue(c)-blue(r)]) for r in colorbar]
        k = argmin(dists)
        float_array[i,j] = values[k]
    end
    return float_array
end


function read2DimageModel(file,cmap;min=0.0,max=1.0,showRecoveredImage=true,Nheight=nothing,Nwidth=nothing)
    
    if cmap isa String
        cmap=getColorPalette(cmap)
    end
    
    rgbpoints=length(cmap)
    values=range(min,max,rgbpoints)
    float_array=read2DimageModel(file; colorbar=cmap,values=values,showRecoveredImage=showRecoveredImage,Nheight,Nwidth)
    return float_array
end


function read2DimageModel(file; Ncolor=256, colorbar = [RGB(0.0, 0.0, 1.0), RGB(0.0, 1.0, 0.0), RGB(1.0, 0.0, 0.0)] ,values = [0.0, 0.5, 1.0], lowestColor=nothing,highestColor=nothing,Nheight=nothing,Nwidth=nothing,showRecoveredImage=true)

    # Read RGB and take the absolute value

    
    image = load(file)

    #rescaling the image
    if Nheight !== nothing || Nwidth !== nothing
        height, width = size(image)[1:2]
        if Nheight === nothing
            Nheight=trunc.(Int,(height-1)*(Nwidth-1)/(width-1)+1)
        end
        if Nwidth === nothing
            Nwidth=trunc.(Int,(width-1)*(Nheight-1)/(height-1)+1)
        end
        image = rescalingImage(image,Nheight,Nwidth)
    end

    if lowestColor !== nothing && highestColor !== nothing
        averageColor=sum(image)/length(image)
        mediumValue = 0.5*(values[1]+values[end])
        values = [values[1],mediumValue,values[end]]
        colorbar=[lowestColor,averageColor,highestColor]
    end

    # Convert the image to float using nearest-color mapping
    float_array=color2Float(image,colorbar,values)

    #@show values[1],values[end]
    #@show minimum(float_array), maximum(float_array)
    # Define colorbar and corresponding float values
    

    if showRecoveredImage
   
        newColorBar,newValues=regenerataionColorMap(colorbar,values,Ncolor)
        scene = heatmap(transpose(float_array[end:-1:1,:]), colormap =  newColorBar,colorrange=(values[1],values[end]))
        display(scene)
    end
   
    return float_array
    
end


greenColorMap=[RGB(1.0,1.0,1.0),RGB(0.0,1.0,0.0)]

#file="DSM1D/data/model/random/colourful.jpg"
#file="DSM1D/data/model/artemis/IMG_6098.jpeg"
#file="DSM1D/data/model/random/marmousi.png"
#file="DSM1D/data/model/random/tmp.png"

#read2DimageModel(file,"jet";Nwidth=101,Nheight=202)

#read2DimageModel(file,"RdYlGn")
#read2DimageModel(file,greenColorMap)

# RdYlGn is not bad for Artemis

#read2DimageModel(file;lowestColor=RGB(0.0,1.0,0.0),highestColor=RGB(0.0,0.0,1.0))
#read2DimageModel(file;lowestColor=RGB(0.0,1.0,0.0),highestColor=RGB(1.0,0.0,0.0))
