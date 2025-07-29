cd(@__DIR__)
Pkg.activate("../..")

ParamFile = "../test/testparam.csv"
include("../src/DSM1D.jl")
using .DSM1D
using DIVAnd,CairoMakie

include("../src/batchStagYY.jl")
include("../src/parameters.jl")
function myDensityFrom1DModel(arrayRadius)

    density  = DSM1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(DSM1D.my1DDSMmodel, arrayRadius*1.e-3)

    return density
end



#dir="/Users/nobuaki/Documents/MantleConvectionTakashi/data2025"
dir="/Users/nobuaki/Documents/MantleConvectionTakashi/op_first_run"


boolFlat = true # we can read but for me it is better to have this information before reading since DIVAnd_rectdom can be applied before reading

# Cartesian grids and interpolation

minX,maxX,nX = -6500e3, 6500e3, 521
minY,maxY,nY = minX,maxX,nX
minZ,maxZ,nZ = minX,maxX,nX
dR = (maxX-minX)/(nX-1) # the interval in X, which we suppose to be the smallest grid interval


correlationLength=(20e3,20e3,20e2) # not yet fully understood this for DIV interpolation
epsilon2 =1.;


if boolFlat
    nZ=1
    minZ=0.0
    maxZ=0.0
    tmpX=correlationLength[1]
    tmpY=correlationLength[2]
    correlationLength=(tmpX,tmpY)

    mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),
                                            range(minY,stop=maxY,length=nY));
else

    mask,(pm,pn,po),(xi,yi,zi) = DIVAnd_rectdom(range(minX,stop=maxX,length=nX),
                                            range(minY,stop=maxY,length=nY),
                                            range(minZ,stop=maxZ,length=nZ));
end


# taking core info from PREM (or something else)



# file types
dir="/Users/nobuaki/Documents/MantleConvectionTakashi/data2025/"
compositionFiles=myListDir(dir; pattern=r"test_c\d");
temperatureFiles=myListDir(dir; pattern=r"test_t\d");
wtrFiles=myListDir(dir; pattern=r"test_wtr\d");
dir="/Users/nobuaki/Documents/MantleConvectionTakashi/op_first_run/"
rhoFiles=myListDir(dir; pattern=r"test_rho\d");



#@show length(rhoFiles)

# iTime definition

iTime = 100

rhoField=nothing
wtrField=nothing



file1=rhoFiles[iTime]
field1, Xnode, Ynode, rcmb = readStagYYFiles(file1)
extendToCoreWithρ!(field1, Xnode, Ynode, rcmb, dR)
#quarterDiskExtrapolationRawGrid!(fi, Xnode, Ynode)
fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field1,correlationLength,epsilon2);
fi = quarterDiskExtrapolation(fi,nX,nY);
fig = Figure()
ax = Axis(fig[1,1],aspect = 1)
hm=heatmap!(ax,fi,colormap=cgrad(:viridis))
Colorbar(fig[:, 2], hm)
display(fig)



fi1,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field1,correlationLength,epsilon2);
rhoField = quarterDiskExtrapolation(fi1,nX,nY);



@show file2=wtrFiles[iTime]
field2, Xnode, Ynode, rcmb = readStagYYFiles(file2)
fi2,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field2,correlationLength,epsilon2);
wtrField = quarterDiskExtrapolation(fi2,nX,nY);


# this is NOT a good equation to estimate electron density!!

fₚ = 0.5 .- 0.25 .* wtrField + 0.5 .* wtrField # Z/A



nₑField=rhoField .* fₚ / (massProton) 
nₑ0Field=rhoField  .* (0.5/massProton)
dnₑField = (nₑField - nₑ0Field)   ./nₑ0Field

#fig = heatmap(wtrField)
#display(fig)
fig = Figure()
ax = Axis(fig[1,1])#

hm=heatmap!(ax,dnₑField,colormap=cgrad(:viridis))
Colorbar(fig[:, 2], hm)
display(fig)
#@show length(wtrFiles)

# below is for making a video

for file in compositionFiles[100:100]
    local field, Xnode, Ynode= readStagYYFiles(file)
    local fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    local fi = quarterDiskExtrapolation(fi,nX,nY);
    local fig = Figure()
    local ax = Axis(fig[1,1],aspect = 1)
   

    local hm=heatmap!(ax,fi,colormap=cgrad(:viridis))
    Colorbar(fig[:, 2], hm)
    display(fig)

end





