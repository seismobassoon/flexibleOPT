using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")
ParamFile = "../test/testparam.csv"
include("../src/DSM1D.jl")
using .DSM1D
using DIVAnd,CairoMakie

include("../src/batchStagYY.jl")
include("../src/parameters.jl")

function extendToCoreWithρ(ρfield, Xnode, Ynode, rcmb, dR)
    # local function here: this requires DSM1D.jl, testparam.csv
    #
    # This function will put the ρ field computed only for the mantle (and the surface) 
    # with (r,ϕ) coordinates

    # the 1D core model will be given by specifying the file to use in testparam.csv 

    ### below is just a recall how to use DSM1D module
    # PREM 
    #arrayRadius, arrayParams=DSM1D.compute1DseismicParamtersFromPolynomialCoefficients(DSM1D.my1DDSMmodel,10)

    #arrayRadius = [0.0, 30., 100., 1000., 3630., 3630., 5971., 6370., 40., 6371., 3480., 3480]

    # note that compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray is 
    # an awesome function but here I use it very naively

    arrayRadius = collect(0:dR:rcmb)
    if arrayRadius[end] != rcmb
        arrayRadius = push!(arrayRadius,rcmb)
    end

arrayRadius, arrayParams  = DSM1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(DSM1D.my1DDSMmodel, arrayRadius, "above")
DSM1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(DSM1D.my1DDSMmodel, arrayRadius.*1.e-3, "above")
    f=Figure()
    #lines(f[1,1],arrayRadius*DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer, arrayParams.ρ, markersize=1,color=:red)
    lines(f[1,1],arrayRadius, arrayParams.ρ,color=:red)
    scatter!(f[1,1],arrayRadius, arrayParams.ρ, markersize=3,color=:blue)
    display(f)
end

#dir="/Users/nobuaki/Documents/MantleConvectionTakashi/data2025"
dir="C:/Users/user/Desktop/stage 2A/données/test_estelle2/data2025"


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
#dir="/Users/nobuaki/Documents/MantleConvectionTakashi/data2025/"
dir="C:/Users/user/Desktop/stage 2A/données/test_estelle2/data2025"

compositionFiles=myListDir(dir; pattern=r"test_c\d");
temperatureFiles=myListDir(dir; pattern=r"test_t\d");
wtrFiles=myListDir(dir; pattern=r"test_wtr\d");
#dir="/Users/nobuaki/Documents/MantleConvectionTakashi/op_first_run/"
dir="C:/Users/user/Desktop/stage 2A/données/test_estelle2/op_first_run"

rhoFiles=myListDir(dir; pattern=r"test_rho\d");

# iTime definition

iTime = 100

rhoField=nothing
wtrField=nothing
file1=rhoFiles[iTime]
field1, Xnode, Ynode, rcmb = readStagYYFiles(file1)

#field1, Xnode, Ynode, rcmb = extendToCoreWithρ(field1, Xnode, Ynode, rcmb, dR)

fi1,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field1,correlationLength,epsilon2);
#rhoField = quarterDiskExtrapolation(fi1,nX,nY);



file2=wtrFiles[iTime]
field2, Xnode, Ynode, rcmb = readStagYYFiles(file2)
fi2,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field2,correlationLength,epsilon2);
#wtrField = quarterDiskExtrapolation(fi2,nX,nY);

fₚ = 0.5 .- 0.25 .* wtrField + 0.5 .* wtrField # Z/A
nₑField=rhoField .* fₚ / (massProton) 
nₑ0Field=rhoField  .* (0.5/massProton)
dnₑField = (nₑField - nₑ0Field)   ./nₑ0Field
#fig = heatmap(wtrField)
#display(fig)
fig = Figure()
ax = Axis(fig[1,1])

hm=heatmap!(ax,dnₑField,colormap=cgrad(:viridis))
Colorbar(fig[:, 2], hm)
#display(fig)


# below is for making a video

for file in rhoFiles[3:3]
    local field, Xnode, Ynode= readStagYYFiles(file)
    local fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    #local fi = quarterDiskExtrapolation(fi,nX,nY);
    local fig = Figure()
    local ax = Axis(fig[1,1],aspect = 1)
   

    local hm=heatmap!(ax,fi,colormap=cgrad(:viridis))
    Colorbar(fig[:, 2], hm)
    #display(fig)


end

#function to plot (rho(r)-rhoprem(r)/rhoprem(r))

field1, Xnode, Ynode, rcmb = readStagYYFiles(file1)
arrayRadius = sqrt.(Xnode.^2 .+ Ynode.^2)


function my_PREM(arrayRadius)
    density  = DSM1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(DSM1D.my1DDSMmodel, arrayRadius*1.e-3)
    return density
end

premDensities = my_PREM.(arrayRadius)

#enlever NaN et Infs
epsilon = 1e-6 
newpremDensities = []
for x in premDensities
    if !isfinite(x) || abs(x) < epsilon
        push!(newpremDensities, epsilon)
    else
        push!(newpremDensities,x)
    end
end

frho = (field1*1.e-3 .- newpremDensities) ./ newpremDensities
fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),frho,correlationLength,epsilon2);


#mat_frho = reshape(frho, 33540,1) #pour avoir une matrice dans heatmap!
#mat_frho_itp = interpolate(mat_frho,(BSpline(Linear()), NoInterp())) #car vecteur a 1 colonne
#mat_frho_itp_resam = Resampler(mat_frho) #trop de données donc resampler

fig2 = Figure()
ax2 = Axis(fig2[1,1],aspect = 1)
hm2=heatmap!(ax2,fi3,colormap=cgrad(:viridis))
Colorbar(fig2[:, 2], hm2)

display(fig2)