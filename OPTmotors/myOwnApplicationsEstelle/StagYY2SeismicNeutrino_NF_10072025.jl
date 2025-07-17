using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")
ParamFile = "../test/testparam.csv"
include("../src/DSM1D.jl")
using .DSM1D
using DIVAnd,CairoMakie

include("../src/batchStagYY.jl")

include("../src/parameters.jl")


function myDensityFrom1DModel(arrayRadius)
    radiusInKilometer = arrayRadius*1.e-3
    density  = DSM1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(DSM1D.my1DDSMmodel, radiusInKilometer)

    return density
end

function myChoiceColormap(fieldname)
    if fieldname ==="temperature"
        colormap=cgrad(:seismic)

    elseif fieldname === "rho" || fieldname === "composition"
        colormap=cgrad(:viridis)

    elseif fieldname === "wtr"
        colormap=cgrad(:blues)

    else 
        colormap=cgrad(:viridis)

    end
end 

function myPlot2DConvectionModel(iTime, fieldname, filename)
#only if the field in DIVandrun is the same as in readStagYYFiles
    file = filename[iTime]
    field, Xnode, Ynode, rcmb = readStagYYFiles(file)
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    fi = quarterDiskExtrapolation(fi,nX,nY)
    
    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)

    colormap = myChoiceColormap(fieldname)
    hm=heatmap!(ax, fi, colormap=colormap)#, colorrange=()) if needed
    Colorbar(fig[:,2], hm)
    display(fig)
end

# below is for making a video
function myAnimation(step,iTime,fieldname,filename)
    fi_list=[]
    for file in filename[1:step:iTime]
        field, Xnode, Ynode= readStagYYFiles(file)
        fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
        fi = quarterDiskExtrapolation(fi,nX,nY)
        push!(fi_list, fi)
    end

    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    data = Observable(fi_list[1])
    colormap = myChoiceColormap(fieldname)
    hm=heatmap!(ax, data, colormap=colormap)#, colorrange=()) if needed
    Colorbar(fig[:, 2], hm)
    record(fig, "animation2D.mp4", 1:length(fi_list); framerate = 2) do i
        data[] = fi_list[i]
    end
end

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
dir="C:/Users/user/Desktop/stage 2A/données/MantleConvectionTakashi/data2025"


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

# file types
dir="C:/Users/user/Desktop/stage 2A/données/MantleConvectionTakashi/data2025/"
dir="C:/Users/user/Desktop/stage 2A/données/MantleConvectionTakashi/op_first_run/"
rhoFiles=myListDir(dir; pattern=r"test_rho\d");
compositionFiles=myListDir(dir; pattern=r"test_c\d");
temperatureFiles=myListDir(dir; pattern=r"test_t\d");
wtrFiles=myListDir(dir; pattern=r"test_wtr\d");


#plot (rho(r)-rhoprem(r)/rhoprem(r))
iTime = 150
file1=rhoFiles[iTime]

field1, Xnode, Ynode, rcmb = readStagYYFiles(file1)
densitiesInGcm3 = field1*1e-3
arrayRadius = sqrt.(Xnode.^2 .+ Ynode.^2)

premDensities = myDensityFrom1DModel.(arrayRadius)
newpremDensities = premDensities
frho = ifelse.(newpremDensities .== 0.0, 0.0, (densitiesInGcm3 .- newpremDensities) ./ newpremDensities)
fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),frho,correlationLength,epsilon2);
fi3 = quarterDiskExtrapolation(fi3,nX,nY)


fig2 = Figure()
ax2 = Axis(fig2[1,1],aspect = 1)
hm2 = heatmap!(ax2,fi3,colormap=cgrad(:viridis),colorrange=(-0.005,0.005))

Colorbar(fig2[:, 2], hm2)

#display(fig2)


function lineDensityElectron1D(positionDetector, NeutrinoSource)
    n_pts = 1000
    x_values = range(NeutrinoSource[1], positionDetector[1] , length=n_pts)
    y_values = range(NeutrinoSource[2], positionDetector[2], length=n_pts)
    z_values = range(NeutrinoSource[3], positionDetector[3], length=n_pts)

    rad = sqrt.(x_values.^2 .+ y_values.^2 .+z_values.^2)
    radInMeter = rad.*1e3
    dens = myDensityFrom1DModel.(radInMeter)
    dist = range(0, sqrt(sum((positionDetector.-NeutrinoSource).^2)), length=n_pts)
    lines!(dist, dens)

end


using Interpolations

function lineDensityElectron2D(positionDetector, NeutrinoSource)
    #draw a line between positionDetector and NeutrinoSource (coordinates) and give the density/distance profile

    file = rhoFiles[iTime]
    field, Xnode, Ynode, rcmb = readStagYYFiles(file)
    densitiesInGcm3 = field*1e-3
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),densitiesInGcm3,correlationLength,epsilon2);
    fi = quarterDiskExtrapolation(fi,nX,nY)
    
    n_pts = 100
    x_values = range(positionDetector[1], NeutrinoSource[1], length=n_pts)
    y_values = range(positionDetector[2], NeutrinoSource[2], length=n_pts)

    itp = interpolate(fi, BSpline(Linear()), OnGrid())
    dens = []
    for i in 1:n_pts
        x = x_values[i]
        y = y_values[i]
        push!(dens, itp(x,y))
    end

    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    hm = heatmap!(ax,fi,colormap=cgrad(:viridis)) 
    lines!(x_values,y_values) 
    Colorbar(fig[:, 2], hm)
    display(fig)


    #dist = range(0, sqrt(sum((positionDetector.-NeutrinoSource).^2)), length=n_pts)
    a = (positionDetector[2]-NeutrinoSource[2])/(positionDetector[1]-NeutrinoSource[1])
    @show a
    if abs(a) <= tan(pi/4)
        dist = range(positionDetector[1], NeutrinoSource[1], length=n_pts)
    else
        dist = range(positionDetector[2], NeutrinoSource[2], length=n_pts)
    end

    fig1 = Figure()
    ax1 = Axis(fig1[1,1], aspect = 1)
    lines!(ax1, dist, dens) #axe dist à modifier
    display(fig1)

end
#==
lineDensityElectron2D([450,100], [350,150])
lineDensityElectron2D([450,100], [410,250])
==#

function detectorcosθ(positionDetector)
    #draw n_vectors (diff θ) for a positionDetector (coordinates)

    file = rhoFiles[200]
    field, Xnode, Ynode, rcmb = readStagYYFiles(file)
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    fi = quarterDiskExtrapolation(fi,nX,nY)
    
    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)

    colormap = myChoiceColormap("rho")
    hm=heatmap!(ax, fi, colormap=colormap)#, colorrange=()) if needed
    Colorbar(fig[:,2], hm)

    n_vectors = 10
    θ_values = range(0.0, 2*pi, length=n_vectors)
    x0 = fill(positionDetector[1], n_vectors)
    y0 = fill(positionDetector[2], n_vectors)

    #pour agrandir les vecteurs
    scale = 500
    dx = scale.*cos.(θ_values)
    dy = scale.*sin.(θ_values)

    #si on depasse de la grille
    X = dx .+ x0
    Y = dy .+ y0
    right = X.>521
    dx[right] .= 521 .-x0[right] 
    left = X.<0
    dx[left] .= 0 .- x0[left]
    above = Y.>521
    dy[above] .= 521 .- y0[above]
    below = Y.<0
    dy[below] .= 0 .- y0[below]

    quiver!(ax, x0, y0, dx, dy, tipwidth=8, tiplength=15, linewidth=1.5)
    display(fig)

end

detectorcosθ([80, 80])
detectorcosθ([80, 450])
detectorcosθ([450, 450])
detectorcosθ([450, 80])


#==
test

myAnimation(20,200,"temperature", temperatureFiles)
myPlot2DConvectionModel(200, "wtr", wtrFiles)
myPlot2DConvectionModel(200, "rho", rhoFiles)
myPlot2DConvectionModel(200, "temperature", temperatureFiles)
myPlot2DConvectionModel(200, "composition", compositionFiles)
==#
