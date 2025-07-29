using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")
ParamFile = "../test/testparam.csv"
include("../src/DSM1D.jl")
using .DSM1D
using DIVAnd,CairoMakie
using Interpolations
using GLMakie
using Colors
using LinearAlgebra



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
    
    diam = maxX - minX
    x = range(0, diam, length=521)
    y = range(0, diam, length=521)

    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)
    colormap = myChoiceColormap(fieldname)
    hm=heatmap!(ax, x, y, fi, colormap=colormap)#, colorrange=()) if needed
    Colorbar(fig[:,2], hm)

    return fig, ax, fi
end

# below is for making a video
function myAnimation(step, iTime, fieldname, filename)
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


function lineDensityElectron1D(positionDetector, NeutrinoSource; n_pts=1000)
    x_values = range(NeutrinoSource[1], positionDetector[1] , length=n_pts)
    y_values = range(NeutrinoSource[2], positionDetector[2], length=n_pts)
    z_values = range(NeutrinoSource[3], positionDetector[3], length=n_pts)

    rad = sqrt.(x_values.^2 .+ y_values.^2 .+z_values.^2)
    radInMeter = rad.*1e3
    dens = myDensityFrom1DModel.(radInMeter)
    dist = range(0, sqrt(sum((positionDetector.-NeutrinoSource).^2)), length=n_pts)
    lines!(dist, dens)

end


function lineDensityElectron2D(positionDetector, NeutrinoSource, colorname, ax1, dR; n_pts = 100, iTime = 200)
    #draw a line between positionDetector and NeutrinoSource (coordinates) and give the density/distance profile
    fig, ax, fi = myPlot2DConvectionModel(iTime, "rho", rhoFiles)
    x_phys = range(positionDetector[1], NeutrinoSource[1], length=n_pts)
    y_phys = range(positionDetector[2], NeutrinoSource[2], length=n_pts)  
    
    lines!(ax, x_phys,y_phys, color=colorname)  # (x,y)_phys in m
    display(GLMakie.Screen(), fig)

    x_grid = x_phys ./dR
    y_grid = y_phys ./dR
    itp = interpolate(fi, BSpline(Linear()), OnGrid())
    exitp = extrapolate(itp, 0.0)

    dens = []
    for i in 1:n_pts
        x = x_grid[i]
        y = y_grid[i]
        push!(dens, exitp(x,y))
    end

    a = (positionDetector[2]-NeutrinoSource[2])/(positionDetector[1]-NeutrinoSource[1])
    if abs(a) <= tan(pi/4)
        dist = range(positionDetector[1], NeutrinoSource[1], length=n_pts)
    else
        dist = range(positionDetector[2], NeutrinoSource[2], length=n_pts)
    end

    lines!(ax1, dist, dens, color=colorname)

end


function interactiveDetector(iTime = 200)
    fig, ax, fi = myPlot2DConvectionModel(iTime, "rho", rhoFiles)
    display(fig)

    clicked_point = Observable(Point2f(NaN, NaN))

    on(events(fig.scene).mousebutton, priority = 0) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                pos = mouseposition(ax.scene)
                println("mouseposition(): $pos")
                clicked_point[] = pos
             end
        end
    end
    return clicked_point, fig, ax, fi
end


function correctedPosition(x,y; center = [6.5e6, 6.5e6])
    dx = x - center[1]
    dy = y - center[2]

    dist_radiale = sqrt(dx^2 + dy^2)
    new_x = center[1] + 6.5e6*dx/dist_radiale #modifier la valeur de 6.5e6 pour la vraie valeur du rayon?
    new_y = center[2] + 6.5e6*dy/dist_radiale
    return new_x, new_y
end


function solveQuadraticEquation(a,b,c)
    Δ = b^2 - 4*a*c

    if Δ>0
        x1 = ((-b - sqrt(Δ))/(2*a))
        x2 = ((-b + sqrt(Δ))/(2*a))
    else
        x1 = -b/(2 *a)
        x2 = x1
    end
    return x1, x2

end

function sourcePosition(center, positionDetector; earthRadius = 6.371e6, n_vectors=10)
    (xc, yc) = center[1], center[2]
    (xd, yd) = positionDetector[1], positionDetector[2]
    XY = []

    cos_θ = range(-1, 0, length = n_vectors)

    for i in 1:n_vectors
        cos_Φα1, cos_Φα2 = solveQuadraticEquation(1, -2 + 2*cos_θ[i], 1 - 2*cos_θ[i]^2)

        if 1 > cos_Φα1^2
            cos_Φα = cos_Φα1
        else
            cos_Φα =cos_Φα2
        end

        if yd > 6.5e6
            sin_Φα = sqrt(1 - cos_Φα^2)
        else
            sin_Φα = - sqrt(1 - cos_Φα^2)
        end

        α = atan((yd-yc)/(xd-xc))
        cos_Φ = cos_Φα *cos(α) - sin_Φα *sin(α)

        if yd > 6.5e6
            sin_Φ = - sqrt(1-cos_Φ^2)
        else
            sin_Φ = sqrt(1-cos_Φ^2)
        end
  
        X = xc + earthRadius*cos_Φ
        Y = yc + earthRadius*sin_Φ
        push!(XY, (X,Y))

    end
    return XY

end


function vectorsFromDetector(n_vectors = 10, scale = 2e7, diam = maxX - minX, center = [6.5e6, 6.5e6])
    #draw n_vectors (diff θ) for a positionDetector (coordinates)

    clicked_point, fig, ax, fi = interactiveDetector()
    println("Choose detector's position")

    while isnan(clicked_point[][1])
        sleep(0.1)
    end

    pos = clicked_point[]
    x, y = pos[1], pos[2]
    new_x, new_y = correctedPosition(x,y)
    XY = sourcePosition((center[1], center[2]), (new_x, new_y))

    segments_pts = []
    for source in XY
        push!(segments_pts, (new_x, new_y))
        push!(segments_pts, source)
    end

    linesegments!(ax, segments_pts)
    display(fig)

    fig1 = Figure()
    ax1 = Axis(fig1[1,1])
    
    for i in 1:n_vectors
        colorname = rand(collect(keys(Colors.color_names)))
        detector = new_x, new_y
        source = XY[i][1], XY[i][2]
        lineDensityElectron2D(detector,source, colorname, ax1, dR)
    end

    display(fig1)
end

vectorsFromDetector()


#==
test

myAnimation(20,200,"temperature", temperatureFiles)
myPlot2DConvectionModel(200, "wtr", wtrFiles)
myPlot2DConvectionModel(200, "rho", rhoFiles)
myPlot2DConvectionModel(200, "temperature", temperatureFiles)
myPlot2DConvectionModel(200, "composition", compositionFiles)
lineDensityElectron2D([450,100], [350,150])
==#