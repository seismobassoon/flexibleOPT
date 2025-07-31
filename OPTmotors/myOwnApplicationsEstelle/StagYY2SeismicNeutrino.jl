using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")
ParamFile = "../test/testparam.csv"
include("../src/DSM1D.jl")
include("../src/batchUseful.jl")
using .DSM1D
using DIVAnd,CairoMakie
using Interpolations
using GLMakie
using Colors
using LinearAlgebra


include("../src/batchStagYY.jl")
#include("testNeurthino.jl")



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
    extendToCoreWithρ!(field, Xnode, Ynode, rcmb, dR, iCheckCoreModel=false)
    quarterDiskExtrapolationRawGrid!(field, Xnode, Ynode)
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    
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
        extendToCoreWithρ!(field, Xnode, Ynode, rcmb, dR)
        quarterDiskExtrapolationRawGrid!(field, Xnode, Ynode)
        fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
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
extendToCoreWithρ!(frho, Xnode, Ynode, rcmb, dR)
quarterDiskExtrapolationRawGrid!(frho, Xnode, Ynode)
fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),frho,correlationLength,epsilon2);

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

    densGrids = []
    for i in 1:n_pts
        x = x_grid[i]
        y = y_grid[i]
        push!(densGrids, exitp(x,y))
    end

    dens=[]

    for i in 1:n_pts-1
        push!(dens, 0.5*(densGrids[i]+densGrids[i+1]))
    end

    segmentLength = sqrt((x_phys[2]-x_phys[1])^2 + (y_phys[2]-y_phys[1])) * 1.e-3 # in km
    sections = segmentLength .* ones(Float64,n_pts-1)
    
    dist = segmentLength*collect(0:1:n_pts-1) #revoir pourquoi pb de dimension

    lines!(ax1, dist, densGrids, color=colorname) #valeurs axe des abscisses à revoir
    return dens, sections
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


function correctedPosition(x,y; center = [6.5e6, 6.5e6], zposition = 2.5e3, earth_radius = 6.371e6)
    dx = x - center[1]
    dy = y - center[2]

    real_pos = earth_radius - zposition
    dist_radiale = sqrt(dx^2 + dy^2)
    new_x = center[1] + real_pos*dx/dist_radiale
    new_y = center[2] + real_pos*dy/dist_radiale
    return new_x, new_y, zposition
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

function sourcePosition(center, positionDetector; zposition=2.5e3, earthRadius = 6.371e6, n_vectors=7)
    (xc, yc) = center[1], center[2]
    (xd, yd) = positionDetector[1], positionDetector[2]
    XY = []

    cos_θ = range(-1, 0, length = n_vectors)
    θ = acos.(cos_θ)
    doubleθ = 2 .*θ
    cos_epi = cos.(doubleθ .- π)
    sin_epi = sin.(doubleθ .- π)
    rotation = cos_epi .+ im .* sin_epi


    for i in eachindex(cos_θ)
        equ = ((xd - xc) + (yd-yc)*im) * rotation[i]
        X = real(equ)+xc
        Y = imag(equ)+yc
        # this is the position with zposition below

        newX = nothing
        newY = nothing
        if cos_θ[i] !== 0.0
            if X-xd != 0.0
                slope = (Y-yd)/(X-xd)

                a = 1+slope^2
                b = -2*xc + 2*Y*slope-2*slope^2*X-2*slope*yc
                c = xc^2 + Y^2 - 2*Y*slope*X + slope^2*X^2 -2*Y*yc + 2*slope*X*yc +yc^2 - earthRadius^2
                sol1,sol2 = solveQuadraticEquation(a,b,c)
                if (sol1-xd)*(X-xd)>0.0
                    newX = sol1
                    newY = Y + slope* (newX - X)
                else
                    newX = sol2
                    newY = Y + slope* (newX - X)
                end

            else
                slope = (X-xd)/(Y-yd)
                
                a = 1+slope^2
                b = -2*yc + 2*X*slope-2*slope^2*Y-2*slope*xc
                c = yc^2 + X^2 - 2*X*slope*Y + slope^2*Y^2 -2*X*xc + 2*slope*Y*xc +xc^2 - earthRadius^2
                sol1,sol2 = solveQuadraticEquation(a,b,c)
                if (sol1-yd)*(Y-yd)>0.0
                    newY = sol1
                    newX = Y + slope* (newY - Y)
                else
                    newY = sol2
                    newX = Y + slope* (newY - Y)
                end

            end

        else
            segmentfromDtoS = sqrt(earthRadius^2-(earthRadius-zposition)^2)
            newX = xd + -(yd-yc)/(earthRadius-zposition)*segmentfromDtoS
            newY = yd + (xd-xc)/(earthRadius-zposition)*segmentfromDtoS
        end
        push!(XY, (newX,newY))

    end
    return XY

end

#== à continuer
function posOrNeg(cos_θ, sign = :positive)
    if sign== :positive
        θ = acos.(cos_θ)
    else
        θ = .- acos.(cos_θ)
    end
    @show θ
end

#theta = posOrNeg(range(-1, 0, 4), :negative)
==#

function vectorsFromDetector(n_vectors = 7, center = [6.5e6, 6.5e6])
    #draw n_vectors (diff θ) for a positionDetector (coordinates)

    clicked_point, fig, ax, fi = interactiveDetector()
    println("Choose detector's position")

    while isnan(clicked_point[][1])
        sleep(0.1)
    end

    pos = clicked_point[]
    x, y = pos[1], pos[2]
    new_x, new_y,zposition = correctedPosition(x,y) 
    XY = sourcePosition((center[1], center[2]), (new_x, new_y); zposition=zposition)

    segments_pts = []
    for source in XY
        push!(segments_pts, (new_x, new_y))
        push!(segments_pts, (source[1], source[2]))
    end

    linesegments!(ax, segments_pts)
    display(fig)

    fig1 = Figure()
    ax1 = Axis(fig1[1,1])

    densities_list = []
    sections_list = []
    for i in 1:n_vectors
        colorname = rand(collect(keys(Colors.color_names)))
        @show detector = new_x, new_y
        @show source = XY[i][1], XY[i][2]
        dens, section = lineDensityElectron2D(detector,source, colorname, ax1, dR)
        push!(densities_list, dens)
        push!(sections_list, section)

    end

    #display(fig1)
    return densities_list, sections_list
end

densities ,sections=vectorsFromDetector()
export densities, sections







#==
function linkWithNeurthino()
    densities_list, sections_list = vectorsFromDetector()
    for i in 1:7
        paths = Path(densities_list[i], sections_list[i])
        energies = 10 .^ range(0, stop=2, length=7);
        prob = Pνν(U, H, energies, paths);
    end
    return energies, prob
end

#energies, prob = linkWithNeurthino()
#@show energies, prob
==#



#==
test

myAnimation(20,200,"temperature", temperatureFiles)
myPlot2DConvectionModel(200, "wtr", wtrFiles)
myPlot2DConvectionModel(200, "rho", rhoFiles)
myPlot2DConvectionModel(200, "temperature", temperatureFiles)
myPlot2DConvectionModel(200, "composition", compositionFiles)
lineDensityElectron2D([450,100], [350,150])
==#