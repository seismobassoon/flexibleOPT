using  Pkg, BenchmarkTools

cd(@__DIR__)
Pkg.activate("../..")
ParamFile = "../test/testparam.csv"
include("../src/DSM1D.jl")
include("../src/batchUseful.jl")
using .DSM1D
using DIVAnd,CairoMakie
using Interpolations
import GLMakie
using Colors
include("../src/batchStagYY.jl")
include("../src_Neurthino/Neurthino.jl")
using .Neurthino
include("usefulFunctionsToPlot.jl")
include("NeurthinoRelated.jl")

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


#affichage des profils densité/distance pour n_vectors à partir d'un détecteur placé par l'utilisateur
iTime = 200
n_pts = 100
n_vectors = 7
zposition = 2.5e3 

densities ,sections = vectorsFromDetector(n_vectors, zposition)
export densities, sections


#==
function creationPaths(n_vectors;depthDetectorInM=2.5e3)

    dens, section = vectorsFromDetector(n_vectors) 

    paths = []
    for i in eachindex(dens)
        Path(dens[i],section[i])
        push!(paths, Path(dens[i],section[i]))
    end
    return paths
end


function linkWithNeurthino()
    osc = OscillationParameters(3)
    setθ!(osc, 1=>2, 0.59)
    setθ!(osc, 1=>3, 0.15)
    setθ!(osc, 2=>3, 0.84)
    setδ!(osc, 1=>3, 3.86)
    setΔm²!(osc, 2=>3, -2.523e-3)
    setΔm²!(osc, 1=>2, -7.39e-5)
    U = PMNSMatrix(osc)
    H = Hamiltonian(osc)


    paths = creationPaths(n_vectors)
    energies = 10 .^ range(0, stop=2, length=n_vectors)   
    prob = [Pνν(U, H, energies, path) for path in paths]


    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)
    hm=heatmap!(ax, energies, cos_θ, prob, colormap=inferno)
    Colorbar(fig[:,2], hm)
    display(fig)

    return energies, prob
end

energies, prob = linkWithNeurthino()
@show energies, prob





test

myAnimation(20,200,"temperature", temperatureFiles)
myPlot2DConvectionModel(200, "wtr", wtrFiles)
myPlot2DConvectionModel(200, "rho", rhoFiles)
myPlot2DConvectionModel(200, "temperature", temperatureFiles)
myPlot2DConvectionModel(200, "composition", compositionFiles)
lineDensityElectron2D([450,100], [350,150])
==#