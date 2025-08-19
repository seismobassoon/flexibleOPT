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
#dir="C:/Users/user/Desktop/stage 2A/données/MantleConvectionTakashi/op_old_full_mars_2025/"
rhoFiles=myListDir(dir; pattern=r"test_rho\d");
compositionFiles=myListDir(dir; pattern=r"test_c\d");
temperatureFiles=myListDir(dir; pattern=r"test_t\d");
wtrFiles=myListDir(dir; pattern=r"test_wtr\d");
#wtrFiles = filter(f -> !occursin(r"/\._", f), wtrFiles) #si données op_old_full_mars_2025



iTime = 200
n_pts = 100
n_vectors = 10
zposition = 2.5e3 

#vectorsFromDetector(n_vectors, zposition)

# Neurthino tests
function creationPaths(n_vectors, zposition)

    dens, section = vectorsFromDetector(n_vectors, zposition) 
    paths = Vector{Vector{Path}}(undef, n_vectors)  


    for i in 1:n_vectors
        paths[i]= [Path(dens[j],section[j]) for j in 1:n_pts-1]
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
    cos_θ = range(-1, 0, length = n_vectors)

    paths = creationPaths(n_vectors, zposition)
    energies = 10 .^ range(0, stop=2, length=n_vectors)
    #

    probsOld = collect(Pνν(U, H, energies, path) for path in paths)
    probs = Pνν(U, H, energies, paths)[:, :, 2, 2]
    #probs = [Pνν(U, H, energies, path)[1, 1, 2, 2] for path in paths]


    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1, xscale=log10, xlabel="Energy (GeV)", ylabel="cos(θ)")
    hm=heatmap!(ax, energies, cos_θ, probs, colormap=cgrad(:inferno))#, colorrange=(0,1))
    Colorbar(fig[:,2], hm, label="Probability")
    display(fig)

    return energies, probs,probsOld, paths
end


energies, probs,probsOld, paths=linkWithNeurthino()
export energies, probs, probsOld, paths






#==
zoa_value = 0.5
zoa = fill(zoa_value, size(fi)[1], size(fi)[2])
densityModified = densities * 2 * zoa
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