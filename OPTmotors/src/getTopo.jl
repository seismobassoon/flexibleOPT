using GMT, Interpolations

include("../src/batchRevise.jl")
myInclude("../src/GeoPoints.jl")

function getTopo(p1::GeoPoint,p2::Geopoint,Δx::Float64,Δz::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64 = 0.0, rightLimit::Float64=(p2-p1).radius)
    # 2D
    x_axis,y_axis,z_axis,R=makeLocalCoordinates(p1,p2) # p1 centred coordinates
    Nx = Int64((rightLimit-leftLimit) ÷ Δx+1) 
    Nz = Int64((altMax-altMin) ÷ Δz + 1 ) 
end


function getTopo(Δx,Δz)
    #2D

    Δy = 100 # in metre as a dummy
    getTopo(Δx,Δy,Δz)
end

function getTopo(Δx,Δy,Δz)
    # 2D and 3D

end