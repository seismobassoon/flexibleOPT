using GMT, Interpolations

include("../src/batchRevise.jl")
myInclude("../src/GeoPoints.jl")

function getTopo(p1::GeoPoint,p2::Geopoint,Δx::Float64,Δz::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64 = 0.0, rightLimit::Float64=(p2-p1).radius)

    # 2D
    Δy = 100.0 # in metre as a dummy
    奥行きMin=0.0 # y axis range
    奥行きMax=0.0
    getTopo(p1,p2,Δx,Δy,Δz,奥行きMin,奥行きMax,altMin,altMax;leftLimit=leftLimit,rightLimit=rightLimit)
end


function getTopo(p1::GeoPoint,p2::GeoPoint,Δx::Float64,Δy::Float64,Δz::Float64,奥行きMin::Float64,奥行きMax::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64=0.0,rightLimit::Float64=(p2-p1).radius)

    x_axis,y_axis,z_axis,R=makeLocalCoordinates(p1,p2) # p1 centred coordinates

    Nx = Int64((rightLimit-leftLimit) ÷ Δx+1) 
    Ny = Int64((奥行きMax-奥行きMin) ÷ Δy +1)
    Nz = Int64((altMax-altMin) ÷ Δz + 1 ) 

    allGridsInGeoPoints=Array{GeoPoint,3}(undef,Nx,Ny,Nz)

    effectiveRadii=Array{Float64,3}(undef,Nx,Ny,Nz)

    allGridsInCartesian=Array{localCoord2D,3}(undef,Nx,Ny,Nz)

    seismicModel=(ρ=zeros(Float64,Nx,Ny,Nz),Vpv=zeros(Float64,Nx,Ny,Nz),Vph=zeros(Float64,Nx,Ny,Nz),Vsv=zeros(Float64,Nx,Ny,Nz),Vsh=zeros(Float64,Nx,Ny,Nz),Qμ=zeros(Float64,Nx,Ny,Nz),Qκ=zeros(Float64,Nx,Ny,Nz),QμPower=zeros(Float64,Nx,Ny,Nz),QκPower=zeros(Float64,Nx,Ny,Nz),η=zeros(Float64,Nx,Ny,Nz))
  


    for iXYZ in CartesianIndices(allGridsInGeoPoints)
        ix, iy, iz = Tuple(iXYZ)
        x = leftLimit+(ix-1)*Δx
        y = 奥行きMin+(iy-1)*Δy
        z = altMin+(iz-1)*Δz 

        tmpGeoPoint=GeoPoint(p_local_to_ECEF(x,z,p1.ecef,R))
        
        allGridsInGeoPoints[iXZ]=tmpGeoPoint

        allGridsInCartesian[iXZ]=localCoord2D(ix,iz,x,z)
    
        effectiveRadii[iXZ]=effectiveRadius(tmpGeoPoint,DSM1D.my1DDSMmodel.averagedPlanetRadiusInKilometer*1.e3 )
    end


    
end





function getTopo(Δx,Δy,Δz)
    # 2D and 3D

end
