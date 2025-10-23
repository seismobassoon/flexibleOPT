# some useful structures and functions to deal with planetary 3D coordinates


using Geodesy, StaticArrays,LinearAlgebra

import Base: +,-,/,*



struct GeoPoint
    lat::Float64 # in degree
    lon::Float64 # in degree
    alt::Float64 # in metre
    ecef::SVector{3,Float64}
    radius::Float64 # in metre
    #effectiveRadius::Float64 # in metre (useful to get values from 1D averaged model)
end

struct localCoord2D
    ix::Int64
    iz::Int64
    x::Float64
    z::Float64
    horizontalVector::SVector{2,Float64}
    normalVector::SVector{2,Float64}
end

struct localCoord3D
    ix::Int64
    iy::Int64
    iz::Int64
    x::Float64
    y::Float64
    z::Float64
    horizontalVector1::SVector{3,Float64}
    horizontalVector2::SVector{3,Float64}
    normalVector::SVector{3,Float64}
end

function localCoord2D(ix,iz,x,z)
    localCoord2D(ix,iz,x,z,(0.0,0.0),(0.0,0.0))
end

function GeoPoint(lat::Float64, lon::Float64; alt=0.0, ell=wgs84, planet="Earth")
    lla = LLA(lat,lon, alt) # be careful LLA uses degrees by default!!
    ecef_coords = ECEF(lla,ell)
    radius = norm([ecef_coords.x,ecef_coords.y,ecef_coords.z])
    GeoPoint(lat, lon, alt, SVector(ecef_coords.x, ecef_coords.y, ecef_coords.z),radius)
end


function GeoPoint(ecef::SVector{3,Float64}; ell=wgs84)
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function +(a::GeoPoint,b::GeoPoint; ell=wgs84)
    ecef=a.ecef + b.ecef
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function -(a::GeoPoint,b::GeoPoint; ell=wgs84)
    ecef=a.ecef - b.ecef
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function /(a::GeoPoint,c::Real; ell=wgs84)
    ecef=a.ecef / c
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end


function *(a::GeoPoint,c::Real; ell=wgs84)
    ecef=a.ecef * c
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end


function effectiveRadius(a::GeoPoint,r0::Float64; ell=wgs84)
    radiusPlanetHere = GeoPoint(a.lat,a.lon).radius 
    ratio = r0/radiusPlanetHere
    return a.radius*ratio
end


makeLocalCoordinates(p1::GeoPoint,p2::GeoPoint;pOrigin::GeoPoint=p1,p0::GeoPoint=(p1+p2)/2.0, p2_1::GeoPoint=p2-p1)=makeLocalCoordinates(p1.ecef,p2.ecef;pOrigin=pOrigin.ecef,p0=((p1+p2)/2.0).ecef,p2_1=(p2-p1).ecef)

makeLocalCoordinates(p1::GeoPoint,p2::GeoPoint;pOrigin::SVector{3,Float64},p0::SVector{3,Float64}, p2_1::SVector{3,Float64})=makeLocalCoordinates(p1.ecef,p2.ecef;pOrigin=pOrigin,p0=p0,p2_1=p2_1)


function makeLocalCoordinates(p1::SVector{3,Float64},p2::SVector{3,Float64};pOrigin::SVector{3,Float64}=p1, p0::SVector{3,Float64}=(p1+p2)/2.0, p2_1::SVector{3,Float64}=p2-p1)

    # this function will define local (x,y,z) coordinates centred at pOrigin
    # 
    # z axis should be parallel to p0 
    # y axis is the normal to the plane that is defined by p0 and p2_1
    # x axis is the vector on the plane that is p0 and p2_1 which can be similar to p2_1 direction
    

    x_axis_tentative = normalize(p2_1)
    z_axis = normalize(p0)
    # y-axis: complete right-handed system
    y_axis = normalize(cross(z_axis, x_axis_tentative))
    x_axis = cross(y_axis, z_axis)  # now perfectly orthogonal
   
    #Rotation matrix
    R = SMatrix{3,3,Float64}(
        x_axis[1], x_axis[2], x_axis[3],
        y_axis[1], y_axis[2], y_axis[3],
        z_axis[1], z_axis[2], z_axis[3]
    )

end


p_2D_to_ECEF(x_2D,z_2D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_2D,0.e0,z_2D)

p_ECEF_to_2D(p_3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = R' * (p_3D - pOrigin)

