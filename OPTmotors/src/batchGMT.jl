# useful functions for GMT.jl

using GMT

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

function GeoPoint(lat::Float64, lon::Float64; alt=0.0, ell=wgs84)
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




