module GeoPoints_hesa

#  GeoPoints_hesa.jl

using Geodesy, StaticArrays, LinearAlgebra
import Base: +, -, /, *

#  Default planet setup (Earth by default) 
DEFAULT_PLANET = Ref(:Earth)

function planet_ellipsoid(name::Symbol)
    if name == :Earth
        a = 6378137.0; f = 1.0 / 298.257223563
    elseif name == :Mars
        a = 3396200.0; f = 1 / 169.8
    elseif name == :Moon
        a = 1737400.0; f = 0.0
    elseif name == :Venus
        a = 6051800.0; f = 0.0
    elseif name == :Mercury
        a = 2439700.0; f = 0.0
    else
        error("Unknown planet: $name")
    end
    e2 = 2f - f^2
    b = a * (1 - f)
    return Ellipsoid(a, b, f, e2, name)
end

const DEFAULT_ELLIPSOID = Ref(planet_ellipsoid(DEFAULT_PLANET[]))

function set_default_planet!(name::Symbol)
    DEFAULT_PLANET[] = name
    DEFAULT_ELLIPSOID[] = planet_ellipsoid(name)
end


# Core structures 

struct GeoPoint
    lat::Float64   # degrees
    lon::Float64   # degrees
    alt::Float64   # meters
    ecef::SVector{3,Float64}
    radius::Float64 # meters
end

struct localCoord2D
    iXZ::SVector{2,Int}
    xz::SVector{2,Float64}
    horizontalVector::SVector{2,Float64}
    normalVector::SVector{2,Float64}
end

struct localCoord3D
    iXYZ::SVector{3,Int}
    xyz::SVector{3,Float64}
    horizontalVector1::SVector{3,Float64}
    horizontalVector2::SVector{3,Float64}
    normalVector::SVector{3,Float64}
end


#  Coordinate conversions 

function localCoord2D(p::localCoord3D)
    iXZ = SVector{2,Int}(p.iXYZ[1], p.iXYZ[3])
    xz  = SVector{2,Float64}(p.xyz[1], p.xyz[3])
    horiz = SVector{2,Float64}(p.horizontalVector1[1], p.horizontalVector1[3])
    normal = SVector{2,Float64}(p.normalVector[1], p.normalVector[3])
    return localCoord2D(iXZ, xz, horiz, normal)
end


# GeoPoint constructors 

function GeoPoint(lat::Float64, lon::Float64; alt=0.0, ell=DEFAULT_ELLIPSOID[])
    lla = LLA(lat, lon, alt)
    ecef = ECEF(lla, ell)
    radius = norm([ecef.x, ecef.y, ecef.z])
    return GeoPoint(lat, lon, alt, SVector(ecef.x, ecef.y, ecef.z), radius)
end

function GeoPoint(ecef::SVector{3,Float64}; ell=DEFAULT_ELLIPSOID[])
    lla = LLA(ECEF(ecef...), ell)
    radius = norm(ecef)
    return GeoPoint(lla.lat, lla.lon, lla.alt, ecef, radius)
end


#  Basic arithmetic for GeoPoints 

function +(a::GeoPoint, b::GeoPoint; ell=DEFAULT_ELLIPSOID[])
    ecef = a.ecef + b.ecef
    lla = LLA(ECEF(ecef...), ell)
    radius = norm(ecef)
    return GeoPoint(lla.lat, lla.lon, lla.alt, ecef, radius)
end

function -(a::GeoPoint, b::GeoPoint; ell=DEFAULT_ELLIPSOID[])
    ecef = a.ecef - b.ecef
    lla = LLA(ECEF(ecef...), ell)
    radius = norm(ecef)
    return GeoPoint(lla.lat, lla.lon, lla.alt, ecef, radius)
end

function /(a::GeoPoint, c::Real; ell=DEFAULT_ELLIPSOID[])
    ecef = a.ecef / c
    lla = LLA(ECEF(ecef...), ell)
    radius = norm(ecef)
    return GeoPoint(lla.lat, lla.lon, lla.alt, ecef, radius)
end

function *(a::GeoPoint, c::Real; ell=DEFAULT_ELLIPSOID[])
    ecef = a.ecef * c
    lla = LLA(ECEF(ecef...), ell)
    radius = norm(ecef)
    return GeoPoint(lla.lat, lla.lon, lla.alt, ecef, radius)
end


# Utilities 

function effectiveRadius(a::GeoPoint, r0::Float64; ell=DEFAULT_ELLIPSOID[])
    radiusPlanetHere = GeoPoint(a.lat, a.lon).radius
    ratio = r0 / radiusPlanetHere
    return a.radius * ratio
end

"""
    makeLocalCoordinates(p1::GeoPoint, p2::GeoPoint)

Return orthogonal (x, y, z) axes and rotation matrix R for the plane between p1 and p2.
"""
function makeLocalCoordinates(p1::GeoPoint, p2::GeoPoint;
                              pOrigin::GeoPoint=p1,
                              p0::GeoPoint=(p1 + p2)/2.0,
                              p2_1::GeoPoint=p2 - p1)
    return makeLocalCoordinates(p1.ecef, p2.ecef;
        pOrigin=pOrigin.ecef,
        p0=((p1+p2)/2.0).ecef,
        p2_1=(p2-p1).ecef)
end


function makeLocalCoordinates(p1::SVector{3,Float64}, p2::SVector{3,Float64};
                              pOrigin::SVector{3,Float64}=p1,
                              p0::SVector{3,Float64}=(p1+p2)/2.0,
                              p2_1::SVector{3,Float64}=p2-p1)
    x_axis_tentative = normalize(p2_1)
    z_axis = normalize(p0)
    y_axis = normalize(cross(z_axis, x_axis_tentative))
    x_axis = cross(y_axis, z_axis)
    R = SMatrix{3,3,Float64}(
        x_axis[1], x_axis[2], x_axis[3],
        y_axis[1], y_axis[2], y_axis[3],
        z_axis[1], z_axis[2], z_axis[3]
    )
    return x_axis, y_axis, z_axis, R
end


#Coordinate transforms 

p_ECEF_to_local(p_3D::SVector{3,Float64}, pOrigin::SVector{3,Float64}, R::SMatrix{3,3,Float64}) =
    R' * (p_3D - pOrigin)

p_local_to_ECEF(x_2D, z_2D, pOrigin::SVector{3,Float64}, R::SMatrix{3,3,Float64}) =
    pOrigin + R * SVector(x_2D, 0.0, z_2D)

p_local_to_ECEF(vec2D::SVector{2,Float64}, pOrigin::SVector{3,Float64}, R::SMatrix{3,3,Float64}) =
    pOrigin + R * SVector(vec2D[1], 0.0, vec2D[2])

p_local_to_ECEF(x_3D, y_3D, z_3D, pOrigin::SVector{3,Float64}, R::SMatrix{3,3,Float64}) =
    pOrigin + R * SVector(x_3D, y_3D, z_3D)

p_local_to_ECEF(vec3D::SVector{3,Float64}, pOrigin::SVector{3,Float64}, R::SMatrix{3,3,Float64}) =
    pOrigin + R * vec3D


# Custom helper examples

function midpoint(p1::GeoPoint, p2::GeoPoint)
    return (p1 + p2) / 2.0
end

function distance_km(p1::GeoPoint, p2::GeoPoint)
    return norm(p1.ecef - p2.ecef) / 1e3
end

end # module GeoPoints_hesa
