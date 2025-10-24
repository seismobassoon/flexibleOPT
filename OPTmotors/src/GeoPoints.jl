# some useful structures and functions to deal with planetary 3D coordinates


using Geodesy, StaticArrays,LinearAlgebra

import Base: +,-,/,*

DEFAULT_PLANET = Ref(:Earth) #set_default_planet! can change this





# Planetary parameters (semi-major axis `a` in meters, flattening `f`)

function planet_ellipsoid(name::Symbol)
    if name == :Earth # WGS84
        a= 6378137.0
        f =1.0/298.257223563
    elseif name == :Mars
        a = 3396200.0
        f = 1 / 169.8
    elseif name == :Moon
        a = 1737400.0
        f = 0.0
    elseif name == :Venus
        a = 6051800.0
        f = 0.0
    elseif name == :Mercury
        a = 2439700.0
        f = 0.0
    else
        error("Unknown planet: $name")
    end
    e2 = 2f - f^2
    b = a * (1 - f)
    return Ellipsoid(a, b, f, e2, name)
end


const DEFAULT_ELLIPSOID = Ref(planet_ellipsoid(DEFAULT_PLANET[]))

# The planet can be changed by this function

function set_default_planet!(name::Symbol)
    DEFAULT_PLANET[]=name
    DEFAULT_ELLIPSOID[] = planet_ellipsoid(name)
end


struct GeoPoint
    lat::Float64 # in degree
    lon::Float64 # in degree
    alt::Float64 # in metre
    ecef::SVector{3,Float64}
    radius::Float64 # in metre
    #effectiveRadius::Float64 # in metre (useful to get values from 1D averaged model)
end

struct localCoord2D
    iXZ::SVector{2,Integer}
    xz::SVector{2,Float64}
    horizontalVector::SVector{2,Float64}
    normalVector::SVector{2,Float64}
end

struct localCoord3D
    iXYZ::SVector{3,Integer}
    xyz::SVector{3,Float64}
    horizontalVector1::SVector{3,Float64}
    horizontalVector2::SVector{3,Float64}
    normalVector::SVector{3,Float64}
end

function localCoord2D(p::localCoord3D)
    iXZ=SVector{2,Integer}(p.iXYZ[1],p.iXYZ[3])
    xz=SVector{2,Float64}(p.xyz[1],p.xyz[3])
    horizontalVector=SVector{2,Float64}(p.horizontalVector1[1],p.horizontalVector1[3])
    normalVector=SVector{2,Float64}(p.normalVector[1],p.normalVector[3])
    localCoord2D(iXZ,xz,horizontalVector,normalVector)
end

function localCoord3D(p::localCoord2D)
    iXYZ=SVector{3,Integer}(p.iXYZ[1],1,p.iXYZ[2])
    xyz=SVector{3,Float64}(p.xyz[1],0.0,p.xyz[2])
    horizontalVector=SVector{3,Float64}(p.horizontalVector1[1],0.0,p.horizontalVector1[2])
    horizontalVector=SVector{3,Float64}(0.0,1.0,0.0)
    normalVector=SVector{3,Float64}(p.normalVector[1],0.0,p.normalVector[2])
    localCoord3D(iXYZ,xyz,horizontalVector1,horizontalVector2,normalVector)
end


function localCoord2D(iX::SVector{2,Integer},X::SVector{2,Float64})
    
end


function GeoPoint(lat::Float64, lon::Float64; alt=0.0, ell=DEFAULT_ELLIPSOID[])
    lla = LLA(lat,lon, alt) # be careful LLA uses degrees by default!!
    ecef_coords = ECEF(lla,ell)
    radius = norm([ecef_coords.x,ecef_coords.y,ecef_coords.z])
    GeoPoint(lat, lon, alt, SVector(ecef_coords.x, ecef_coords.y, ecef_coords.z),radius)
end

function GeoPoint(x::Real,y::Real,z::Real; ell=DEFAULT_ELLIPSOID[])
    ecef = SVector{3,Float64}(x,y,z)
    GeoPoint(ecef;ell=ell)
end

function GeoPoint(ecef::SVector{3,Float64}; ell=DEFAULT_ELLIPSOID[])
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function +(a::GeoPoint,b::GeoPoint; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef + b.ecef
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function -(a::GeoPoint,b::GeoPoint; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef - b.ecef
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function /(a::GeoPoint,c::Real; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef / c
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end


function *(a::GeoPoint,c::Real; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef * c
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end


function effectiveRadius(a::GeoPoint,r0::Float64; ell=DEFAULT_ELLIPSOID[])
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
    return x_axis,y_axis,z_axis,R
end

p_ECEF_to_local(p_3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = R' * (p_3D - pOrigin)

p_local_to_ECEF(x_2D,z_2D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_2D,0.e0,z_2D)

p_local_to_ECEF(vec2D::SVector{2,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(vec2D[1],0.e0,vec2D[2])

p_local_to_ECEF(x_3D,y_3D,z_3D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_3D,y_3D,z_3D)

p_local_to_ECEF(vec3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(vec3D[1],vec3D[2],vec3D[3])