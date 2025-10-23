# this is generated completely by chatGPT


using Geodesy


"""
    moon, mars, venus, mercury :: Ellipsoid

IAU 2015 reference ellipsoids (mean radii or spheroids) for Solar System bodies.
Use them with Geodesy.jl like:

```julia

p = LLA(10.0, 45.0, 1000.0)           # latitude, longitude, altitude [m]
ecef = ECEFfromLLA(p, mars)
lla_back = LLAfromECEF(ecef, mars)
"""


# Planetary parameters (semi-major axis `a` in meters, flattening `f`)
mars_a = 3396200.0      # m
mars_f = 1 / 169.8

moon_a = 1737400.0      # m
moon_f = 0.0            # essentially spherical

venus_a = 6051800.0     # m
venus_f = 0.0           # very small flattening

mercury_a = 2439700.0   # m
mercury_f = 0.0         # spherical

function planet_ellipsoid(name::Symbol)
    if name == :Mars
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