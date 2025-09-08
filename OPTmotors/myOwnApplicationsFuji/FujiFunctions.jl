# FujiFunctions.jl

using Images, FileIO
using CairoMakie
using Interpolations
using GLMakie

"""
    loadFujiGrid(filename; vmin=2.0, vmax=6.0)
"""
function loadFujiGrid(filename; vmin=2.0, vmax=6.0)
    img = load(filename)
    gray = Gray.(img)
    vals = Float32.(gray)

    Vp = vmin .+ (vmax - vmin) .* vals
    Vs = 0.55 .* Vp
    ρ  = 2200 .+ 500 .* (Vp ./ vmax)

    nx, nz = size(Vp)
    x = range(0, stop=100e3, length=nx)
    z = range(0, stop=100e3, length=nz)

    return Vp, Vs, ρ, x, z
end


using CairoMakie

"""
    plotFujiModel2D(field, x, z; label="Vp (km/s)", cmap=:viridis)

Plot a Fuji 2D model (Vp, Vs, or ρ) as a heatmap.

- `field`: 2D array (nx × nz)
- `x`, `z`: coordinate vectors in meters
"""
function plotFujiModel2D(field, x, z; label="Vp (km/s)", cmap=:viridis)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="X (km)", ylabel="Z (km)", aspect=1)

    hm = heatmap!(ax, x ./ 1e3, z ./ 1e3, field', colormap=cmap)
    Colorbar(fig[:,2], hm, label=label)

    return fig, ax
end





using Interpolations, CairoMakie

"""
    FujiLineProfile2D(field, x, z, source, receiver; n_pts=200, label="Vp (km/s)", color=:red)

Sample a field (Vp, Vs, or ρ) along a straight line between source and receiver.

- `field`: 2D array (nx × nz)
- `x`, `z`: coordinate vectors (in meters)
- `source`, `receiver`: tuples (x,z) in meters
- `n_pts`: number of sampling points

Returns (values, sections, dist).
"""
function FujiLineProfile2D(field, x, z, source::Tuple, receiver::Tuple;
                           n_pts=200, label="Value", color=:red)

    # Interpolator
    itp = interpolate((x, z), field, Gridded(Linear()))

    # Discretize the line
    x_line = range(source[1], receiver[1], length=n_pts)
    z_line = range(source[2], receiver[2], length=n_pts)

    # Sample field values along the line
    values = [itp(xi, zi) for (xi, zi) in zip(x_line, z_line)]

    # Segment length in km
    seg_len_km = sqrt((x_line[2]-x_line[1])^2 + (z_line[2]-z_line[1])^2) * 1e-3
    sections = fill(seg_len_km, n_pts-1)
    dist = collect(0:seg_len_km:(n_pts-1)*seg_len_km)

    # Plot
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Distance (km)", ylabel=label)
    lines!(ax, dist, values, color=color)
    display(fig)

    return values, sections, dist
end




"""
    generateReceiversAtSurface(x_range; n_receivers=10, z=0.0)

Generate receiver positions evenly spaced along the surface (z=0 by default).
Returns a vector of (x,z) tuples in meters.
"""
function generateReceiversAtSurface(x_range; n_receivers=10, z=0.0)
    xs = range(first(x_range), last(x_range), length=n_receivers)
    receivers = [(xi, z) for xi in xs]
    return receivers
end


"""
    generateSourcesAtDepth(x_range, z_depth; n_sources=3)

Generate source positions evenly spaced along the bottom of the model.
"""
function generateSourcesAtDepth(x_range, z_depth; n_sources=3)
    xs = range(first(x_range), last(x_range), length=n_sources)
    sources = [(xi, z_depth) for xi in xs]
    return sources
end


"""
    FujiPathsFromPairs(field, x, z, pairs; n_pts=200, label="Vp (km/s)")

Compute profiles along multiple source–receiver pairs.
- `field`: 2D array (Vp, Vs, or ρ)
- `x`, `z`: coordinate vectors (meters)
- `pairs`: vector of ((xsrc,zsrc), (xrec,zrec)) tuples
- Returns (values_list, sections_list, dist_list)
"""
function FujiPathsFromPairs(field, x, z, pairs; n_pts=200, label="Vp (km/s)")
    values_list = []
    sections_list = []
    dist_list = []

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Distance (km)", ylabel=label)

    for (src, rec) in pairs
        vals, sections, dist = FujiLineProfile2D(field, x, z, src, rec;
                                                 n_pts=n_pts, label=label, color=:red)
        push!(values_list, vals)
        push!(sections_list, sections)
        push!(dist_list, dist)
    end

    display(fig)
    return values_list, sections_list, dist_list
end


struct SeismicPath
    values::Vector{Float64}   # sampled field values (Vp, Vs, ρ)
    sections::Vector{Float64} # segment lengths in km
    dist::Vector{Float64}     # cumulative distance along path
    source::Tuple{Float64,Float64}
    receiver::Tuple{Float64,Float64}
end

function FujiCreatePaths(field, x, z, pairs; n_pts=200)
    paths = SeismicPath[]
    for (src, rec) in pairs
        vals, sections, dist = FujiLineProfile2D(field, x, z, src, rec; n_pts=n_pts, label="Value")
        push!(paths, SeismicPath(vals, sections, dist, src, rec))
    end
    return paths
end


using GLMakie

"""
    interactiveFujiDemo(field, x, z; label="Vp (km/s)", cmap=:viridis)

Interactive demo:
- Left-click: place or move source (red star)
- Right-click: add receivers (blue triangles)
- Draws lines between source and receivers
- Returns Observables: (source, receivers)
"""
function interactiveFujiDemo(field, x, z; label="Vp (km/s)", cmap=:viridis)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="X (km)", ylabel="Z (km)", aspect=1)

    hm = heatmap!(ax, x ./ 1e3, z ./ 1e3, field', colormap=cmap)
    Colorbar(fig[:,2], hm, label=label)

    # Observables for source and receivers
    source = Observable(Point2f0(x[end] / 2e3, z[end] / 2e3)) # default center
    receivers = Observable(Point2f0[])

    # Markers
    scatter!(ax, lift(s -> [s], source), color=:red, marker=:star5, markersize=18)
    scatter!(ax, lift(r -> r, receivers), color=:blue, marker=:utriangle, markersize=14)

    # Draw connecting lines (source → each receiver)
    on(receivers) do recs
        delete!(ax.plots, 3:length(ax.plots)) # clear old lines
        for r in recs
            lines!(ax, [source[][1], r[1]], [source[][2], r[2]], color=:white, linewidth=1.5)
        end
    end

    # Mouse interactions
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            source[] = Point2f0(event.data[1], event.data[2]) # left click = source
        elseif event.button == Mouse.right && event.action == Mouse.press
            push!(receivers[], Point2f0(event.data[1], event.data[2])) # right click = add receiver
            notify(receivers)
        end
    end

    display(fig)
    return source, receivers
end
