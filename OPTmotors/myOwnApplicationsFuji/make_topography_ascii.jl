using CSV, DataFrames, Statistics, Plots

# ---- parameters: MUST match your Fortran model ----
csvfile = "/Users/hessiemohammadi/Documents/FUJI/events/surface_data_new_offsets.csv"
outfile = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/myOwnApplicationsFuji/topography_profile_ix_iz.dat"
plotfile = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/myOwnApplicationsFuji/topography_profile_ix_iz.png"

dx = 0.296     # km
dz = 0.296     # km
nx = 337       # horizontal (ix) — elevation_m
nz = 160       # vertical (iz) — x_offset_km
xmin, xmax = -50.0, 50.0   # range for x_offset_km filtering

# ---- load and filter data ----
df = CSV.read(csvfile, DataFrame)
filter!(row -> xmin <= row.x_offset_km <= xmax, df)

# convert elevation (m, upwards) to depth (km, positive down)
df.depth_km = -df.elevation_m ./ 1000.0

# keep only relevant columns and ensure sorted by x_offset_km
prof = unique(df[:, [:x_offset_km, :depth_km]])
sort!(prof, :x_offset_km)

# ---- regular grid along x_offset_km (vertical direction now) ----
x_min = minimum(prof.x_offset_km)
x_max = maximum(prof.x_offset_km)
xs = range(x_min, x_max, length=nz)   # nz vertical indices

# interpolate or map elevation/depth to the same x grid
function nearest_depth(xq)
    i = findmin(abs.(prof.x_offset_km .- xq))[2]
    return prof.depth_km[i]
end
z_on_rows = [nearest_depth(x) for x in xs]

# ---- regular grid along elevation (horizontal) ----
zmin = minimum(z_on_rows)
zmax = maximum(z_on_rows)
zs = range(zmin, zmax, length=nx)   # nx horizontal columns

# ---- map elevation (depth_km) to horizontal index ix
function nearest_elev(zq)
    i = findmin(abs.(z_on_rows .- zq))[2]
    return prof.x_offset_km[i]
end
x_on_cols = [nearest_elev(z) for z in zs]

# ---- convert to indices ----
ix = round.(Int, (z_on_rows .- zmin) ./ dx) .+ 1
iz = round.(Int, (xs .- x_min) ./ dz) .+ 1
ix = clamp.(ix, 1, nx)
iz = clamp.(iz, 1, nz)

# ---- write (ix, iz) table ----
open(outfile, "w") do io
    for i in 1:length(xs)
        println(io, ix[i], " ", iz[i])
    end
end
println("Wrote ", outfile, " with nz=", nz, " entries (one iz per ix).")

# ---- plot the topography ----
p = plot(
    prof.depth_km, prof.x_offset_km,
    xlabel = "Depth (km, positive down)",
    ylabel = "Horizontal distance (km)",
    title = "Topography (ix=elevation, iz=x_offset_km)",
    lw = 3, color = :brown, label = "Surface",
    xlim = (minimum(prof.depth_km), maximum(prof.depth_km)),
    ylim = (xmin, xmax),
    legend = :topleft,
    fillrange = maximum(prof.x_offset_km),
    fillalpha = 0.3,
    fillcolor = :brown,
    aspect_ratio = :equal
)

savefig(p, plotfile)
println("Saved plot: ", plotfile)
