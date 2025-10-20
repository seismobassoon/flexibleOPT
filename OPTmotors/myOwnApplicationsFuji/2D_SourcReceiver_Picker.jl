using GLMakie, Colors
using Makie: Axis
using Makie: Top
using DelimitedFiles
using CSV, DataFrames 

"""
2D interactive figure to pick Sources (S) and Receivers (R)
 Click to add S (red circles) and R (green inverted triangles, z ≤ 0)
 Next button: switch from S to R, then save the figure and export points
 Toggle Pair mode: connect S→R by clicking near them (nearest point is picked)
 Background: velocity model heatmap (Vp)
"""
fig = Figure(size = (800, 500))
ax  = Axis(fig[1, 1], xlabel = "x (km)", ylabel = "z (km)")
ax.aspect = DataAspect()  # Ensures equal scale in x and z

xlims!(ax, -50.0, 50.0)
ylims!(ax,  50.0, -10.0)     # surface near top

# binary model files (from usefulImageReaderForHesaneh.jl output)

datadir = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/usefulScriptsForFriends"

nz, nx = (338, 461)   # model dimensions
#nz, nx = (389, 206)
#nz, nx = (427, 461)
npts = nz * nx
  

# Read topography file 
topofile = "/Users/hessiemohammadi/Documents/FUJI/events/surface_data_new_offsets.csv"
df = CSV.read(topofile, DataFrame)

x_topo = df.x_offset_km
z_topo = df.elevation_m ./ 1000.0   # convert m → km

# Extend topography to match x-axis limits
xlims_plot = (-50.0, 50.0)


# Overlay topography curve
lines!(ax, x_topo, z_topo; color=RGB(0.55, 0.27, 0.07), linewidth=2, label="Topography")



function read_binary_matrix(path, nz, nx)
    vec = Array{Float32}(undef, nz*nx)
    open(path, "r") do io
        read!(io, vec)
    end
    
    return reshape(vec, nx, nz)'  
    
end


rho = read_binary_matrix(joinpath(datadir, "volcano_f_vp300.rho"), nz, nx)
vp  = read_binary_matrix(joinpath(datadir, "volcano_f_vp300.vp"),  nz, nx)
vs  = read_binary_matrix(joinpath(datadir, "volcano_f_vp300.vs"),  nz, nx)

# Coordinates
xs = range(50, -50, length=nx)     # horizontal (km)
zs = range(start = -3.776, stop = 50, length = nz)
     # vertical (km)


# custom gradient (low → high)
my_colormap = cgrad([
    RGB(1.0, 1.0, 1.0),     # Air (white)  
    RGB(1.0, 0.0, 0.0),  # red = low velocity
    RGB(1.0, 0.65, 0.0),     # orange = medium velocity
    RGB(0.55, 0.27, 0.07)       # brown = high velocity
], [0.0, 0.5, 1.0])          # control points (low, mid, high)

# Background heatmap (Vp)
#hm = heatmap!(ax, xs, zs, reverse(vp, dims=1); colormap=my_colormap, alpha=0.8)
hm = heatmap!(ax, xs, zs, reverse(vp, dims=1); 
              colormap=my_colormap, alpha=0.8)


Colorbar(fig[1, 2], hm, label = "Vp (m/s)", height = Relative(0.6))


#  Data stores 
sources   = Observable(Point2f[])   # red circles
receivers = Observable(Point2f[])   # green inverted triangles
s_count   = Ref(0)
r_count   = Ref(0)
pairs     = Tuple{Int,Int}[]        # store S–R pairs

#plots
scatter!(ax, sources;   color=:red,   markersize=15, marker=:circle)
scatter!(ax, receivers; color=:green, markersize=16, marker=:dtriangle)

#  Modes / state
mode      = Ref(:S)      # :S => add Sources, :R => add Receivers
saved     = Ref(false)   # set after figure is saved (pair mode still works)
pair_mode = Ref(false)   # connect S→R when true

# Files
outfile = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/myOwnApplicationsFuji/clicked_points.txt"
pngfile = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/myOwnApplicationsFuji/final_points.png"

#status overlay
status = Observable("Mode: S (Sources). Click to add | Next → Receivers (z ≤ 0)")
Label(fig[1, 1, Top()], status; fontsize=12, halign=:center, padding=(6,10,6,10), tellwidth=false)

# Controls (buttons)
fig[2, 1] = grid = GridLayout()
btn_next  = Button(grid[1, 1], label = "Next (S→R→Save)")
tgl_pair  = Toggle(grid[1, 2], active = false)
Label(grid[1, 3], "Pair mode")


function nearest_index_any(pts::Vector{Point2f}, q::Point2f)
    isempty(pts) && return 0, Inf
    best = 1
    dmin = hypot(pts[1][1]-q[1], pts[1][2]-q[2])
    @inbounds for i in 2:length(pts)
        pt = pts[i]
        d  = hypot(pt[1]-q[1], pt[2]-q[2])
        if d < dmin
            dmin = d; best = i
        end
    end
    return best, dmin
end

function write_full_export!(path::AbstractString)
    open(path, "w") do io
        println(io, "# Sources:")
        for i in 1:length(sources[])
            p = sources[][i]
            println(io, "S$i: x=$(round(p[1]; digits=2)), z=$(round(p[2]; digits=2))")
        end
        println(io, "")
        println(io, "# Receivers (z ≤ 0):")
        for j in 1:length(receivers[])
            p = receivers[][j]
            println(io, "R$j: x=$(round(p[1]; digits=2)), z=$(round(p[2]; digits=2))")
        end
        println(io, "")
        println(io, "# Pairs:")
        for (i, j) in pairs
            println(io, "PAIR: S$i -> R$j")
        end
    end
end

function draw_pair!(i::Int, j::Int)
    S = sources[][i]; R = receivers[][j]
    linesegments!(ax, [S, R]; color=:black, linewidth=2)
    mid = Point2f((S[1]+R[1])/2, (S[2]+R[2])/2)
    text!(ax, "S$(i)→R$(j)", position=mid, align=(:center, :bottom), offset=(0,6))
    push!(pairs, (i, j))
    write_full_export!(outfile)  
    status[] = "Pair added: S$(i)→R$(j). Toggle Pair mode to add more."
end

# Selection markers
selS_idx = Ref(0)
selR_idx = Ref(0)
selS_pt  = Observable([Point2f(NaN, NaN)])
selR_pt  = Observable([Point2f(NaN, NaN)])
scatter!(ax, selS_pt; color=:transparent, marker=:circle,    markersize=26, strokewidth=3, strokecolor=:orange)
scatter!(ax, selR_pt; color=:transparent, marker=:dtriangle, markersize=28, strokewidth=3, strokecolor=:cyan)

click_probe = Observable([Point2f(NaN, NaN)])
scatter!(ax, click_probe; color=:black, markersize=8, marker=:x)

function clear_selection_markers!()
    selS_pt[] = [Point2f(NaN, NaN)]
    selR_pt[] = [Point2f(NaN, NaN)]
    notify(selS_pt); notify(selR_pt)
end

# Mouse clicks 
on(events(fig).mousebutton) do ev
    if ev.button == Mouse.left && ev.action == Mouse.press
        pos = mouseposition(ax)
        pos === nothing && return
        p = Point2f(Float32(pos[1]), Float32(pos[2]))
        click_probe[] = [p]; notify(click_probe)

        if pair_mode[]
            if isempty(sources[]) || isempty(receivers[])
                status[] = "Add at least one Source and one Receiver first."
                return
            end
            if selS_idx[] == 0
                i, dS = nearest_index_any(sources[], p)
                selS_idx[] = i
                selS_pt[] = [sources[][i]]; notify(selS_pt)
                status[] = "Selected S$(i) (≈$(round(dS; digits=2)) km). Now click near a Receiver."
            elseif selR_idx[] == 0
                j, dR = nearest_index_any(receivers[], p)
                selR_idx[] = j
                selR_pt[] = [receivers[][j]]; notify(selR_pt)
                draw_pair!(selS_idx[], selR_idx[])
                selS_idx[] = 0
                selR_idx[] = 0
                clear_selection_markers!()
            end
            return
        end

        if saved[]
            return
        end
        if mode[] == :S
            s_count[] += 1
            name = "S$(s_count[])"
            push!(sources[], p); notify(sources)
            text!(ax, name, position=p, align=(:left, :bottom), offset=(5,5))
            write_full_export!(outfile)
            status[] = "Added $name  |  Next → Receivers (z ≤ 0)"
        elseif mode[] == :R
            if p[2] > 0
                status[] = "Receiver rejected: z=$(round(p[2]; digits=2)) > 0 km (must be ≤ 0)."
                return
            end
            r_count[] += 1
            name = "R$(r_count[])"
            push!(receivers[], p); notify(receivers)
            text!(ax, name, position=p, align=(:left, :bottom), offset=(5,5))
            write_full_export!(outfile)
            status[] = "Added $name  |  Next → SAVE (then Pair mode)"
        end
    end
end

# Button actions 
on(btn_next.clicks) do _
    if mode[] == :S
        mode[] = :R
        status[] = "Mode: R (Receivers). Click to add (z ≤ 0) | Next → SAVE"
    elseif mode[] == :R
        save(pngfile, fig)
        write_full_export!(outfile)
        saved[] = true
        status[] = "Figure saved → $pngfile and log exported → $(outfile). Toggle Pair mode to connect S→R."
    end
end

on(tgl_pair.active) do v
    pair_mode[] = v
    selS_idx[] = 0
    selR_idx[] = 0
    clear_selection_markers!()
    if v
        status[] = "PAIR MODE: click near a Source, then near a Receiver."
    else
        status[] = "Exited Pair mode."
    end
end

fig
