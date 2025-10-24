using Statistics, Printf, DelimitedFiles, FileIO, Colors
import GLMakie
using GLMakie: Figure, Axis, Point2f, hidespines!, xlims!, ylims!, axislegend, rowgap!, cgrad
using Interpolations  # NEW

GLMakie.activate!()
try
    GLMakie.set_window_config!(; samples=8)  # MSAA for smoother edges (if your GL version supports it)
catch
    @warn "Could not set samples=8; falling back to default anti-aliasing."
end

#  read (nz,nx) in Fortran order 
function read_binary_matrix(path, nz, nx)
    vec = Array{Float32}(undef, nz * nx)
    open(path, "r") do io
        read!(io, vec)
    end
    reshape(vec, nx, nz)'  # (nz, nx)
end

# smooth upsample to a finer (z,x) grid
function upsample_field(field::AbstractMatrix{T}, z::AbstractVector, x::AbstractVector; factor::Int=2) where {T<:Real}
    factor ≤ 1 && return field, z, x
    zfine = range(first(z), last(z); length=length(z)*factor)
    xfine = range(first(x), last(x); length=length(x)*factor)
    itp = interpolate((z, x), field, Gridded(Linear()))
    # evaluate on fine grid
    F = [itp(zz, xx) for zz in zfine, xx in xfine]
    return F, zfine, xfine
end


function main()

    
    dir = "/Users/hessiemohammadi/Documents/github/OptimallyAccurate2D/volcano_simu1/rock_magma_chamber_topo_air/bin_files"
    datadir = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/usefulScriptsForFriends"

    # Visualization parameters
    UPSAMPLE_FACTOR = 3          # 1 = off; 2–4 looks great (higher = slower)
    FIG_RES = (1600, 850)        # pixel size
    PX_PER_UNIT = 3              
    VIDEO_FPS = 30               # output frame rate
    VIDEO_CRF = 18               # video quality (lower = better; 18–23 is good range)
    WAVE_ALPHA = 0.85

    # Model grid sizes 
    nx_full, nz_full = 261, 438
    nx_crop, nz_crop = 161, 338
    n_expected = nx_full * nz_full

    #  Physical coordinates (km) 
    x_min, x_max = -50.0, 50.0
    z_min, z_max = -10.0, 43.70
    x_full = range(x_min, x_max; length=nx_full)
    z_full = range(z_min, z_max; length=nz_full)

    # Cropping indices (center crop) 
    ix_start = Int(floor((nx_full - nx_crop) / 2)) + 1
    ix_end   = ix_start + nx_crop - 1
    iz_start = Int(floor((nz_full - nz_crop) / 2)) + 1
    iz_end   = iz_start + nz_crop - 1

    #  Cropped axes 
    x = x_full[ix_start:ix_end]
    z_data_min = -3.776
    z = range(z_data_min, z_max; length=nz_crop)

    #  Load background Vp 
    vp_file = joinpath(datadir, "rock_magma_chamber_topo_air.vp")
    println("Loading velocity model: $vp_file")
    vp = read_binary_matrix(vp_file, nz_crop, nx_crop)
    vp_min, vp_max = extrema(vp)

    # Color maps 
    my_colormap = cgrad([
        RGB(1.0, 1.0, 1.0),   # air
        RGB(1.0, 0.0, 0.0),   # magma
        RGB(1.0, 0.65, 0.0),  # mush
        RGB(0.55, 0.27, 0.07) # rock
    ], [0.0, 0.33, 0.66, 1.0])

    my_wave_colormap = cgrad([RGB(0,0,1), RGB(1,1,1), RGB(1,0,0)])

    #  List binary wavefield files 
    files = sort(filter(f -> endswith(f, ".bin"), readdir(dir; join=true)))

    # Source location (indices and coordinates)
    iz_src, ix_src = 133, 60
    x_src, z_src = x[ix_src], z[iz_src]
    println("Source indices (ix, iz) = ($ix_src, $iz_src)")
    println(" Source coordinates (x, z) = ($(round(x_src, digits=2)), $(round(z_src, digits=2))) km")

    # Receivers
    receiver_coords = [(-10.0, 1.0), (-5.0, 1.0), (0.0, 1.0), (10.0, 1.0), (20.0, 1.0)]
    receiver_indices = [(argmin(abs.(x .- xr)), argmin(abs.(z .- zr))) for (xr, zr) in receiver_coords]
    time_series = [Float32[] for _ in receiver_indices]

    # Global amplitude color range (using compressed field) 
    global_max = 0.0
    for file in files
        raw = read(file)
        data = reinterpret(Float32, raw)
        if length(data) == n_expected
            data = reshape(data, (nx_full, nz_full))'
            data = data[iz_start:iz_end, ix_start:ix_end]
            c = 1e-4
            B = sign.(data) .* log1p.(abs.(data) ./ c)
            global_max = max(global_max, maximum(abs, B))
        end
    end
    println(" Global color scale range: ±$(round(global_max, digits=4))")

    #  Pre-upsample the static Vp for a crisp background 
    vp_plot, z_plot, x_plot = upsample_field(vp, z, x; factor=UPSAMPLE_FACTOR)

    #  Process each wavefield snapshot 
    for (i, file) in enumerate(files)
        println("Processing: $(i)/$(length(files)) - $(basename(file))")

        raw = read(file)
        data = reinterpret(Float32, raw)
        if length(data) != n_expected
            @warn "Skipping $(basename(file)): length $(length(data)) ≠ $n_expected"
            continue
        end

        data = reshape(data, (nx_full, nz_full))'
        data = data[iz_start:iz_end, ix_start:ix_end]

        # amplitude compression
        c = 1e-4
        B = sign.(data) .* log1p.(abs.(data) ./ c)

        # upsample the wavefield for smooth display
        B_plot, zf, xf = upsample_field(B, z, x; factor=UPSAMPLE_FACTOR)

        # figure
        fig = Figure(resolution=FIG_RES)
        rowgap!(fig.layout, 0)
        ax = Axis(fig[1, 1], xlabel="x (km)", ylabel="z (km)",
                  xgridvisible=false, ygridvisible=false)
        hidespines!(ax)
        xlims!(ax, first(x), last(x))
        ylims!(ax, 43.70, -10)

        # z ticks
        custom_ticks = sort(unique(vcat(collect(-10:10:50), [z_data_min])))
        ax.yticks = (custom_ticks, string.(round.(custom_ticks, digits=3)))

        # air layer
        GLMakie.poly!(ax, Point2f[(x_min, -10.0), (x_max, -10.0), (x_max, z_data_min), (x_min, z_data_min)], color=:white)

        # Vp (upsampled) — use interpolate=true for smooth look
        GLMakie.heatmap!(ax, x_plot, z_plot, vp_plot;
            colormap=my_colormap, colorrange=(vp_min, vp_max), alpha=1.0, interpolate=true)

        # wavefield overlay (upsampled)
        hm = GLMakie.heatmap!(ax, xf, zf, B_plot;
            colormap=my_wave_colormap, colorrange=(-global_max, global_max),
            alpha=WAVE_ALPHA, interpolate=true)

        # source marker (reproject to fine grid coords visually)
        GLMakie.scatter!(ax, [x_src], [z_src];
            color=:red, marker=:star5, markersize=12, strokewidth=1.25, label="Source")

        # receivers
        for (ri, (xr, zr)) in enumerate(receiver_coords)
            GLMakie.scatter!(ax, [xr], [zr];
                color=:green, marker=:utriangle, rotation=π,
                markersize=13, strokewidth=1.2, label="Receiver $(ri)")
            GLMakie.text!(ax, xr + 0.5, zr, text=string(ri),
                align=(:left, :center), color=:black, fontsize=11)
        end

        # colorbars
        GLMakie.Colorbar(fig[1, 2], hm, label="Amplitude",
            colorrange=(-global_max, global_max))
        transparent_colormap = [RGBA(c.r, c.g, c.b, 0.7) for c in my_colormap.colors]
        GLMakie.Colorbar(fig[1, 0];
            limits=(vp_min, vp_max),
            colormap=GLMakie.cgrad(transparent_colormap),
            label="Vp (km/s)")

        axislegend(ax; position=:rb, orientation=:vertical, patchsize=(8,8), labelsize=10)

        # save PNG (HiDPI)
        out_png = replace(file, ".bin" => ".png")
        GLMakie.save(out_png, fig; px_per_unit=PX_PER_UNIT)
        println("Saved snapshot to: $out_png")

        # accumulate native (non-upsampled) receiver samples
        for (ri, (ix_r, iz_r)) in enumerate(receiver_indices)
            push!(time_series[ri], data[iz_r, ix_r])
        end
    end

    println("All snapshots saved")

    # Save time series 
    ts_dir = joinpath(dir, "time_series")
    isdir(ts_dir) || mkpath(ts_dir)
    Δt = 0.005
    t = collect(0:Δt:(length(files)-1)*Δt)

    for (ri, (xr, zr)) in enumerate(receiver_coords)
        series = time_series[ri]
        out_ascii = joinpath(ts_dir, @sprintf("receiver_%02d_x%.1f_z%.1f.txt", ri, xr, zr))
        writedlm(out_ascii, hcat(t, series))
        println("Saved time series to: $out_ascii")

        fig = Figure(resolution=(900, 360))
        ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Amplitude",
                  title=@sprintf("Receiver %d (x=%.1f km, z=%.1f km)", ri, xr, zr))
        GLMakie.lines!(ax, t, series; linewidth=1.6)
        GLMakie.save(replace(out_ascii, ".txt" => ".png"), fig; px_per_unit=2)
    end

    println("All receiver time series saved in $ts_dir")

    #  Make high-quality video 
    output_video = joinpath(dir, "wavefield.mp4")
    cd(dir) do
        cmd = `ffmpeg -framerate $VIDEO_FPS -pattern_type glob -i "*.png" \
            -c:v libx264 -preset slow -crf $VIDEO_CRF -pix_fmt yuv420p \
            -vf "fps=$VIDEO_FPS,scale=trunc(iw/2)*2:trunc(ih/2)*2" $output_video`
        run(cmd)
        println("Video saved to: $output_video")
    end
end

main()
