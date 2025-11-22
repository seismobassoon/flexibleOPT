using Statistics, Printf, DelimitedFiles, FileIO, Colors
import GLMakie
using GLMakie: Figure, Axis, xlims!, ylims!, axislegend, rowgap!, cgrad
using Interpolations

GLMakie.activate!()
try
    GLMakie.set_window_config!(; samples=8)
catch
    @warn "Could not set samples=8; falling back to default anti-aliasing."
end

# --- read (nx,nz) exactly as written (no transpose) ---
function read_binary_matrix(path, nx, nz)
    vec = Array{Float32}(undef, nx * nz)
    open(path, "r") do io
        read!(io, vec)
    end
    reshape(vec, nx, nz)   # shape = (ix, iz)
end

# --- smooth upsample to a finer (x,z) grid ---
function upsample_field(field::AbstractMatrix{T}, x::AbstractVector, z::AbstractVector; factor::Int=2) where {T<:Real}
    factor ≤ 1 && return field, x, z
    xfine = range(first(x), last(x); length=length(x)*factor)
    zfine = range(first(z), last(z); length=length(z)*factor)
    itp = interpolate((x, z), field, Gridded(Linear()))
    F = [itp(xx, zz) for xx in xfine, zz in zfine]
    return F, xfine, zfine
end

function main()
    # === I/O ===
    dir     = "/Users/hessiemohammadi/Documents/github/OptimallyAccurate2D/volcano_simu1/snapshots_rock1/bin_files"
    frames  = joinpath(dir, "frames_highres")
    ts_dir  = joinpath(dir, "time_series")
    isdir(frames) || mkpath(frames)
    isdir(ts_dir) || mkpath(ts_dir)

    # === Visualization params ===
    UPSAMPLE_FACTOR = 5             # higher = finer detail
    FIG_RES   = (2400, 1200)        # higher pixel resolution
    PX_PER_UNIT = 4
    VIDEO_FPS = 30
    VIDEO_CRF = 18
    WAVE_ALPHA = 1.0

    # === Model grid sizes ===
    nx_full, nz_full = 438, 261
    nx_crop, nz_crop = 337, 168
    n_expected = nx_full * nz_full

    # === Physical coordinates (km) ===
    x_min, x_max = 0.0, 100.0
    z_min, z_max = -8.0, 43.70
    x_full = range(x_min, x_max; length=nx_full)
    z_full = range(z_min, z_max; length=nz_full)

    # === Cropping indices (ix first) ===
    ix_start = Int(floor((nx_full - nx_crop) / 2)) + 1
    ix_end   = ix_start + nx_crop - 1
    iz_start = Int(floor((nz_full - nz_crop) / 2)) + 1
    iz_end   = iz_start + nz_crop - 1

    # === Cropped axes ===
    x = x_full[ix_start:ix_end]
    z_data_min = -3.776
    z = range(z_data_min, z_max; length=nz_crop)

    # === Files ===
    files = sort(filter(f -> endswith(f, ".bin"), readdir(dir; join=true)))
    @assert !isempty(files) "No .bin files found in $dir"

    # === Source & receivers ===
    ix_src, iz_src = 125, 65
    x_src, z_src = x[ix_src], z[iz_src]
    println("Source indices (ix, iz) = ($ix_src, $iz_src)")
    println("Source coords (x, z) = ($(round(x_src, digits=2)), $(round(z_src, digits=2))) km")

    receiver_coords = [(22.0, 2.0), (22.0, 5.0), (21.90, 3.0), (21.90, 1.5), (21.90, 4.0)]
    receiver_indices = [(argmin(abs.(x .- xr)), argmin(abs.(z .- zr))) for (xr, zr) in receiver_coords]
    time_series = [Float32[] for _ in receiver_indices]

    # === Global amplitude range ===
    global_max = 0.0
    for file in files
        raw = read(file)
        data = reinterpret(Float32, raw)
        if length(data) == n_expected
            data = reshape(data, (nx_full, nz_full))
            data = data[ix_start:ix_end, iz_start:iz_end]
            c = 5e-4
            B = sign.(data) .* log1p.(abs.(data) ./ c)
            global_max = max(global_max, maximum(abs, B))
        end
    end
    println("Global color scale range: ±$(round(global_max, digits=4))")

    my_wave_colormap = cgrad([RGB(0,0,0.8), RGB(1,1,1), RGB(0.9,0,0)], scale=:linear)

    # === Loop ===
    for (i, file) in enumerate(files)
        raw = read(file)
        data = reinterpret(Float32, raw)
        if length(data) != n_expected
            @warn "Skipping $(basename(file)): length $(length(data)) ≠ $n_expected"
            continue
        end

        data = reshape(data, (nx_full, nz_full))
        data = data[ix_start:ix_end, iz_start:iz_end]

        c = 5e-4
        B = sign.(data) .* log1p.(abs.(data) ./ c)

        # upsample
        B_plot, xf, zf = upsample_field(B, x, z; factor=UPSAMPLE_FACTOR)

        # figure
        fig = Figure(resolution=FIG_RES, fontsize=16)
        rowgap!(fig.layout, 0)
        ax = Axis(fig[1, 1],
                  xlabel="x (km)", ylabel="z (km)", ylabelsize=18, xlabelsize=18,
                  xgridvisible=true, ygridvisible=true,
                  xgridcolor=(:gray, 0.3), ygridcolor=(:gray, 0.3),
                  spinewidth=1.4)

        xlims!(ax, first(x), last(x))
        ylims!(ax, z_max, z_min)

        custom_ticks = sort(unique(vcat(collect(-10:10:50), [z_data_min])))
        ax.yticks = (custom_ticks, string.(round.(custom_ticks, digits=3)))

        hm = GLMakie.heatmap!(ax, xf, zf, B_plot;
            colormap=my_wave_colormap,
            colorrange=(-global_max, global_max),
            alpha=WAVE_ALPHA, interpolate=false)  # <- finer texture

        # Source & receivers
        GLMakie.scatter!(ax, [x_src], [z_src];
            color=:red, marker=:star5, markersize=14, strokewidth=1.5, label="Source")

        for (ri, (xr, zr)) in enumerate(receiver_coords)
            GLMakie.scatter!(ax, [xr], [zr];
                color=:lime, marker=:utriangle, rotation=π,
                markersize=16, strokewidth=1.2, label="Receiver $(ri)")
            GLMakie.text!(ax, xr + 0.5, zr, text=string(ri),
                align=(:left, :center), color=:black, fontsize=15)
        end

        GLMakie.Colorbar(fig[1, 2], hm, label="Amplitude", labelsize=16,
            colorrange=(-global_max, global_max))
        axislegend(ax; position=:rb, orientation=:vertical, patchsize=(10,10), labelsize=14)

        out_png = joinpath(frames, @sprintf("frame_%05d.png", i))
        GLMakie.save(out_png, fig; px_per_unit=PX_PER_UNIT)
        println("Saved high-res snapshot: $out_png")

        # receiver samples
        for (ri, (ix_r, iz_r)) in enumerate(receiver_indices)
            push!(time_series[ri], data[ix_r, iz_r])
        end
    end

    println(" All snapshots saved in: $frames")

    # === Save time series ===
    Δt = 0.005
    t = collect(0:Δt:(length(files)-1)*Δt)
    for (ri, (xr, zr)) in enumerate(receiver_coords)
        series = time_series[ri]
        out_ascii = joinpath(ts_dir, @sprintf("receiver_%02d_x%.1f_z%.1f.txt", ri, xr, zr))
        writedlm(out_ascii, hcat(t, series))
        fig = Figure(resolution=(900, 360))
        ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Amplitude",
                  title=@sprintf("Receiver %d (x=%.1f km, z=%.1f km)", ri, xr, zr))
        GLMakie.lines!(ax, t, series; linewidth=1.6)
        GLMakie.save(replace(out_ascii, ".txt" => ".png"), fig; px_per_unit=2)
    end

    # === Make video ===
    output_video = joinpath(dir, "wavefield_highres.mp4")
    cd(frames) do
        cmd = `ffmpeg -framerate $VIDEO_FPS -i frame_%05d.png \
            -c:v libx264 -preset slow -crf $VIDEO_CRF -pix_fmt yuv420p \
            -vf "fps=$VIDEO_FPS,scale=trunc(iw/2)*2:trunc(ih/2)*2" $output_video`
        run(cmd)
    end
    println(" High-resolution video saved to: $output_video")
end

main()
