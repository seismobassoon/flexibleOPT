using GLMakie, Statistics, Printf, DelimitedFiles, FileIO, Colors

# read 2D binary file in Fortran order 
function read_binary_matrix(path, nz, nx)
    vec = Array{Float32}(undef, nz * nx)
    open(path, "r") do io
        read!(io, vec)
    end
    return reshape(vec, nx, nz)'  # reshape and transpose to get (nz, nx)
end

function main()
    
    dir = "/Users/hessiemohammadi/Documents/github/OptimallyAccurate2D/volcano_simu1/snapshots/bin_files1"
    datadir = "/Users/hessiemohammadi/Documents/FUJI/Github/flexibleDSM/OPTmotors/usefulScriptsForFriends"

    # Grid sizes
    nx_full, nz_full = 261, 438        # original full grid including margins
    nx_crop, nz_crop = 161, 338        # cropped central region
    n_expected = nx_full * nz_full

    # Physical coordinates (in km)
    x_min, x_max = -50.0, 50.0
    z_min, z_max = -10.0, 50.0
    x_full = range(x_min, x_max; length=nx_full)
    z_full = range(z_min, z_max; length=nz_full)

    # Cropping indices for center region 
    ix_start = Int(floor((nx_full - nx_crop) / 2)) + 1
    ix_end   = ix_start + nx_crop - 1
    iz_start = Int(floor((nz_full - nz_crop) / 2)) + 1
    iz_end   = iz_start + nz_crop - 1
    #println("Cropping indices: x[$ix_start:$ix_end], z[$iz_start:$iz_end]")

    # Cropped coordinate arrays
    x = x_full[ix_start:ix_end]

    # Data starts at -3.775 km, but figure still goes from -10 km
    z_data_min = -3.775
    z = range(z_data_min, z_max; length=nz_crop)

    # Load background velocity model for overlay
    vp_file = joinpath(datadir, "volcano_rock_topo_air_vp300.vp")
    println("Loading velocity model: $vp_file")
    vp = read_binary_matrix(vp_file, nz_crop, nx_crop)
    vp_min, vp_max = extrema(vp)
    #println("Vp range: $(round(vp_min, digits=2)) – $(round(vp_max, digits=2)) km/s")

    # colormap for velocity background 
    my_colormap = cgrad([
        RGB(1.0, 1.0, 1.0),    # Air (white)
        RGB(1.0, 0.0, 0.0),    # Red = low velocity
        RGB(1.0, 0.65, 0.0),   # Orange = medium velocity
        RGB(0.55, 0.27, 0.07)  # Brown = high velocity
    ], [0.0, 0.33, 0.66, 1.0])

    # Find all .bin files 
    files = sort(filter(f -> endswith(f, ".bin"), readdir(dir; join=true)))
    #println(" $(length(files)) binary files")

    # Detect source from first file
    first_file = first(files)
    

    raw_first = read(first_file)
    data_first = reinterpret(Float32, raw_first)
    if length(data_first) != n_expected
        error("Unexpected length $(length(data_first)) vs $n_expected")
    end

    field_first = reshape(data_first, (nx_full, nz_full))'
    field_first = field_first[iz_start:iz_end, ix_start:ix_end]  # crop center

    max_val = maximum(abs.(field_first))
    src_index = findfirst(x -> abs(x) == max_val, field_first)
    iz_src, ix_src = Tuple(src_index)
    println("Source indices (ix, iz) = ($ix_src, $iz_src)")
    println(" Amplitude at source = $max_val")

    x_src = x[ix_src]
    z_src = z[iz_src]
    println(" Source coordinates (x, z) = ($(round(x_src, digits=2)), $(round(z_src, digits=2))) km")

    # Define Receivers (in km)
    receiver_coords = [
        (-10.0, 1.0),
        (-5.0,  1.0),
        (0.0,   1.0),
        (10.0,  1.0),
        (20.0,   1.0),
    ]
    receiver_indices = [(argmin(abs.(x .- xr)), argmin(abs.(z .- zr))) for (xr, zr) in receiver_coords]
    time_series = [Float32[] for _ in receiver_indices]

    #  Compute color scale from all snapshots
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

    # Loop through snapshots 
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

        c = 1e-4
        B = sign.(data) .* log1p.(abs.(data) ./ c)

        # Figure 
        fig = Figure(size=(1000, 550))
        rowgap!(fig.layout, 0)  # tighter layout

        ax = Axis(fig[1, 1],
            xlabel="x (km)",
            ylabel="z (km)",
            xgridvisible = false,
            ygridvisible = false,
            )
        #ax.aspect = DataAspect()
        hidespines!(ax)

        # Full figure frame from -10 to 50
        xlims!(ax, minimum(x), maximum(x))
        ylims!(ax, 50.0, -10.0)

        # Show -3.776 km label on z-axis
        custom_ticks = sort(unique(vcat(collect(-10:10:50), [-3.776])))
        ax.yticks = (custom_ticks, string.(round.(custom_ticks, digits=3)))

        # Fill air layer (white) above z = -3.775 
        poly!(
            ax,
            Point2f[
                (x_min, -10.0),
                (x_max, -10.0),
                (x_max, z_data_min),
                (x_min, z_data_min)
            ],
            color = :white
        )

        #  Velocity background 
        heatmap!(ax, x, z, vp;
            colormap = my_colormap,
            colorrange = (vp_min, vp_max),
            alpha = 0.6)

        
        #  Wavefield overlay 
        # Custom colormap for wavefield
        my_wave_colormap = cgrad([
            RGB(0.0, 0.0, 1.0),    # blue for negative
            RGB(1.0, 1.0, 1.0),    # white at zero
            RGB(1.0, 0.0, 0.0)     # red for positive
       ])

        hm = heatmap!(ax, x, z, B;
            colormap = my_wave_colormap,
            colorrange = (-global_max, global_max),
            alpha = 0.85)

        #  Source marker 
        scatter!(ax, [x_src], [z_src];
            color = :red, marker = :star5, markersize = 10, strokewidth = 1.0, label = "Source")

        #  Receivers 
        for (ri, (xr, zr)) in enumerate(receiver_coords)
            scatter!(ax, [xr], [zr];
                color = :green, marker = :utriangle, rotation = π,
                markersize = 12, strokewidth = 1, label = "Receiver $(ri)")
            text!(ax, xr + 0.5, zr, text = string(ri),
                  align = (:left, :center), color = :black, fontsize = 11)
        end

        #  Colorbars
        Colorbar(fig[1, 2], hm, label="Amplitude", colorrange=(-global_max, global_max))
        Colorbar(fig[1, 0], limits=(vp_min, vp_max), colormap=my_colormap, label="Vp (km/s)")
        axislegend(ax;
            position = :rb,
            orientation = :vertical,
            patchsize = (8, 8),      
            labelsize = 10,          
            rowgap = 0, colgap = 2,  
            padding = (2, 2, 2, 2))  


        # Save snapshot
        out_png = replace(file, ".bin" => ".png")
        save(out_png, fig)
        println("Saved snapshot to: $out_png")

        # Accumulate time series at receivers 
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

        fig = Figure(size=(700, 300))
        ax = Axis(fig[1, 1],
            xlabel="Time (s)", ylabel="Amplitude",
            title=@sprintf("Receiver %d  (x=%.1f km, z=%.1f km)", ri, xr, zr))
        lines!(ax, t, series; color=:black, linewidth=1.2)
        save(replace(out_ascii, ".txt" => ".png"), fig)
    end

    println("All receiver time series saved in $ts_dir")

    #  Create video
    output_video = joinpath(dir, "wavefield.mp4")
    cd(dir) do
        cmd = `ffmpeg -framerate 25 -pattern_type glob -i "*.png" -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" $output_video`
    
        run(cmd)
        println("Video saved to: $output_video")
    end
end

main()
