using CairoMakie, Symbolics

include("../src/batchNewSymbolics.jl")

# this is a home-made B-spline functions' plot in order to understand the validity of truncation at left and right extremeties

@variables x
∂x = Differential(x)

# plot parameters

pointPerSegment = 20

# input parameters

# the left and right indices
νₗ = 0
νᵣ = 10

# δy for plotting, Δy for the real discretisation
Δy = 1.0

# maximum order of B-spline

maximumOrder = 5 + 1



# let's start

numberNodes = νᵣ - νₗ + 1 # [νₗ,νₗ + 1), ⋯, [νᵣ - 1 to νᵣ), {νᵣ}

nodeIndices = collect(νₗ:1:νᵣ)

nodes = Δy .* nodeIndices
append!(nodes, nodes[end]) # so that the last segment be {νᵣ}

# colour bar

colors = [get(Makie.colorschemes[:gist_rainbow], (i - 1) / (numberNodes - 1)) for i in 1:numberNodes]

# b-splines

b = zeros(Num, numberNodes, numberNodes, maximumOrder)



for ι in 0:1:maximumOrder-1
    local fig = Figure()
    local ax = Axis(fig[1, 1]; title="B-spline function of $ι-th order")

    if ι === 0
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
        
            b[tmpν, tmpν, ι+1] = 1
        end
    else

        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1

            local neighbour = Int((1 - (-1)^ι) / 2)
            floor = Int((ι - neighbour) / 2)
            ceiling = Int((ι + neighbour) / 2)

            # the denominator for the ι>0

            denominator = ι * Δy

            # for the upgoing part

            if ν - neighbour <= νᵣ && ν - neighbour >= νₗ
                @show rightlimit = minimum((numberNodes, tmpν + floor))
                @show leftlimit = maximum((1, tmpν - ceiling))
                for νSegment in leftlimit:1:rightlimit
                    tmpνSegment = νSegment - νₗ + 1
                    numerator = x - (ν - ceiling) * Δy
                    b[tmpνSegment, tmpν, ι+1] += numerator / denominator * b[tmpνSegment, tmpν-neighbour, ι]
                end
            end

            #for the downgoing part

            if ν - neighbour + 1 <= νᵣ && ν - neighbour + 1 >= νₗ
                rightlimit = minimum((numberNodes, tmpν + floor + 1))
                leftlimit = maximum((1, tmpν - ceiling + 1))
                for νSegment in leftlimit:1:rightlimit
                    tmpνSegment = νSegment - νₗ + 1
                    numerator = (ν + floor + 1) * Δy - x
                    b[tmpνSegment, tmpν, ι+1] += numerator / denominator * b[tmpνSegment, tmpν-neighbour+1, ι]
                end

            end
        end
    end

    for ν in nodeIndices # this will run for all the ν related to nodes
        tmpν = ν - νₗ + 1
        for νSegment in nodeIndices
            tmpνSegment = νSegment - νₗ + 1
            local xs = range(nodes[tmpνSegment], nodes[tmpνSegment+1], length=pointPerSegment)
            local ys = similar(xs)
            for ix in eachindex(xs)
                @show b[tmpνSegment, tmpν, ι+1],xs[ix]
                ys[ix] = Symbolics.eval(b[tmpνSegment, tmpν, ι+1], Dict(x => xs[ix]))
            end
            lines!(ax, xs, ys, color=color = colors[tmpν])
        end
    end
    #axislegend(ax)

    display(fig)
end
