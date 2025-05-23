# this programme should be well saved somewhere 

using CairoMakie, Symbolics,Pkg
cd(@__DIR__)
Pkg.activate("../..")
include("../src/batchNewSymbolics.jl")

# this is a home-made B-spline functions' plot in order to understand the validity of truncation at left and right extremeties

function makeBspline()

@variables x Δx ξ
∂x = Differential(x)

# plot parameters

pointPerSegment = 20

# input parameters

# the left and right indices
νₗ = 0
νᵣ = 20

# δy for plotting, Δy for the real discretisation
Δy = 1.0

# maximum order of B-spline

maximumOrder = 3 + 1



# let's start

numberNodes = νᵣ - νₗ + 1 # [νₗ,νₗ + 1), ⋯, [νᵣ - 1 to νᵣ), {νᵣ}

nodeIndices = collect(νₗ:1:νᵣ)

nodesSymbolic = Δx .* nodeIndices
append!(nodesSymbolic,nodesSymbolic[end])

nodes = Δy .* nodeIndices
append!(nodes, nodes[end]) # so that the last segment be {νᵣ}


# colour bar

colors = [get(Makie.colorschemes[:Paired_5], (i - 1) / (numberNodes - 1)) for i in 1:numberNodes]

# b-splines

b = zeros(Num, numberNodes, numberNodes, maximumOrder)
bX = zeros(Num,numberNodes, 2, maximumOrder)

# b-splines derivatives

b_deriv = zeros(Num,numberNodes, numberNodes,maximumOrder,maximumOrder)
local figBig = Figure(size = (500*maximumOrder, 500*maximumOrder))

for ι in 0:1:maximumOrder-1
    
    local neighbour = Int((1 - (-1)^ι) / 2)
    floor = Int((ι - neighbour) / 2)
    ceiling = Int((ι + neighbour) / 2)
    if ι === 0
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
        
            b[tmpν, tmpν, ι+1] = 1
        end
    else

        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
            # the denominator for the ι>0

            denominator = BigInt(ι) * Δx

            # for the upgoing part

            if ν - neighbour <= νᵣ && ν - neighbour >= νₗ
                rightlimit = minimum((numberNodes, tmpν + floor))
                leftlimit = maximum((1, tmpν - ceiling))
                for νSegment in leftlimit:1:rightlimit
                    tmpνSegment = νSegment #- νₗ + 1
                    numerator = x - BigInt(ν - ceiling) * Δx
                    b[tmpνSegment, tmpν, ι+1] += mySimplify(numerator / denominator * b[tmpνSegment, tmpν-neighbour, ι])
                end
            end

            #for the downgoing part

            if ν - neighbour + 1 <= νᵣ && ν - neighbour + 1 >= νₗ
                rightlimit = minimum((numberNodes, tmpν + floor + 1))
                leftlimit = maximum((1, tmpν - ceiling + 1))
                for νSegment in leftlimit:1:rightlimit
                    tmpνSegment = νSegment #- νₗ + 1
                    numerator = BigInt(ν + floor + 1) * Δx - x
                    b[tmpνSegment, tmpν, ι+1] += mySimplify(numerator / denominator * b[tmpνSegment, tmpν-neighbour+1, ι])
                end

            end

        end
    end

    # unfortunately the νₗ and νᵣ related functions should recover the 'lost' amplitudes due to the trunctation 
    # we need to re-build based on the residual

    for ν in nodeIndices
        tmpν = ν - νₗ + 1
        
        #b[:,tmpν,ι+1] .= 0
        if ν===νₗ  # this loop runs only for the extremeties!!!!!
            # the nodes involved at most are: νₗ, νₗ + 1, ⋯, νₗ + ceiling 
            bX[1:ceiling+1,1,ι+1] .=1
            for νSegment in 1:1:ceiling+1
                for νFunction in 2:1:numberNodes
                    bX[νSegment,1,ι+1] -= b[νSegment,νFunction,ι+1]
                end
            end
        elseif ν===νᵣ
            bX[numberNodes-ceiling-1:numberNodes,2,ι+1] .=1
            for νSegment in numberNodes-ceiling-1:1:numberNodes
                for νFunction in 1:1:numberNodes-1
                    bX[νSegment,2,ι+1] -= b[νSegment,νFunction,ι+1]
                end
            end
        end

    end

    b[:,1,ι+1]=bX[:,1,ι+1]
    b[:,end,ι+1]=bX[:,end,ι+1]

    
    # computing the derivatives 

    for i in 0:1:maximumOrder-1 
        if i === 0
            b_deriv[:,:,i+1,ι+1] = b[:,:,ι+1]
        else
            b_deriv[:,:,i+1,ι+1] = mySimplify.(∂x.(b_deriv[:,:,i,ι+1]))
        
        end

    end

   
    

    for i in 0:1:maximumOrder-1 
        local fig = Figure()
        local ax = Axis(fig[1, 1]; title="B-spline function of $ι-th order, $i-th derivative")
        local axBig=Axis(figBig[ι+1,i+1]; title="B-spline function of $ι-th order, $i-th derivative")
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
            for νSegment in nodeIndices
                tmpνSegment = νSegment - νₗ + 1
                local xs = range(nodes[tmpνSegment], nodes[tmpνSegment+1], length=pointPerSegment)
                local ys = similar(xs)
                # Precompute symbolic expression and turn into function
                expr = substitute(b_deriv[tmpνSegment, tmpν, i+1,ι+1],Dict(Δx=>Δy))
                #if ν === νₗ
                #    expr = bX[tmpνSegment,1,ι+1]
                #elseif ν === νᵣ
                #    expr = bX[tmpνSegment,2,ι+1]
                #end
                f_expr = Symbolics.build_function(expr, x; expression = false)
                local f = eval(f_expr)  # Now this works
                ys =f.(xs)
                
                lines!(ax, xs, ys, color= colors[tmpν])
                lines!(axBig,xs,ys,color=colors[tmpν])
            end
        end
        #axislegend(ax)

        display(fig)

    end

end
display(figBig)

b_deriv_ξ = similar(b_deriv)
for ι in 0:1:maximumOrder-1
    for i in 0:1:maximumOrder-1
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
            for νSegment in nodeIndices
                tmpνSegment = νSegment - νₗ + 1
                b_deriv_ξ[tmpνSegment,tmpν,i+1,ι+1]=substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>ξ+Δx*(tmpν-1)))
            end
        end
    end
end


integral_b= zeros(Num,numberNodes)
gvec = (@variables (g(x))[1:maximumOrder])[1]

for ι in 0:1:maximumOrder-1
    for i in 0:1:maximumOrder-1
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
            for νSegment in nodeIndices
                tmpνSegment = νSegment - νₗ + 1
                
                diff = Δx * (tmpνSegment - tmpν)
                tmpDic = Dict(x=>diff)
                integral_b[tmpν] -= (-1)^(i)*substitute(gvec[i+1],tmpDic)*substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>nodesSymbolic[tmpνSegment]))

                diff = Δx * (tmpνSegment - tmpν+1)
                tmpDic = Dict(x=>diff)
                integral_b[tmpν] += (-1)^(i)*substitute(gvec[i+1],tmpDic)*substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>nodesSymbolic[tmpνSegment+1]*Δx))
                
            end
        end
    end
end

#display.(b_deriv)
b_deriv_ξ=mySimplify.(b_deriv_ξ)
#b_deriv_extremes=mySimplify.(b_deriv_extremes)
integral_b=mySimplify.(integral_b)

# for special g(x) = C ξ^N 
@variables C N


dictionaryForSubstitute =Dict()


taylorNum = 1

for i in 1:1:maximumOrder
    @show taylorNum *= N+i
    dictionaryForSubstitute[gvec[i]]=x^(N+i)/taylorNum
    
end
@show dictionaryForSubstitute
integral_b_polys= zeros(Num,numberNodes,maximumOrder)
for ι in 0:1:maximumOrder-1
    for i in 0:1:maximumOrder-1
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
            for νSegment in nodeIndices
                tmpνSegment = νSegment - νₗ + 1
                
                tmp_b_deriv = b_deriv[tmpνSegment,tmpν,i+1,ι+1]
                tmpG = substitute(gvec[i+1],dictionaryForSubstitute)

                diff = nodesSymbolic[tmpνSegment]-nodesSymbolic[tmpν]
                tmpDic = Dict(x=>diff)
                if tmpνSegment !== tmpν
                    integral_b_polys[tmpν,ι+1] -= (-1)^(i)*substitute(tmpG,tmpDic)*substitute(tmp_b_deriv,Dict(x=>nodesSymbolic[tmpνSegment]))
                end

                diff = nodesSymbolic[tmpνSegment+1]-nodesSymbolic[tmpν]
                tmpDic = Dict(x=>diff)
                if tmpνSegment+1 !== tmpν
                    integral_b_polys[tmpν,ι+1] += (-1)^(i)*substitute(tmpG,tmpDic)*substitute(tmp_b_deriv,Dict(x=>nodesSymbolic[tmpνSegment+1]))
                end
                
            end
        end
    end
end
integral_b_polys=simplify.(integral_b_polys)

return b_deriv_ξ, b_deriv,integral_b,integral_b_polys
end


b_deriv_ξ, b_deriv,integral_b,integral_b_polys=makeBspline()
export b_deriv_ξ, b_deriv,integral_b,integral_b_polys