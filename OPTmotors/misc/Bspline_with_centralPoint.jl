using CairoMakie


# this is a home-made B-spline functions' plot in order to understand the validity of truncation at left and right extremeties

# input parameters

# the left and right indices
νₗ = 0
νᵣ =10

# δy for plotting, Δy for the real discretisation

δy = 0.05
Δy = 1.0

# maximum order of B-spline

maximumOrder=5+1



# let's start

ratioΔδ = Int(Δy/δy) 

if ratioΔδ*δy !== Δy  
    @error "this ratio should be an integer!" 
end


numberFunctions=νᵣ-νₗ+1
numberNodes = νᵣ - νₗ + 1

nodeIndices = collect(νₗ:1:νᵣ)
nodes =  Δy .* nodeIndices

discretisedIndices = collect(νₗ*ratioΔδ:1:νᵣ*ratioΔδ)
discretisedY = δy .* discretisedIndices
NdiscretisedY = length(discretisedY) 

discretisedNodeIndices = nodeIndices .* ratioΔδ


b= zeros(Float64,NdiscretisedY,numberFunctions,maximumOrder)


                  
            
# starting from ι = 0

ι = 0


fig =Figure()
ax=Axis(fig[1,1]; title="B-spline function of $ι-th order")
       

for ν in nodeIndices # this will run for all the ν related to nodes
    tmpν = ν - νₗ + 1
    if ν+1 <= νᵣ
        for tmpIndex in discretisedNodeIndices[tmpν]-discretisedNodeIndices[1]+1:1:discretisedNodeIndices[tmpν+1]-discretisedNodeIndices[1]
            b[tmpIndex,tmpν,ι+1] = 1.0
        end
    else
        tmpIndex = discretisedNodeIndices[tmpν]-discretisedNodeIndices[1]+1
        b[tmpIndex,tmpν,ι+1] = 1.0
    end
    scatter!(ax,discretisedY,b[:,tmpν,ι+1],label="$ν-th function")  
end

axislegend(ax)
display(fig)


# for ι > 0

for ι in 1:1:maximumOrder-1
    local fig =Figure()
    local ax=Axis(fig[1,1]; title="B-spline function of $ι-th order")
        
    for ν in nodeIndices # this will run for all the ν related to nodes
        tmpν = ν - νₗ + 1
        
        local neighbour= Int((1-(-1)^ι)/2)
        floor = Int((ι-neighbour)/2)
        ceiling = Int((ι+neighbour)/2)

        if ν-neighbour <= νᵣ && ν-neighbour >= νₗ
            denominator = Δy*ι
            rightlimit = discretisedNodeIndices[minimum((numberNodes,tmpν+floor))]-discretisedNodeIndices[1]+1
            leftlimit = discretisedNodeIndices[maximum((1,tmpν-ceiling))]-discretisedNodeIndices[1]+1

            for tmpIndex in leftlimit:1:rightlimit
                numerator = δy * (tmpIndex - ((ν-ceiling)*ratioΔδ-discretisedNodeIndices[1]+1))
                #@show ν, tmpν, neighbour
                b[tmpIndex,tmpν,ι+1] += numerator/denominator*b[tmpIndex,tmpν-neighbour,ι]
            end
        end
        
        if ν-neighbour+1 <= νᵣ && ν-neighbour+1 >= νₗ
            denominator = Δy*ι
            rightlimit = discretisedNodeIndices[minimum((numberNodes,tmpν+floor+1))]-discretisedNodeIndices[1]+1
            leftlimit = discretisedNodeIndices[maximum((1,tmpν-ceiling+1))]-discretisedNodeIndices[1]+1
            for tmpIndex in leftlimit:1:rightlimit
                @show numerator = δy * (((ν+floor+1)*ratioΔδ-discretisedNodeIndices[1]+1)-tmpIndex)
                b[tmpIndex,tmpν,ι+1] += numerator/denominator*b[tmpIndex,tmpν-neighbour+1,ι]
            end

        end

        scatter!(ax,discretisedY,b[:,tmpν,ι+1],label="$ν-th function")  
    end
    #axislegend(ax)
    
    display(fig)
end
