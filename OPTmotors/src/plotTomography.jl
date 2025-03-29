using CairoMakie

function myPcolor(A)
    # This function reads 2D array and plot with the aid of CairoMakie

    if length(size(A)) != 2
        @error "The array is not 2D"
    end
    
    nx=size(mymatrix)[1]
    ny=size(mymatrix)[2]

    x=1:1:nx
    y=1:1:ny

    myValue(x,y) = A[x,y]

    fig=Figure()

    ax = Axis(fig[1, 1], title = "Matrix elements", xlabel = "X", ylabel = "Y")
    hm = heatmap!(ax,x,y,value,colormap=Reverse(:deep))
    Colorbar(fig[1,2],hm,label="color scale")
    return fig
end
