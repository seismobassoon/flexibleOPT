dir="C:/Users/user/Desktop/stage 2A/données/MantleConvectionTakashi/op_first_run/"
rhoFiles=myListDir(dir; pattern=r"test_rho\d");

#function to plot (rho(r)-rhoprem(r)/rhoprem(r))
iTime = 100
file1=rhoFiles[iTime]
field1, Xnode, Ynode, rcmb = readStagYYFiles(file1)
arrayRadius = sqrt.(Xnode.^2 .+ Ynode.^2)


premDensities = myDensityFrom1DModel.(arrayRadius)
newpremDensities = premDensities
frho = ifelse.(newpremDensities .== 0.0, 0.0, (field1 .* 1e-3 .- newpremDensities) ./ newpremDensities)


fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),frho,correlationLength,epsilon2);
fig2 = Figure()
ax2 = Axis(fig2[1,1],aspect = 1)
hm2=heatmap!(ax2,fi3,colormap=cgrad(:viridis),colorrange=(-0.005,0.005))
Colorbar(fig2[:, 2], hm2)

display(fig2)


iTime=200

function myPlot2DConvectionModel(fieldname, filename)

    file = filename[iTime]
    field, Xnode, Ynode, rcmb = readStagYYFiles(file)
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)


    if fieldname ==="temperature"
        colormap=cgrad(:seismic)

    elseif fieldname === "rho" || fieldname === "composition"
        colormap=cgrad(:viridis)

    elseif fieldname === "wtr"
        colormap=cgrad(:blues)

    else 
        colormap=cgrad(:viridis)

    end
    
    hm=heatmap!(ax, fi, colormap=colormap, colorrange=(-0.005,0.005))
    Colorbar(fig[:,2], hm)
    display(fig)

end

myPlot2DConvectionModel("rho", rhoFiles)
