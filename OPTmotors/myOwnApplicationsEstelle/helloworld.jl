

dir="C:/Users/user/Desktop/stage 2A/données/test_estelle2/op_first_run"
rhoFiles=myListDir(dir; pattern=r"test_rho\d");

#function to plot (rho(r)-rhoprem(r)/rhoprem(r))

for file3 in rhoFiles[1:20]
#file3 = rhoFiles[iTime]
    local field3, Xnode, Ynode, rcmb = readStagYYFiles(file3)
    local arrayRadius = sqrt.(Xnode.^2 .+ Ynode.^2)


    function my_PREM(arrayRadius)
        density  = DSM1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(DSM1D.my1DDSMmodel, arrayRadius*1.e-3)
        return density
    end

    local premDensities = my_PREM.(arrayRadius)

#enlever NaN et Infs
    epsilon = 1e-6 
    newpremDensities = []
    for x in premDensities
        if !isfinite(x) || abs(x) < epsilon
            push!(newpremDensities, epsilon)
        else
            push!(newpremDensities,x)
        end
    end

    local frho = (field3*1.e-3 .- newpremDensities) ./ newpremDensities
    local fi3,s = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),frho,correlationLength,epsilon2);



#mat_frho = reshape(frho, 33540,1) #pour avoir une matrice dans heatmap!
#mat_frho_itp = interpolate(mat_frho,(BSpline(Linear()), NoInterp())) #car vecteur a 1 colonne
#mat_frho_itp_resam = Resampler(mat_frho) #trop de données donc resampler

    local fig2 = Figure()
    local ax2 = Axis(fig2[1,1],aspect = 1)
    local hm2=heatmap!(ax2,fi3,colormap=cgrad(:viridis))
    Colorbar(fig2[:, 2], hm2)
    display(fig2)
end 
