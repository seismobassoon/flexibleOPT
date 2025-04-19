include("../src/batchDiffTools.jl")

# I kind of use what Tibo did in ExactSolutions.jl

function timeMarchingScheme(opt, Nt, Δnum;sourceType="Ricker",t₀=50,f₀=0.03,initialCondition=0.0)

    #region reading data and source time functions

    operators = opt["numOperators"]
    costfunctions,fieldLHS,fieldRHS,champsLimité = operators
    #@show size(costfunctions),size(fieldLHS[1,1]),size(fieldRHS)
    @show champsLimité[1,1]

    NField = size(costfunctions)[1]
    NCostfunctions=size(costfunctions)[2]
    if size(fieldLHS)[1]
    end

    itVec=collect(1:1:Nt)
    t=(itVec.-1).*Δnum[end] # time vector # if it's not time marching t will give you just 0.0 regardless of Δnum[end]

    #@show costfunctions[1,431] # check if we can see the source terms

    #specificication of parameters such as Nt (which can be 1 for no time-marching scheme) and source time function 


    if sourceType === "Ricker"
        # making a small source

        # if you want to plot:
        # myRicker(x)=Ricker(x,20,0.03)
        # scatter(0:1:100,myRicker) or lines 

       
        myRicker(x)=Ricker(x,t₀,f₀)
        
        sourceTime = myRicker.(t)
    end

    #endregion

    #region initial condition

    tmpField = similar(fieldLHS)
    tmpSource = similar(fieldRHS)

    


    #endregion

    
    #region

    for it in itVec

    end


    #endregion



end

#=
function dummy()
# Time loop
for it=1:nt
    U00 .= U0
    U0  .= U
    t   += Δt 
    noisy==true ? @printf("########### Step %06d ###########\n", it) : nothing

    for i in eachindex(f)
        sol   = Analytics(problem, @SVector([xc[i]; t-Δt])) # ???
        fC    = sol.s 
        sol   = Analytics(problem, @SVector([xc[i]-Δx; t-Δt]))
        fW    = sol.s 
        sol   = Analytics(problem, @SVector([xc[i]+Δx; t-Δt]))
        fE    = sol.s 
        sol   = Analytics(problem, @SVector([xc[i]; t-2*Δt]))
        fS    = sol.s 
        sol   = Analytics(problem, @SVector([xc[i]; t-0*Δt]))
        fN    = sol.s 
        f[i]  = fC - 1*( 4*fC - fW - fE - fS - fN )/12

    end

    r1 = 1.0
    for iter=1:10

        # Residual evaluation: T is found if F = 0
        Residual!(F, U, U0, U00, f, G, Gc, ρ, Δx, Δt, x, t)
        Res_closed! = (F, U) -> Residual!(F, U, U0, U00, f, G, Gc, ρ, Δx, Δt, x, t)
        r = norm(F)/ncx
        if iter==1 r1 = r; end
        noisy==true ? @printf("## Iteration %06d: r/r1 = %1.2e ##\n", iter, r/r1) : nothing
        if r/r1 < 1e-8 break end
            
        # Jacobian assembly
        forwarddiff_color_jacobian!(J, Res_closed!, U, colorvec = colors)

        # Solve
        δU   .= .-J\F

        # update
        U    .+= δU
    end
end
end
=#