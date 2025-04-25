using CairoMakie, JLD2,Symbolics


cases=[]
cases = push!(cases,(name="λᵤ",u=cos(x),β=sin(x)+2))
cases = push!(cases,(name="2λᵤ",u=cos(x),β=sin(x/2) + 2))
cases = push!(cases,(name="λᵤ, shifted",u=cos(x),β=sin(x+π/3) + 2))
cases = push!(cases,(name="λᵤ/2",u=cos(x),β=cos(x).^2 + 1))
cases = push!(cases,(name="quandratic",u=cos(x),β=x^2+ 1))
cases = push!(cases,(name="homo",u=cos(x),β=1.0))
#

L = 10.0*π 
fig=Figure()
ax=Axis(fig[1,1];title="κ model")
for iCase in eachindex(cases)
    @unpack name,u,β = cases[iCase]
    Ncases=length(cases)+1
    Nx=10001
    Δx=L/Nx
    X=[(i-1)*Δx for i ∈ range(1,Nx)]
    κ=[Symbolics.value(substitute(β,Dict(x=>X[i]))) for i ∈ range(1,Nx)]
    colors = [get(Makie.colorschemes[:viridis], (i - 1) / (Ncases - 1)) for i in 1:Ncases]    
    lines!(ax,X,κ,color=colors[iCase],label=cases[iCase].name)
end

display(fig)

fig =Figure()
ax=Axis(fig[1,1]; title="Misfit")




@load "tmp_misfit.jld2" misfit




fig =Figure()
ax=Axis(fig[1,1]; title="Misfit")
N=length(cases)+1
colors = [get(Makie.colorschemes[:viridis], (i - 1) / (N - 1)) for i in 1:N]
for iCase in eachindex(cases)
    scatter!(ax,logsOfHinverse,log.(misfit[:,iCase,1]),color=colors[iCase],label=cases[iCase].name)
end


O_2=.-2.0*logsOfHinverse
O_4=.-4.0*logsOfHinverse
#O_8=log.(misfit[1,1,1]).-1.0*logsOfHinverse
#lines!(ax,logsOfHinverse,O_1,color=:black,label="O1")
lines!(ax,logsOfHinverse,O_2,color=:black,label="O2")
lines!(ax,logsOfHinverse,O_4,color=:black,label="O4")
#lines!(ax,logsOfHinverse,O_8,color=:black,label="O8")
#axislegend(ax,position=:lb)
display(fig)