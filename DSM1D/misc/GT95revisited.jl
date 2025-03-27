using Symbolics
using CairoMakie
@variables  Δz::Real 

subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))


N=50 # number of nodes
K=3 # maximum order of spline
Nresolution =10 # number of points between nodes
Xspline=Array{Symbolics.Real,3}(undef,K+1+1,N+1,N*Nresolution+1);
XsplinePrime=Array{Symbolics.Real,3}(undef,K+1+1,N+1,N*Nresolution+1);
u₀array=Array{Symbolics.Real,1}(undef,N+1);
∂u₀array=Array{Symbolics.Real,1}(undef,N+1);
∂²u₀array=Array{Symbolics.Real,1}(undef,N+1);
∂³u₀array=Array{Symbolics.Real,1}(undef,N+1);
∂⁴u₀array=Array{Symbolics.Real,1}(undef,N+1);
uarray=Array{Symbolics.Real,2}(undef,N+1,N*Nresolution+1);
nodes=[Δz*(i-1) for i in 1:N+1]

xarray=[Symbolics.simplify(Δz*(i-1)/Nresolution) for i in 1:N*Nresolution+1]


function Spline(K,N,nodes,x)
    for k = 0:K # k is shifted after when stored (but not now in this routine) with respect to the ordering of spline functions
        for j = 1:N+1 # j is shifted by 1 compared to the paper
          for xi= 1:N*Nresolution+1 # xi is shifted by 1 
            x=xarray[xi]
            X_local=0
            X_local_prime=0
            if k>0 
                if j+k <= N+1 
                    
                    X_local += (x-nodes[j])//(nodes[j+k]-nodes[j])*Xspline[k-1+1,j,xi]
                    X_local_prime +=  (x-nodes[j])/(nodes[j+k]-nodes[j])*XsplinePrime[k-1+1,j,xi] + Xspline[k-1+1,j,xi]/(nodes[j+k]-nodes[j])
                end
                if j+k+1 <= N+1
                    X_local +=  (nodes[j+k+1]-x)/(nodes[j+k+1]-nodes[j+1])*Xspline[k-1+1,j+1,xi]
                    X_local_prime += -Xspline[k-1+1,j+1,xi]/(nodes[j+k+1]-nodes[j+1]) + (nodes[j+k+1]-x)/(nodes[j+k+1]-nodes[j+1])*XsplinePrime[k-1+1,j+1,xi]
                end      
            elseif k==0
                if j>=1 && j+1 <= N+1 && nodes[j]/Δz <= x/Δz < nodes[j+1]/Δz
                    X_local = 1
                end
            end
            Xspline[k+1,j,xi]=Symbolics.simplify(X_local)
            XsplinePrime[k+1,j,xi]=Symbolics.simplify(X_local_prime)
          end
        end
    end
    

for j = 1:N+1
    u₀array[j]="u₀"*subscript(j)
    ∂u₀array[j]="∂u₀"*subscript(j)
    ∂²u₀array[j]="∂²u₀"*subscript(j)
    ∂³u₀array[j]="∂³u₀"*subscript(j)
    ∂⁴u₀array[j]="∂⁴u₀"*subscript(j)
end

for j = 1:N+1
    for xi= 1:N*Nresolution+1
        x=xarray[xi]
        uarray[j,xi]=u₀array[j]+(x-nodes[j])*∂u₀∂zarray[j]+1//2*(x-nodes[j])^2*∂²u₀∂z²array[j]+1//6*(x-nodes[j])^3*∂³u₀∂z³array[j]+1//24*(x-nodes[j])^4*∂⁴u₀∂z⁴array[j]
    end
end

for k = 0:K # k is shifted after when stored (but not now in this routine) with respect to the ordering of spline functions
    for j = 1:N+1 # j is shifted by 1 compared to the paper
      for xi= 1:N*Nresolution+1 # xi is shifted by 1 
        x=xarray[xi]
        X_local=0
        X_local_prime=0
        if k>0 
            if j+k <= N+1 
                
                X_local += (x-nodes[j])//(nodes[j+k]-nodes[j])*Xspline[k-1+1,j,xi]
                X_local_prime +=  (x-nodes[j])/(nodes[j+k]-nodes[j])*XsplinePrime[k-1+1,j,xi] + Xspline[k-1+1,j,xi]/(nodes[j+k]-nodes[j])
            end
            if j+k+1 <= N+1
                X_local +=  (nodes[j+k+1]-x)/(nodes[j+k+1]-nodes[j+1])*Xspline[k-1+1,j+1,xi]
                X_local_prime += -Xspline[k-1+1,j+1,xi]/(nodes[j+k+1]-nodes[j+1]) + (nodes[j+k+1]-x)/(nodes[j+k+1]-nodes[j+1])*XsplinePrime[k-1+1,j+1,xi]
            end      
        elseif k==0
            if j>=1 && j+1 <= N+1 && nodes[j]/Δz <= x/Δz < nodes[j+1]/Δz
                X_local = 1
            end
        end
        Xspline[k+1,j,xi]=Symbolics.simplify(X_local)
        XsplinePrime[k+1,j,xi]=Symbolics.simplify(X_local_prime)
      end
    end
end



#XsplinePrime=Symbolics.substitute(XsplinePrime,Dict([Δz=>1]))
#Xspline=Symbolics.substitute(Xspline,Dict([Δz=>1]))
println(XsplinePrime[1,3,:])
println(Xspline[1,3,:])

println(XsplinePrime[2,3,:])
println(Xspline[2,3,:])

println(XsplinePrime[3,3,:])
println(Xspline[3,3,:])