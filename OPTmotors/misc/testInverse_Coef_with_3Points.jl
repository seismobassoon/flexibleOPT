using  Pkg,  Symbolics, UnPack,Symbolics

cd(@__DIR__)
Pkg.activate("../..")


#include("../src/OPTwrappers.jl") 
include("../src/batchnewSymbolics.jl")

@variables dx
matrixA=[1 -dx dx^2//2; 1 0 0 ; 1 dx dx^2//2]
@show invA=inv(matrixA)


matrixleftA=[1 -2*dx 4*dx^2//2; 1 -dx dx^2//2 ; 1 0 0]
@show invLeftA=inv(matrixleftA)



staggeredleftA = [1 -(dx/2) (dx/2)^2//2; 1 dx/2 (dx/2)^2//2 ; 1 (3*dx/2) (3*dx/2)^2//2 ]
@show invStaggeredLeftA=inv(staggeredleftA)


matrixleftAreduced=[1 -dx dx^2//2 -dx^3//6 dx^4//12; 1 0 0 0 0]
invA=mySimplify(inv((matrixleftAreduced*transpose(matrixleftAreduced))))
@show a=mySimplify(transpose(matrixleftAreduced)*invA)