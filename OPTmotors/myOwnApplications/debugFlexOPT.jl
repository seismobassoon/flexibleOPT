### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 6c0d05f8-b66f-11f0-84c5-f7e749a4268e
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "../.."))
    using Pluto
    #import Pkg
end

# ╔═╡ ec4d52e6-c9f1-46c4-8362-73d231fcb0e3
using BenchmarkTools

# ╔═╡ cc967e2f-235a-41b0-8da9-75bbcf3e3035
# below are the tools to debug the code
#using Revise # if we use Revise, myInclude will be Revise.myIncludet
using Profile, StatProfilerHTML

# ╔═╡ b13b3658-3d66-4b77-82c5-d250db494e96
include("../src/batchRevise.jl")

# ╔═╡ 6edc09ed-7b03-4429-b866-c3c8c0b70816
cd(@__DIR__)

# ╔═╡ c0ed6199-bd7c-4c2e-a5fe-545ec140a76a
dirOPTMotors=@__DIR__

# ╔═╡ 35b5cb0f-b775-470a-a285-7cb4b4915369
@show joinpath(dirOPTMotors,"../src/batchRevise.jl")

# ╔═╡ f762ef8c-1fcf-4978-9169-746292aaab2f

myInclude("../src/batchUseful.jl")

# ╔═╡ 7073cea1-7df3-434d-8480-1dbe500f6a6b
myInclude(joinpath(dirOPTMotors,"../src/imageReader.jl")) # read 2D images for models

# ╔═╡ c1918676-e17e-4b36-8489-da11770d1e0a
myInclude(joinpath(@__DIR__,"../src/OPTwrappers.jl"))

# ╔═╡ 81a17cd5-6f9b-4260-bb48-0cabf640cff4
myInclude("../src/OPTnewEngines.jl")  # I do this to test

# ╔═╡ Cell order:
# ╠═6c0d05f8-b66f-11f0-84c5-f7e749a4268e
# ╠═6edc09ed-7b03-4429-b866-c3c8c0b70816
# ╠═c0ed6199-bd7c-4c2e-a5fe-545ec140a76a
# ╠═ec4d52e6-c9f1-46c4-8362-73d231fcb0e3
# ╠═cc967e2f-235a-41b0-8da9-75bbcf3e3035
# ╠═35b5cb0f-b775-470a-a285-7cb4b4915369
# ╠═b13b3658-3d66-4b77-82c5-d250db494e96
# ╠═f762ef8c-1fcf-4978-9169-746292aaab2f
# ╠═7073cea1-7df3-434d-8480-1dbe500f6a6b
# ╠═c1918676-e17e-4b36-8489-da11770d1e0a
# ╠═81a17cd5-6f9b-4260-bb48-0cabf640cff4
