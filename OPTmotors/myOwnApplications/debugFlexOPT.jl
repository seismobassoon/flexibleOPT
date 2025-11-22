### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 6c0d05f8-b66f-11f0-84c5-f7e749a4268e
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "../.."))   # Path to your Project.toml
    Pkg.instantiate()
end

# ╔═╡ ec4d52e6-c9f1-46c4-8362-73d231fcb0e3
using BenchmarkTools

# ╔═╡ cc967e2f-235a-41b0-8da9-75bbcf3e3035
# below are the tools to debug the code
#using Revise # if we use Revise, myInclude will be Revise.myIncludet
using Profile, StatProfilerHTML

# ╔═╡ 9af87ac1-2331-4b9a-a352-24381cc0a7b7
include(joinpath(@__DIR__,"../src/batchImages.jl"))

# ╔═╡ 48d5c0c6-391c-4e38-b58d-00bf4b306694
include(joinpath(@__DIR__,"../src/imageReader.jl"))

# ╔═╡ 82d0aed5-5829-448a-9f97-dd3303b1d890
begin
include(joinpath(@__DIR__,"../src/batchNewSymbolics.jl"))
end

# ╔═╡ ccee7ed3-dc2f-4b97-bbab-4662120e57aa
include(joinpath(@__DIR__,"../src/batchUseful.jl"))

# ╔═╡ 071f028c-f123-4e3b-8c80-7495b21113af
include(joinpath(@__DIR__,"../src/OPTwrappers.jl"))

# ╔═╡ 24b257aa-bd3f-42be-8a18-169ae0cc27d9
include("../src/OPTnewEngines.jl")  # I do this to test

# ╔═╡ 6edc09ed-7b03-4429-b866-c3c8c0b70816
cd(@__DIR__)

# ╔═╡ c0ed6199-bd7c-4c2e-a5fe-545ec140a76a
dirOPTMotors=@__DIR__

# ╔═╡ 6aa8b4da-be22-4ad2-ad0c-a7cda219fb3f
include(joinpath(@__DIR__,"../src/OPTnewEngines.jl")) 
include(joinpath(@__DIR__,"../src/famousSourceFunctions.jl"))
include(joinpath(@__DIR__,"../src/famousEquations.jl"))
include(joinpath(@__DIR__,"../src/timeMarchingSchemes.jl"))

# ╔═╡ 262cf4fa-4dc1-41ba-9a89-bb7e29963240


# ╔═╡ 3d4f92bd-818a-4eed-91a5-c7c3f8fcb837


# ╔═╡ Cell order:
# ╠═6c0d05f8-b66f-11f0-84c5-f7e749a4268e
# ╠═6edc09ed-7b03-4429-b866-c3c8c0b70816
# ╠═c0ed6199-bd7c-4c2e-a5fe-545ec140a76a
# ╠═ec4d52e6-c9f1-46c4-8362-73d231fcb0e3
# ╠═cc967e2f-235a-41b0-8da9-75bbcf3e3035
# ╠═9af87ac1-2331-4b9a-a352-24381cc0a7b7
# ╠═48d5c0c6-391c-4e38-b58d-00bf4b306694
# ╠═82d0aed5-5829-448a-9f97-dd3303b1d890
# ╠═ccee7ed3-dc2f-4b97-bbab-4662120e57aa
# ╠═6aa8b4da-be22-4ad2-ad0c-a7cda219fb3f
# ╠═262cf4fa-4dc1-41ba-9a89-bb7e29963240
# ╠═071f028c-f123-4e3b-8c80-7495b21113af
# ╠═24b257aa-bd3f-42be-8a18-169ae0cc27d9
# ╠═3d4f92bd-818a-4eed-91a5-c7c3f8fcb837
