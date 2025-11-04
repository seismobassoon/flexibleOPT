"""
    myInclude(file::AbstractString)

Include or track a file depending on whether `Revise` is active.

Revise will be automatically loaded (if the package exists in the environment)

This lets you develop interactively with Revise and
run in production mode without modification.
"""




# now I don't use Revise.jl (incompatible with Pluto.jl)
function myInclude(file::AbstractString)
    absfile = joinpath(@__DIR__, file)
    Base.include_dependency(absfile)  # so Pluto tracks it
    return Base.include(Main, absfile) # load into Main (so visible across cells)
end
