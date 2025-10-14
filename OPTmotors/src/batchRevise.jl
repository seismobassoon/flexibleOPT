"""
    myInclude(file::AbstractString)

Include or track a file depending on whether `Revise` is active.

Revise will be automatically loaded (if the package exists in the environment)

This lets you develop interactively with Revise and
run in production mode without modification.
"""


# Try to activate Revise automatically
using_revise = try
    @eval using Revise
    true
catch
    @warn "Revise not found; falling back to myInclude."
    false
end

"Smart myInclude: uses Revise.includet if Revise is active, otherwise include()."
function myInclude(file::AbstractString)
    if using_revise #&& hasproperty(Main, :Revise)
        @info "Including with Revise: $file"
        @eval Main Revise.includet($file)
    else
        include(file)
    end
end

