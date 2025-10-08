"""
    myInclude(file::AbstractString)

Include or track a file depending on whether `Revise` is active.

- If Revise is loaded, calls `Revise.includet(file)`
- Otherwise, calls plain `include(file)`

This lets you develop interactively with Revise and
run in production mode without modification.
"""

module BatchRevise

    export myInclude, using_revise

    const using_revise = isdefined(Main, :Revise)

    function myInclude(file::AbstractString)
        if using_revise
            # Ensure Revise is available and callable
            if hasproperty(Main, :Revise)
                @info "Including with Revise: $file"
                Main.Revise.includet(file)
            else
                @warn "Revise not found; falling back to include for $file"
                include(file)
            end
        else
            include(file)
        end
    end

end 