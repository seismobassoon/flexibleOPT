"""
    myInclude(file::AbstractString)

Include or track a file depending on whether `Revise` is active.

- If Revise is loaded, calls `Revise.includet(file)`
- Otherwise, calls plain `myInclude(file)`

This lets you develop interactively with Revise and
run in production mode without modification.
"""
#module BatchRevise

    export myInclude, using_revise

    # Try to activate Revise automatically
    const using_revise = try
        @eval using Revise
        true
    catch
        @warn "Revise not found; falling back to myInclude."
        false
    end

    "Smart myInclude: uses Revise.includet if Revise is active, otherwise myInclude()."
    function myInclude(file::AbstractString)
        if using_revise && hasproperty(Main, :Revise)
            @info "Including with Revise: $file"
            @eval Main Revise.includet($file)
        else
            include(file)
        end
    end

#end # module
