using DrWatson

using SHA

function hash_parameters(params)
    # Serialize as string and hash
    filtered = Dict(k => v for (k,v) in params)
    str = repr(filtered)
    return bytes2hex(sha1(str))[1:8]  # Short hash
end

function myProduceOrLoad(functionName,paramDict,directoryName::String)
    output=myProduceOrLoad(functionName,paramDict,directoryName,directoryName)
    return output
end

function myProduceOrLoad(functionName,paramDict,directoryName::String,prefixName::String)

    # this is myFunction strategy for DrWatson

    hash_id = hash_parameters(paramDict)
    

    newDict = Dict{String,Any}(paramDict)
    newDict["hash_id"] = hash_id
    
    output, _ = produce_or_load(functionName,newDict,datadir(directoryName);filename = config -> savename(prefixName,"",hash_id))
    
    return output

end