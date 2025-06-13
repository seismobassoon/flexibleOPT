using DrWatson

using SHA

function hash_parameters(params)

    # Serialize as string and hash
    str = repr(filtered)
    return bytes2hex(sha1(str))[1:8]  # Short hash
end
