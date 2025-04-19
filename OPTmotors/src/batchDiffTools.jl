using SparseDiffTools

function myForwarddiff_color_jacobian!(J, Res_closed!, U, colorvec = colors)
    return forwarddiff_color_jacobian!(J, Res_closed!, U, colorvec = colors)
end