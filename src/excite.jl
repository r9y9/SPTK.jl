### Source excitation generation ###

function excite(pitch::StridedVector{Cdouble}, hopsize=80;
                interp_period::Int=1,
                gaussian::Bool=false,
                seed::Int=1)
    expectedlen = trunc(hopsize * (length(pitch) - 1))
    excitation = Array{Cdouble}(expectedlen)

    ccall((:excite, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Bool, Cint),
          pitch, length(pitch), excitation, hopsize, interp_period,
          gaussian, seed)

    excitation
end
