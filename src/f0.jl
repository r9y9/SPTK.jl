# F0 estimation

function swipe(x::Vector{Cdouble}, samplerate, hopsize=80;
               min::Float64=50.0, max::Float64=800.0,
               st::Float64=0.3, otype::Int=1)
    expectedlen = div(length(x), hopsize) + 1
    f0 = zeros(expectedlen)
    ccall((:swipe, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble, Cdouble,
           Cdouble, Cint),
          x, f0, length(x), samplerate, hopsize,
          min, max, st, otype)
    f0
end
