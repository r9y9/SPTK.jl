# MFCC

function mfcc(x::Vector{Cdouble}, order=20, samplerate=16000;
              α::Float64=0.97,
              eps::Float64=1.0, numfilterbunks::Int=20, cepslift::Int=22,
              usedft::Bool=false, usehamming::Bool=true,
              czero::Bool=false, power::Bool=false)

    if order > numfilterbunks - 1
        throw(ArgumentError("order must be larger than num filterbanks-1"))
    end

    # order of MFCC + 0-th + power
    cc = zeros(order+2)

    ccall((:mfcc, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble,
           Cint, Cint, Cint, Cint, Cint, Bool, Bool),
          x, cc, samplerate, α, eps,
          length(x), length(x), order+1, numfilterbunks, cepslift,
          usedft, usehamming)

    # after ccall we get
    # mfcc[0], mfcc[1], mfcc[2], ... mfcc[m-1], E(C0), Power

    if !czero && power
        cc[endof(cc)-1] = cc[endof(cc)]
    end

    if !power
        cc = cc[1:endof(cc)-1]
    end

    if !czero
        cc = cc[1:endof(cc)-1]
    end

    cc
end
