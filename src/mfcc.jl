# MFCC

function mfcc(x::StridedVector{Cdouble}, order=13, samplerate=16000;
              α::Float64=0.97,
              eps::Float64=1.0,
              windowlen::Int=length(x),
              framelen::Int=length(x),
              numfilterbunks::Int=20,
              cepslift::Int=22,
              usedft::Bool=true,
              usehamming::Bool=true,
              czero::Bool=false,
              power::Bool=false)
    if order+1 > numfilterbunks
        throw(ArgumentError("order+1 must be less than or equal to the number of filterbanks"))
    end

    # order of MFCC + 0-th + power
    cc = zeros(order+2)

    ccall((:mfcc, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble,
           Cint, Cint, Cint, Cint, Cint, Bool, Bool),
          x, cc, samplerate, α, eps,
          windowlen, framelen, order+1, numfilterbunks, cepslift,
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
