# const libSPTK = "libSPTK"

# b2mc converts MLSA filter coefficients to Mel-Cepstrum.
function b2mc(b::Vector{Float64}; alpha::Float64=0.41)
    order = length(b)-1
    mc = zeros(length(b))
    ccall((:b2mc, libSPTK), Void, (Ptr{Float64}, Ptr{Float64}, Int, Float64),
          b, mc, order, alpha)
    return mc
end

# c2ir converts cepstrum to impulse response.
function c2ir(c::Vector{Float64}, len::Int)
    h = zeros(len)
    ccall((:c2ir, libSPTK), Void, (Ptr{Float64}, Int, Ptr{Float64}, Int),
          c, length(c), h, len)
    return h
end

# gc2gc performs conversion between generalized cepstrum.
function gc2gc(c1::Vector{Float64}, gamma1::Float64, m2::Int, gamma2::Float64)
    m1 = length(c1) - 1
    c2 = zeros(m2*+1)
    ccall((:gc2gc, libSPTK), Void, (Ptr{Float64}, Int, Float64,
                                      Ptr{Float64}, Int, Float64),
          c1, m1, gamma1, c2, m2, gamma2)
    return c2
end

# gnorm performs cepstrum gain normailzation
function gnorm(c::Vector{Float64}, gamma::Float64)
    normalizedC = zeros(length(c))
    m = length(c)-1
    ccall((:gnorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                      Int, Float64),
          c, normalizedC, m, gamma)
    return normalizedC
end

# ignorm performs inverse cepstrum gain normailzation
function ignorm(normalizedC::Vector{Float64}, gamma::Float64)
    c = zeros(length(normalizedC))
    m = length(normalizedC)-1
    ccall((:gnorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                      Int, Float64),
         cnormalizedC, c, m, gamma)
    return c
end

# freqt performs frequency tranformation on cepstrum. It can be used to
# convert linear frequency cepstrum to mel frequency cepstrum.
function freqt(c::Vector{Float64}, desiredOrder::Int, alpha::Float64)
    originalOrder = length(c)-1
    transformed = zeros(desiredOrder+1)
    ccall((:freqt, libSPTK), Void,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Float64),
          c, originalOrder, transformed, desiredOrder, alpha)
    return transformed
end

# mc2b converts mel-cepstrum to MLSA filter coefficients.
function mc2b(mc::Vector{Float64}; alpha::Float64=0.41)
    order = length(mc)-1
    b = zeros(length(mc))
    ccall((:mc2b, libSPTK), Void, (Ptr{Float64}, Ptr{Float64}, Int, Float64),
          mc, b, order, alpha)
    return b
end

# mcep preforms Mel-Cepstrum analysis.
function mcep(x::Vector{Float64};
              order::Int=40, alpha::Float64=0.41, iter1::Int=2, iter2::Int=30,
              dd::Float64=0.001, etype::Int=0, e::Float64=0.0,
              f::Float64=0.0001, itype::Int=0)
    mc = zeros(order+1)
    ccall((:mcep, libSPTK), Int,
          (Ptr{Float64}, Int, Ptr{Float64}, Int,
           Float64, Int, Int, Float64, Int, Float64, Float64, Int),
          x, length(x), mc, order, alpha,
          iter1, iter2, dd, etype, e, f, itype)
    return mc
end

# mgc2mgc converts between mel log-generalized cesptrum.
function mgc2mgc(c1::Vector{Float64}, alpha1::Float64, gamma1::Float64,
                 m2::Int, alpha2::Float64, gamma2::Float64)
    c2 = zeros(m2+1)
    m1 = length(c1)-1
    ccall((:mgc2mgc, libSPTK), Void, (Ptr{Float64}, Int, Float64, Float64,
                                        Ptr{Float64}, Int, Float64, Float64),
          c1, m1, alpha1, gamma1, c2, m2, alpha2, gamma2)
    return c2
end

# Mel log-generalized cepstrum analysis
function mgcep(x::Vector{Float64};
               order::Int=40, alpha::Float64=0.41, gamma::Float64=0.0,
               n::Int=length(x)-1,
               iter1::Int=2, iter2::Int=30,
               dd::Float64=0.001, etype::Int=0, e::Float64=0.0,
               f::Float64=0.0001, itype::Int=0, otype::Int=0)
    mgc = zeros(order+1)
    ccall((:mgcep, libSPTK), Int,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Float64, Float64,
           Int, Int, Int, Float64, Int, Float64, Float64, Int),
          x, length(x), mgc, order, alpha, gamma, n,
          iter1, iter2, dd, etype, e, f, itype)

    if otype == 0 || otype == 1 || otype == 2 || otype == 4
        ccall((:ignorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                           Int, Float64),
              mgc, mgc, order, gamma)
    end

    if otype == 0 || otype == 2 || otype == 4
        ccall((:b2mc, libSPTK), Void,
              (Ptr{Float64}, Ptr{Float64}, Int, Float64),
              mgc, mgc, order, alpha)
    end

    if otype == 2 || otype == 4
        ccall((:gnorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                          Int, Float64),
              mgc, mgc, order, gamma)
    end

    if otype == 4 || otype == 5
        mgc = [mgc[1], mgc[2:end]*gamma]
    end

    return mgc
end

# swipe performs fundamental frequency (f0) estimation based on
# SWIPE - A Saw-tooth Waveform Inspired Pitch Estimation
function swipe(x::Vector{Float64}, samplerate::Int;
               hopsize::Int=80,
               min::Float64=50.0, max::Float64=800.0,
               st::Float64=0.3, otype::Int=1)
    expected_len = length(x)/hopsize + 1
    f0 = zeros(int(expected_len))
    ccall((:swipe, libSPTK), Void,
          (Ptr{Float64}, Ptr{Float64}, Int, Int, Int, Float64, Float64,
           Float64, Int),
          x, f0, length(x), samplerate, hopsize,
          min, max, st, otype)
    return f0
end

# MFCC computes Mel-Frequency Cepstrum Coefficients.
function mfcc(x::Vector{Float64};
              order::Int=20,
              samplerate::Float64=16000.0, alpha::Float64=0.97,
              eps::Float64=1.0, numfilterbunks::Int=20, cepslift::Int=22,
              usedft::Bool=false, usehamming::Bool=true,
              czero::Bool=false, power::Bool=false)

    @assert order+1 <= numfilterbunks

    # order of MFCC + 0-th + power
    cc = zeros(order+2)

    ccall((:mfcc, libSPTK), Void,
          (Ptr{Float64}, Ptr{Float64}, Float64, Float64, Float64,
           Int, Int, Int, Int, Int, Bool, Bool),
          x, cc, samplerate, alpha, eps,
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

    return cc
end
