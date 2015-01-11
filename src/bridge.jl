# Defined in deps.jl
@assert isdefined(:libSPTK)

# mcep preforms Mel-Cepstrum analysis.
function mcep(x::Vector{Float64}, order::Int=40, α::Float64=0.41;
              miniter::Int=2, maxiter::Int=30,
              dd::Float64=0.001, etype::Int=0, e::Float64=0.0,
              f::Float64=0.0001, itype::Int=0)
    mc = zeros(order+1)
    ccall((:mcep, libSPTK), Int,
          (Ptr{Float64}, Int, Ptr{Float64}, Int,
           Float64, Int, Int, Float64, Int, Float64, Float64, Int),
          x, length(x), mc, order, α,
          miniter, maxiter, dd, etype, e, f, itype)
    mc
end

# gcep performs generalized cesptrum analysis.
function gcep(x::Vector{Float64}, order::Int=40, γ::Float64=0.0;
              miniter::Int=2, maxiter::Int=30, d::Float64=0.001,
              etype::Int=0, e::Float64=0.0, f::Float64=0.000001,
              itype::Int=0)
    gc = zeros(order+1)
    ccall((:gcep, libSPTK), Int,
          (Ptr{Float64}, Int, Ptr{Float64}, Int,
           Float64, Int, Int, Float64, Int, Float64, Float64, Int),
          x, length(x), gc, order, γ,
          miniter, maxiter, d, etype, e, f, itype)
    gc
end

# mgcep performs Mel log-generalized cepstrum analysis.
function mgcep(x::Vector{Float64}, order::Int=40,
               α::Float64=0.41, γ::Float64=0.0;
               n::Int=length(x)-1,
               miniter::Int=2, maxiter::Int=30,
               dd::Float64=0.001, etype::Int=0, e::Float64=0.0,
               f::Float64=0.0001, itype::Int=0, otype::Int=0)
    mgc = zeros(order+1)
    ccall((:mgcep, libSPTK), Int,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Float64, Float64,
           Int, Int, Int, Float64, Int, Float64, Float64, Int),
          x, length(x), mgc, order, α, γ, n,
          miniter, maxiter, dd, etype, e, f, itype)

    if otype == 0 || otype == 1 || otype == 2 || otype == 4
        ccall((:ignorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                           Int, Float64),
              mgc, mgc, order, γ)
    end

    if otype == 0 || otype == 2 || otype == 4
        ccall((:b2mc, libSPTK), Void,
              (Ptr{Float64}, Ptr{Float64}, Int, Float64),
              mgc, mgc, order, α)
    end

    if otype == 2 || otype == 4
        ccall((:gnorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                          Int, Float64),
              mgc, mgc, order, γ)
    end

    if otype == 4 || otype == 5
        mgc = [mgc[1], mgc[2:end]*γ]
    end

    mgc
end

# uels performs unbiased estimation of target log spectrum.
function uels(x::Vector{Float64}, order::Int;
              miniter::Int=2, maxiter::Int=30, dd::Float64=0.001,
              etype::Int=0, e::Float64=0.0, f::Float64=0.0001, itype::Int=0)
    c = zeros(order+1)
    ccall((:uels, libSPTK), Int,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Int, Int, Float64, Int,
           Float64, Float64, Int),
          x, length(x), c, order,
          miniter, maxiter, dd, etype, e, f, itype)
    c
end

# fftcep computes cepstrum from log spectrum (that can be computed using FFT).
function fftcep(logsp::Vector{Float64}, order::Int;
                itr::Int=0, accelerationf::Float64=0.0)
    c = zeros(order+1)

    ccall((:fftcep, libSPTK), Void,
          (Ptr{Float64}, Int64, Ptr{Float64}, Int64, Int64, Float64),
          logsp, length(logsp), c, length(c), itr, accelerationf)
    c
end

# Matrix input version of mfcc. For given audio frames, this function returns
# mfcc spectrogram.
function mfcc(x::Matrix{Float64}, order::Int=20, samplerate::Int=16000;
              α::Float64=0.97,
              eps::Float64=1.0, numfilterbunks::Int=20, cepslift::Int=22,
              usedft::Bool=false, usehamming::Bool=true,
              czero::Bool=false, power::Bool=false)
    rows = order
    if czero
        rows += 1
    end
    if power
        rows += 1
    end

    const T = size(x, 2)
    ccs = Array(Float64, rows, T)
    for t=1:T
        ccs[:,t] = mfcc(x[:,t], order, samplerate,
                        α=α, eps=eps, numfilterbunks=numfilterbunks,
                        cepslift=cepslift, usedft=usedft,
                        usehamming=usehamming, czero=czero, power=power)
    end
    ccs
end

# mfcc computes Mel-Frequency Cepstrum Coefficients using DCT.
function mfcc(x::Vector{Float64}, order::Int=20, samplerate::Int=16000;
              α::Float64=0.97,
              eps::Float64=1.0, numfilterbunks::Int=20, cepslift::Int=22,
              usedft::Bool=false, usehamming::Bool=true,
              czero::Bool=false, power::Bool=false)

    @assert order+1 <= numfilterbunks

    # order of MFCC + 0-th + power
    cc = zeros(order+2)

    ccall((:mfcc, libSPTK), Void,
          (Ptr{Float64}, Ptr{Float64}, Float64, Float64, Float64,
           Int, Int, Int, Int, Int, Bool, Bool),
          x, cc, float64(samplerate), α, eps,
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

# mc2b converts mel-cepstrum to MLSA filter coefficients.
function mc2b(mc::Vector{Float64}, α::Float64=0.41)
    order = length(mc)-1
    b = zeros(length(mc))
    ccall((:mc2b, libSPTK), Void, (Ptr{Float64}, Ptr{Float64}, Int, Float64),
          mc, b, order, α)
    b
end

# b2mc converts MLSA filter coefficients to Mel-Cepstrum.
function b2mc(b::Vector{Float64}, α::Float64=0.41)
    order = length(b)-1
    mc = zeros(length(b))
    ccall((:b2mc, libSPTK), Void, (Ptr{Float64}, Ptr{Float64}, Int, Float64),
          b, mc, order, α)
    mc
end

# c2ir converts cepstrum to impulse response.
function c2ir(c::Vector{Float64}, len::Int)
    h = zeros(len)
    ccall((:c2ir, libSPTK), Void, (Ptr{Float64}, Int, Ptr{Float64}, Int),
          c, length(c), h, len)
    h
end


# c2ndps converts cepstrum to negative derivative of phase spectrum.
function c2ndps(c::Vector{Float64}, fftlen::Int)
    ndps = zeros(fftlen)
    m = length(c)-1
    ccall((:c2ndps, libSPTK), Void, (Ptr{Float64}, Int, Ptr{Float64}, Int),
          c, m, ndps, fftlen)
    ndps[1:fftlen>>1+1]
end

# c2ndps converts negative derivative of phase spectrum to cepstrum.
function ndps2c(ndps::Vector{Float64}, order::Int)
    fftlen = (length(ndps)-1)*2 # assuming the length of npds is fftsize/2+1
    c = zeros(order+1)
    ccall((:ndps2c, libSPTK), Void, (Ptr{Float64}, Int, Ptr{Float64}, Int),
          ndps, fftlen, c, order)
    c
end

# gc2gc performs conversion between generalized cepstrum.
function gc2gc(c1::Vector{Float64}, γ₁::Float64, m2::Int, γ₂::Float64)
    m1 = length(c1) - 1
    c2 = zeros(m2+1)
    ccall((:gc2gc, libSPTK), Void, (Ptr{Float64}, Int, Float64,
                                      Ptr{Float64}, Int, Float64),
          c1, m1, γ₁, c2, m2, γ₂)
    c2
end

# gnorm performs cepstrum gain normailzation
function gnorm(c::Vector{Float64}, γ::Float64)
    normalizedC = zeros(length(c))
    m = length(c)-1
    ccall((:gnorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                      Int, Float64),
          c, normalizedC, m, γ)
    normalizedC
end

# ignorm performs inverse cepstrum gain normailzation
function ignorm(normalizedC::Vector{Float64}, γ::Float64)
    c = zeros(length(normalizedC))
    m = length(normalizedC)-1
    ccall((:ignorm, libSPTK), Void, (Ptr{Float64}, Ptr{Float64},
                                      Int, Float64),
         normalizedC, c, m, γ)
    c
end

# freqt performs frequency tranformation on cepstrum. It can be used to
# convert linear frequency cepstrum to mel frequency cepstrum.
function freqt(c::Vector{Float64}, order::Int, α::Float64)
    org_order = length(c)-1
    transformed = zeros(order+1)
    ccall((:freqt, libSPTK), Void,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Float64),
          c, org_order, transformed, order, α)
    transformed
end

function frqtr(c::Vector{Float64}, order::Int, α::Float64)
    org_order = length(c)-1
    transformed = zeros(order+1)
    ccall((:frqtr, libSPTK), Void,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Float64),
          c, org_order, transformed, order, α)
    transformed
end

# mgc2mgc converts between mel log-generalized cesptrum.
function mgc2mgc(c1::Vector{Float64}, α₁::Float64, γ₁::Float64,
                 m2::Int, α₂::Float64, γ₂::Float64)
    c2 = zeros(m2+1)
    m1 = length(c1)-1
    ccall((:mgc2mgc, libSPTK), Void, (Ptr{Float64}, Int, Float64, Float64,
                                      Ptr{Float64}, Int, Float64, Float64),
          c1, m1, α₁, γ₁, c2, m2, α₂, γ₂)
    c2
end

# mgc2sp converts mel generalized cepstrum to log spectrum
function mgc2sp(mgc::Vector{Float64}, α::Float64, γ::Float64, fftlen::Int)
    order = length(mgc)-1
    sp = Array(Complex{Float64}, fftlen>>1+1)
    sp_r = zeros(Float64, fftlen)
    sp_i = zeros(Float64, fftlen)
    ccall((:mgc2sp, libSPTK), Void,
          (Ptr{Float64}, Int, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Int),
          mgc, order, α, γ, sp_r, sp_i, fftlen)
    for i=1:length(sp)
        @inbounds sp[i] = Complex(sp_r[i], sp_i[i])
    end
    sp
end

# mgclsp2sp converts mgc-lsp to spectrum.
function mgclsp2sp(lsp::Vector{Float64}, α::Float64, γ::Float64,
                   fftlen::Int;
                   gain::Bool=true)
    sp = zeros(fftlen>>1 + 1)
    m = gain ? length(lsp)-1 : length(lsp)
    ccall((:mgclsp2sp, libSPTK), Void, (Float64, Float64, Ptr{Float64}, Int,
                                        Ptr{Float64}, Int, Int),
          α, γ, lsp, m, sp, length(sp), int(gain))
    sp
end

# swipe performs fundamental frequency (f0) estimation based on
# SWIPE - A Saw-tooth Waveform Inspired Pitch Estimation
function swipe(x::Vector{Float64}, samplerate::Int, hopsize::Int=80;
               min::Float64=50.0, max::Float64=800.0,
               st::Float64=0.3, otype::Int=1)
    expectedlen = div(length(x), hopsize) + 1
    f0 = zeros(expectedlen)
    ccall((:swipe, libSPTK), Void,
          (Ptr{Float64}, Ptr{Float64}, Int, Int, Int, Float64, Float64,
           Float64, Int),
          x, f0, length(x), samplerate, hopsize,
          min, max, st, otype)
    f0
end

# mlsadf performs Mel Log Spectrum Approximation (MLSA) digital filtering.
function mlsadf(x::Float64, b::Vector{Float64}, α::Float64, pd::Int,
                delay::Vector{Float64})
    ccall((:mlsadf, libSPTK), Float64,
          (Float64, Ptr{Float64}, Int, Float64, Int, Ptr{Float64}),
          x, b, length(b)-1, α, pd, delay)
end

# see mlsadf.c in original SPTK for this magic allocation
mlsadf_delay(order::Int, pd::Int) = zeros(3*(pd+1) + pd*(order+2))

# mglsadf performs Mel Generalized Log Spectrum Approximation (MGLSA) digital
# filtering.
function mglsadf(x::Float64, b::Vector{Float64}, α::Float64, stage::Int,
                 delay::Vector{Float64})
    ccall((:mglsadf, libSPTK), Float64,
          (Float64, Ptr{Float64}, Int, Float64, Int, Ptr{Float64}),
          x, b, length(b)-1, α, stage, delay)
end

# see mglsadf.c in original SPTK for this magic allocation
mglsadf_delay(order::Int, stage::Int) = zeros((order+1)*stage)

immutable WindowType
    w::Int
end

function _window!(wtype::WindowType, x::Vector{Float64}; normalize::Int=0)
    # normalize: Int
    #    0 : don't normalize
    #    1 : normalize by power
    #    2 : normalize by magnitude
    @assert normalize == 0 || normalize == 1 || normalize == 2
    g = ccall((:window, libSPTK), Float64,
              (WindowType, Ptr{Float64}, Int, Int), wtype, x, length(x),
              normalize)
    x
end

for (f, wtype) in [(:blackman, 0),
                   (:hamming, 1),
                   (:hanning, 2),
                   (:bartlett, 3),
                   (:trapezoid, 4),
                   (:rectangular, 5)]
    @eval begin
        function $f(n::Integer; normalize::Int=0)
            y = ones(Float64, n)
            _window!(WindowType($wtype), y; normalize=normalize)
        end
    end
end

function theq(t::Vector{Float64}, h::Vector{Float64}, a::Vector{Float64},
              b::Vector{Float64}, n::Int, e::Float64)
    ccall((:theq, libSPTK), Int,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int,
           Float64), t, h, a, b, n, e)
end

function b2c(c::Vector{Float64}, order::Int, α::Float64)
    org_order = length(c)-1
    transformed = zeros(order+1)
    ccall((:b2c, libSPTK), Void,
          (Ptr{Float64}, Int, Ptr{Float64}, Int, Float64),
          c, org_order, transformed, order, α)
    transformed
end

# extend real-vector to real-vector transformation for matrix input
for f in [:mcep,
          :gcep,
          :mgcep,
          :uels,
          :fftcep,
          :mfcc,
          :mc2b,
          :b2mc,
          :c2ir,
          :c2ndps,
          :ndps2c,
          :gc2gc,
          :gnorm,
          :ignorm,
          :freqt,
          :mgc2mgc,
          :mgclsp2sp,
          ]
    @eval begin
        function $f(x::Matrix{Float64}, args...; kargs...)
            r = $f(x[:, 1], args...; kargs...)
            ret = Array(eltype(x), size(r, 1), size(x, 2))
            for i = 1:length(r)
                @inbounds ret[i, 1] = r[i]
            end
            for i = 2:size(x, 2)
                @inbounds ret[:, i] = $f(x[:, i], args...; kargs...)
            end
            ret
        end
    end
end

# extend vector to complex vector transformation for matrix input
for f in [:mgc2sp]
    @eval begin
        function $f(x::Matrix{Float64}, args...; kargs...)
            r = $f(x[:, 1], args...; kargs...)
            ret = Array(Complex{eltype(x)}, size(r, 1), size(x, 2))
            for i = 1:length(r)
                @inbounds ret[i, 1] = r[i]
            end
            for i = 2:size(x, 2)
                @inbounds ret[:, i] = $f(x[:, i], args...; kargs...)
            end
            ret
        end
    end
end
