function agexp(r, x, y)
    ccall((:agexp, libSPTK), Cdouble, (Cdouble, Cdouble, Cdouble), r, x, y)
end

gexp(r, x) = ccall((:gexp, libSPTK), Cdouble, (Cdouble, Cdouble), r, x)
glog(r, x) = ccall((:glog, libSPTK), Cdouble, (Cdouble, Cdouble), r, x)

function cholesky(c::Vector{Cdouble}, a::Vector{Cdouble}, b::Vector{Cdouble};
                  eps::Float64=1.0e-6)
    if length(c) != length(a) || length(c) != length(b)
        error("input vectors should have same length")
    end
    n = length(c)
    ret = ccall((:cholesky, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
                c, a, b, n, eps)
    if ret != 0
        error("failed to compute choleskey decomposition")
    end
end

# mcep preforms Mel-Cepstrum analysis.
function mcep(x::Vector{Cdouble}, order=40, α=0.41;
              miniter::Int=2, maxiter::Int=30,
              dd::Float64=0.001, etype::Int=0, e::Float64=0.0,
              f::Float64=0.0001, itype::Int=0)
    mc = zeros(order+1)
    ccall((:mcep, libSPTK), Cint,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint,
           Cdouble, Cint, Cint, Cdouble, Cint, Cdouble, Cdouble, Cint),
          x, length(x), mc, order, α,
          miniter, maxiter, dd, etype, e, f, itype)
    mc
end

# gcep performs generalized cesptrum analysis.
function gcep(x::Vector{Cdouble}, order=40, γ=0.0;
              miniter::Int=2, maxiter::Int=30, d::Float64=0.001,
              etype::Int=0, e::Float64=0.0, f::Float64=0.000001,
              itype::Int=0)
    gc = zeros(order+1)
    ccall((:gcep, libSPTK), Cint,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint,
           Cdouble, Cint, Cint, Cdouble, Cint, Cdouble, Cdouble, Cint),
          x, length(x), gc, order, γ,
          miniter, maxiter, d, etype, e, f, itype)
    gc
end

# mgcep performs Mel log-generalized cepstrum analysis.
function mgcep(x::Vector{Cdouble}, order=40, α=0.41, γ=0.0;
               n::Int=length(x)-1,
               miniter::Int=2, maxiter::Int=30,
               dd::Float64=0.001, etype::Int=0, e::Float64=0.0,
               f::Float64=0.0001, itype::Int=0, otype::Int=0)
    mgc = zeros(order+1)
    ccall((:mgcep, libSPTK), Cint,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble, Cdouble,
           Cint, Cint, Cint, Cdouble, Cint, Cdouble, Cdouble, Cint),
          x, length(x), mgc, order, α, γ, n,
          miniter, maxiter, dd, etype, e, f, itype)

    if otype == 0 || otype == 1 || otype == 2 || otype == 4
        ccall((:ignorm, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble},
                                           Cint, Cdouble),
              mgc, mgc, order, γ)
    end

    if otype == 0 || otype == 2 || otype == 4
        ccall((:b2mc, libSPTK), Void,
              (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
              mgc, mgc, order, α)
    end

    if otype == 2 || otype == 4
        ccall((:gnorm, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble},
                                          Cint, Cdouble),
              mgc, mgc, order, γ)
    end

    if otype == 4 || otype == 5
        mgc = [mgc[1], mgc[2:end]*γ]
    end

    mgc
end

# uels performs unbiased estimation of target log spectrum.
function uels(x::Vector{Cdouble}, order;
              miniter::Int=2, maxiter::Int=30, dd::Float64=0.001,
              etype::Int=0, e::Float64=0.0, f::Float64=0.0001, itype::Int=0)
    c = zeros(order+1)
    ccall((:uels, libSPTK), Cint,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble, Cint,
           Cdouble, Cdouble, Cint),
          x, length(x), c, order,
          miniter, maxiter, dd, etype, e, f, itype)
    c
end

# fftcep computes cepstrum from log spectrum (that can be computed using FFT).
function fftcep(logsp::Vector{Cdouble}, order;
                itr::Int=0, accelerationf::Float64=0.0)
    c = zeros(order+1)

    ccall((:fftcep, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Cdouble),
          logsp, length(logsp), c, length(c), itr, accelerationf)
    c
end

# mfcc computes Mel-Frequency Cepstrum Coefficients using DCT.
function mfcc(x::Vector{Cdouble}, order=20, samplerate=16000;
              α::Float64=0.97,
              eps::Float64=1.0, numfilterbunks::Int=20, cepslift::Int=22,
              usedft::Bool=false, usehamming::Bool=true,
              czero::Bool=false, power::Bool=false)

    @assert order+1 <= numfilterbunks

    # order of MFCC + 0-th + power
    cc = zeros(order+2)

    ccall((:mfcc, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble,
           Cint, Cint, Cint, Cint, Cint, Bool, Bool),
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
function mc2b(mc::Vector{Cdouble}, α=0.41)
    order = length(mc)-1
    b = zeros(length(mc))
    ccall((:mc2b, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          mc, b, order, α)
    b
end

# b2mc converts MLSA filter coefficients to Mel-Cepstrum.
function b2mc(b::Vector{Cdouble}, α=0.41)
    order = length(b)-1
    mc = zeros(length(b))
    ccall((:b2mc, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          b, mc, order, α)
    mc
end

# c2ir converts cepstrum to impulse response.
function c2ir(c::Vector{Cdouble}, len)
    h = zeros(len)
    ccall((:c2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, length(c), h, len)
    h
end


# c2ndps converts cepstrum to negative derivative of phase spectrum.
function c2ndps(c::Vector{Cdouble}, fftlen)
    ndps = zeros(fftlen)
    m = length(c)-1
    ccall((:c2ndps, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, m, ndps, fftlen)
    ndps[1:fftlen>>1+1]
end

# c2ndps converts negative derivative of phase spectrum to cepstrum.
function ndps2c(ndps::Vector{Cdouble}, order)
    fftlen = (length(ndps)-1)*2 # assuming the length of npds is fftsize/2+1
    c = zeros(order+1)
    ccall((:ndps2c, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          ndps, fftlen, c, order)
    c
end

# gc2gc performs conversion between generalized cepstrum.
function gc2gc(c1::Vector{Cdouble}, γ₁, m2, γ₂)
    m1 = length(c1) - 1
    c2 = zeros(m2+1)
    ccall((:gc2gc, libSPTK), Void, (Ptr{Cdouble}, Cint, Cdouble,
                                    Ptr{Cdouble}, Cint, Cdouble),
          c1, m1, γ₁, c2, m2, γ₂)
    c2
end

function lpc(x::Vector{Cdouble}, order;
             f::Float64=1e-6)
    a = Array(Cdouble, order+1)
    ccall((:lpc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          x, length(x), a, order, f)
    a
end

function lpc2c(a::Vector{Cdouble})
    order = length(a) - 1
    c = Array(Cdouble, order+1)
    ccall((:lpc2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), a, order, c, order)
    c
end

function lpc2lsp(lpc::Vector{Cdouble}, order;
                 numsp::Int=128,
                 maxiter::Int=4,
                 e::Float64=1e-6,
                 loggain::Bool=true,
                 otype::Int=0,
                 fs=nothing
    )
    lsp = zeros(Cdouble, order+1)
    ccall((:lpc2lsp, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble),
          lpc, sub(lsp, 2:length(lsp)), order, numsp, maxiter, e)

    if otype == 0
        for i=2:length(lsp)
            @inbounds lsp[i] *= 2π
        end
    elseif otype == 2 || otype == 3
        fs == nothing && error("fs must be specified when otype == 2 or 3")
        for i=2:length(lsp)
            @inbounds lsp[i] *= fs
        end
    end

    # this is really bad ...
    if otype == 3
        for i=2:length(lsp)
            @inbounds lsp[i] *= 1000.0
        end
    end

    if loggain
        lsp[1] = log(lpc[1])
    else
        lsp[1] = lpc[1]
    end

    lsp
end

function lpc2par(lpc::Vector{Cdouble})
    par = similar(lpc)
    ccall((:lpc2par, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
          lpc, par, length(lpc)-1)
    par
end

function lsp2sp(lsp::Vector{Cdouble}, fftlen)
    # assume lsp has loggain at lsp[1]
    sp = Array(Cdouble, fftlen>>1+1)
    ccall((:lsp2sp, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint),
          lsp, length(lsp)-1, sp, length(sp), 1)
    sp
end

function lspcheck(lpc::Vector{Cdouble})
    r = ccall((:lspcheck, libSPTK),
              Cint, (Ptr{Cdouble}, Cint), lpc, length(lpc)-1)::Cint
    ifelse(r == 0, true, false)
end

# gnorm performs cepstrum gain normailzation
function gnorm(c::Vector{Cdouble}, γ)
    normalizedC = zeros(length(c))
    m = length(c)-1
    ccall((:gnorm, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble},
                                      Cint, Cdouble),
          c, normalizedC, m, γ)
    normalizedC
end

# ignorm performs inverse cepstrum gain normailzation
function ignorm(normalizedC::Vector{Cdouble}, γ)
    c = zeros(length(normalizedC))
    m = length(normalizedC)-1
    ccall((:ignorm, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble},
                                      Cint, Cdouble),
         normalizedC, c, m, γ)
    c
end

# freqt performs frequency tranformation on cepstrum. It can be used to
# convert linear frequency cepstrum to mel frequency cepstrum.
function freqt(c::Vector{Cdouble}, order, α)
    org_order = length(c)-1
    transformed = zeros(order+1)
    ccall((:freqt, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          c, org_order, transformed, order, α)
    transformed
end

function frqtr(c::Vector{Cdouble}, order, α)
    org_order = length(c)-1
    transformed = zeros(order+1)
    ccall((:frqtr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          c, org_order, transformed, order, α)
    transformed
end

# mgc2mgc converts between mel log-generalized cesptrum.
function mgc2mgc(c1::Vector{Cdouble}, α₁, γ₁, m2, α₂, γ₂)
    c2 = zeros(m2+1)
    m1 = length(c1)-1
    ccall((:mgc2mgc, libSPTK), Void, (Ptr{Cdouble}, Cint, Cdouble, Cdouble,
                                      Ptr{Cdouble}, Cint, Cdouble, Cdouble),
          c1, m1, α₁, γ₁, c2, m2, α₂, γ₂)
    c2
end

# mgc2sp converts mel generalized cepstrum to log spectrum
function mgc2sp(mgc::Vector{Cdouble}, α, γ, fftlen)
    order = length(mgc)-1
    sp = Array(Complex{Cdouble}, fftlen>>1+1)
    sp_r = zeros(Cdouble, fftlen)
    sp_i = zeros(Cdouble, fftlen)
    ccall((:mgc2sp, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Cint),
          mgc, order, α, γ, sp_r, sp_i, fftlen)
    for i=1:length(sp)
        @inbounds sp[i] = Complex(sp_r[i], sp_i[i])
    end
    sp
end

# mgclsp2sp converts mgc-lsp to spectrum.
function mgclsp2sp(lsp::Vector{Cdouble}, α, γ, fftlen;
                   gain::Bool=true)
    sp = zeros(fftlen>>1 + 1)
    m = gain ? length(lsp)-1 : length(lsp)
    ccall((:mgclsp2sp, libSPTK), Void, (Cdouble, Cdouble, Ptr{Cdouble}, Cint,
                                        Ptr{Cdouble}, Cint, Cint),
          α, γ, lsp, m, sp, length(sp), int(gain))
    sp
end

# swipe performs fundamental frequency (f0) estimation based on
# SWIPE - A Saw-tooth Waveform Inspired Pitch Estimation
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

function poledf(x::Cdouble, a::Vector{Cdouble}, delay::Vector{Cdouble})
    ccall((:poledf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, a, length(a)-1, delay)
end

poledf_delay(order) = zeros(order)

function lmadf(x::Cdouble, b::Vector{Cdouble}, pd, delay::Vector{Cdouble})
    ccall((:lmadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, pd, delay)
end

lmadf_delay(order, pd) = zeros(2pd*(order+1))

# mlsadf performs Mel Log Spectrum Approximation (MLSA) digital filtering.
function mlsadf(x::Cdouble, b::Vector{Cdouble}, α, pd, delay::Vector{Cdouble})
    ccall((:mlsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, α, pd, delay)
end

# see mlsadf.c in original SPTK for this magic allocation
mlsadf_delay(order, pd) = zeros(3*(pd+1) + pd*(order+2))

# mglsadf performs Mel Generalized Log Spectrum Approximation (MGLSA) digital
# filtering.
function mglsadf(x::Cdouble, b::Vector{Cdouble}, α, stage,
                 delay::Vector{Cdouble})
    ccall((:mglsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, α, stage, delay)
end

# see mglsadf.c in original SPTK for this magic allocation
mglsadf_delay(order, stage) = zeros((order+1)*stage)

immutable WindowType
    w::Int
end

function _window!(wtype::WindowType, x::Vector{Cdouble}; normalize::Int=0)
    # normalize: Int
    #    0 : don't normalize
    #    1 : normalize by power
    #    2 : normalize by magnitude
    @assert normalize == 0 || normalize == 1 || normalize == 2
    g = ccall((:window, libSPTK), Cdouble,
              (WindowType, Ptr{Cdouble}, Cint, Cint), wtype, x, length(x),
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

function theq(t::Vector{Cdouble}, h::Vector{Cdouble}, a::Vector{Cdouble},
              b::Vector{Cdouble}, n, e)
    ccall((:theq, libSPTK), Cint,
          (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint,
           Cdouble), t, h, a, b, n, e)
end

function b2c(c::Vector{Cdouble}, order, α)
    org_order = length(c)-1
    transformed = zeros(order+1)
    ccall((:b2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          c, org_order, transformed, order, α)
    transformed
end
