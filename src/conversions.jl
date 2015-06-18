# Conversions

function mc2b(mc::Vector{Cdouble}, α=0.41)
    order = length(mc)-1
    b = zeros(length(mc))
    ccall((:mc2b, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          mc, b, order, α)
    b
end

function b2mc(b::Vector{Cdouble}, α=0.41)
    order = length(b)-1
    mc = zeros(length(b))
    ccall((:b2mc, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          b, mc, order, α)
    mc
end

function c2acr(c::Vector{Cdouble}, m2, fftlen)
    r = zeros(m2 + 1)
    ccall((:c2acr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint),
          c, length(c) - 1, r, m2, fftlen)
    r
end

function c2ir(c::Vector{Cdouble}, len)
    h = zeros(len)
    ccall((:c2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, length(c), h, len)
    h
end

function ic2ir(h::Vector{Cdouble}, order)
    c = zeros(order+1)
    ccall((:ic2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          h, length(h), c, length(c))
    c
end

function c2ndps(c::Vector{Cdouble}, fftlen)
    ndps = zeros(fftlen)
    m = length(c)-1
    ccall((:c2ndps, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, m, ndps, fftlen)
    ndps[1:fftlen>>1+1]
end

function ndps2c(ndps::Vector{Cdouble}, order)
    fftlen = (length(ndps)-1)*2 # assuming the length of npds is fftsize/2+1
    c = zeros(order+1)
    ccall((:ndps2c, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          ndps, fftlen, c, order)
    c
end

function gc2gc(c1::Vector{Cdouble}, γ₁, m2, γ₂)
    m1 = length(c1) - 1
    c2 = zeros(m2+1)
    ccall((:gc2gc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Ptr{Cdouble}, Cint, Cdouble),
          c1, m1, γ₁, c2, m2, γ₂)
    c2
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
                 fs=nothing)
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

    # this is really ugly...
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

function gnorm(c::Vector{Cdouble}, γ)
    normalizedC = zeros(length(c))
    m = length(c)-1
    ccall((:gnorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), c, normalizedC, m, γ)
    normalizedC
end

function ignorm(normalizedC::Vector{Cdouble}, γ)
    c = zeros(length(normalizedC))
    m = length(normalizedC)-1
    ccall((:ignorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), normalizedC, c, m, γ)
    c
end

function freqt(c::Vector{Cdouble}, order, α)
    org_order = length(c)-1
    converted = zeros(order+1)
    ccall((:freqt, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          c, org_order, converted, order, α)
    converted
end

function frqtr(c::Vector{Cdouble}, order, α)
    org_order = length(c)-1
    converted = zeros(order+1)
    ccall((:frqtr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          c, org_order, converted, order, α)
    converted
end

function mgc2mgc(c1::Vector{Cdouble}, α₁, γ₁, m2, α₂, γ₂)
    c2 = zeros(m2+1)
    m1 = length(c1)-1
    ccall((:mgc2mgc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Cdouble, Ptr{Cdouble}, Cint, Cdouble,
           Cdouble),
          c1, m1, α₁, γ₁, c2, m2, α₂, γ₂)
    c2
end

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

function mgclsp2sp(lsp::Vector{Cdouble}, α, γ, fftlen;
                   gain::Bool=true)
    sp = zeros(fftlen>>1 + 1)
    m = gain ? length(lsp)-1 : length(lsp)
    ccall((:mgclsp2sp, libSPTK), Void, (Cdouble, Cdouble, Ptr{Cdouble}, Cint,
                                        Ptr{Cdouble}, Cint, Cint),
          α, γ, lsp, m, sp, length(sp), convert(Int, gain))
    sp
end

function b2c(c::Vector{Cdouble}, order, α)
    org_order = length(c)-1
    converted = zeros(order+1)
    ccall((:b2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          c, org_order, converted, order, α)
    converted
end
