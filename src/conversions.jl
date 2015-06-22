# Conversions

## LPC, LSP, PARCOR conversions

function lpc2c!(dst_ceps::Vector{Cdouble}, src_lpc::Vector{Cdouble})
    src_order = length(src_lpc) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:lpc2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          src_lpc, src_order, dst_ceps, dst_order)
    dst_ceps
end

function lpc2c(src_lpc::Vector{Cdouble}, dst_order=length(src_lpc)-1)
    src_order = length(src_lpc) - 1
    dst_ceps = Array(Cdouble, dst_order+1)
    lpc2c!(dst_ceps, src_lpc)
end

function lpc2lsp!(lsp::Vector{Cdouble}, lpc::Vector{Cdouble};
                  numsp::Int=128,
                  maxiter::Int=4,
                  eps::Float64=1e-6,
                  loggain::Bool=true,
                  otype::Int=0,
                  fs=nothing)
    dst_order = length(lsp) - 1
    ccall((:lpc2lsp, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble),
          lpc, sub(lsp, 2:length(lsp)), dst_order, numsp, maxiter, eps)

    if otype == 0
        for i=2:length(lsp)
            @inbounds lsp[i] *= 2π
        end
    elseif otype ∈ 2:3
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

function lpc2lsp(lpc::Vector{Cdouble}, dst_order; kargs...)
    lsp = zeros(Cdouble, dst_order+1)
    lpc2lsp!(lsp, lpc; kargs...)
end

function lpc2par!(par::Vector{Cdouble}, lpc::Vector{Cdouble})
    if length(par) != length(lpc)
        throw(DimensionMismatch("inconsistent dimentions"))
    end
    ccall((:lpc2par, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
          lpc, par, length(lpc)-1)
    par
end

function lpc2par(lpc::Vector{Cdouble})
    par = similar(lpc)
    lpc2par!(par, lpc)
end

# assume lsp has loggain at lsp[1]
function lsp2sp!(sp::Vector{Cdouble}, lsp::Vector{Cdouble})
    fftlen = (length(sp) - 1)<<1
    assert_fftlen(fftlen)
    ccall((:lsp2sp, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint),
          lsp, length(lsp)-1, sp, length(sp), 1)
    sp
end

function lsp2sp(lsp::Vector{Cdouble}, fftlen)
    assert_fftlen(fftlen)
    sp = Array(Cdouble, fftlen>>1+1)
    lsp2sp!(sp, lsp)
end

## Mel-generalized cepstrum conversions

function mc2b!(b::Vector{Cdouble}, mc::Vector{Cdouble}, α=0.41)
    if length(b) != length(mc)
        throw(DimensionMismatch("inconstent dimensions"))
    end
    order = length(b) - 1
    ccall((:mc2b, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          mc, b, order, α)
    b
end

function mc2b(mc::Vector{Cdouble}, α=0.41)
    b = similar(mc)
    mc2b!(b, mc, α)
end

function b2mc!(mc::Vector{Cdouble}, b::Vector{Cdouble}, α=0.41)
    if length(mc) != length(b)
        throw(DimensionMismatch("inconstent dimensions"))
    end
    order = length(b) - 1
    ccall((:b2mc, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          b, mc, order, α)
    mc
end

function b2mc(b::Vector{Cdouble}, α=0.41)
    mc = similar(b)
    b2mc!(mc, b, α)
end

function b2c!(dst_ceps::Vector{Cdouble}, src_b::Vector{Cdouble}, α)
    src_order = length(src_b) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:b2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          src_b, src_order, dst_ceps, dst_order, α)
    dst_ceps
end

function b2c(src_b::Vector{Cdouble}, dst_order, α)
    dst_ceps = Array(Cdouble, dst_order + 1)
    b2c!(dst_ceps, src_b, α)
end

function c2acr!(r::Vector{Cdouble}, c::Vector{Cdouble}, fftlen=256)
    assert_fftlen(fftlen)
    designed_order = length(r) - 1
    ccall((:c2acr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint),
          c, length(c) - 1, r, designed_order, fftlen)
    r
end

function c2acr(c::Vector{Cdouble}, order=length(c)-1, fftlen=256)
    r = Array(Cdouble, order + 1)
    c2acr!(r, c, fftlen)
end

function c2ir!(h::Vector{Cdouble}, c::Vector{Cdouble})
    order = length(c) # NOT length(c) - 1
    ccall((:c2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, order, h, length(h))
    h
end

function c2ir(c::Vector{Cdouble}, len)
    h = Array(Cdouble, len)
    c2ir!(h, c)
end

function ic2ir!(c::Vector{Cdouble}, h::Vector{Cdouble})
    ccall((:ic2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          h, length(h), c, length(c))
    c
end

function ic2ir(h::Vector{Cdouble}, order)
    c = Array(Cdouble, order+1)
    ic2ir!(c, h)
end

function c2ndps!(ndps::Vector{Cdouble}, c::Vector{Cdouble})
    fftlen = (length(ndps) - 1)<<1
    assert_fftlen(fftlen)
    buf = zeros(fftlen)
    m = length(c)-1
    ccall((:c2ndps, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, m, buf, fftlen)
    for i = 1:length(ndps)
        @inbounds ndps[i] = buf[i]
    end

    ndps
end

function c2ndps(c::Vector{Cdouble}, fftlen)
    assert_fftlen(fftlen)
    ndps = Array(Cdouble, fftlen>>1 + 1)
    c2ndps!(ndps, c)
end

function ndps2c!(dst_ceps::Vector{Cdouble}, ndps::Vector{Cdouble})
    order = length(dst_ceps) - 1
    fftlen = (length(ndps) - 1)<<1 # assuming the length of npds is fftsize/2+1
    assert_fftlen(fftlen)
    ccall((:ndps2c, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          ndps, fftlen, dst_ceps, order)
    dst_ceps
end

function ndps2c(ndps::Vector{Cdouble}, order)
    dst_ceps = Array(Cdouble, order + 1)
    ndps2c!(dst_ceps, ndps)
end

function gc2gc!(dst_ceps::Vector{Cdouble}, dst_γ,
                src_ceps::Vector{Cdouble}, src_γ)
    dst_order = length(dst_ceps) - 1
    src_order = length(src_ceps) - 1
    ccall((:gc2gc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, src_order, src_γ, dst_ceps, dst_order, dst_γ)
    dst_ceps
end

function gc2gc(src_ceps::Vector{Cdouble}, src_γ, dst_order, dst_γ)
    src_order = length(src_ceps) - 1
    dst_ceps = Array(Cdouble, dst_order + 1)
    gc2gc!(dst_ceps, dst_γ, src_ceps, src_γ)
end

function gnorm!(c::Vector{Cdouble}, γ)
    m = length(c) - 1
    ccall((:gnorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), c, c, m, γ)
    c
end

function gnorm(c::Vector{Cdouble}, γ)
    normalizedC = similar(c)
    m = length(c) - 1
    ccall((:gnorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), c, normalizedC, m, γ)
    normalizedC
end

function ignorm!(c::Vector{Cdouble}, γ)
    m = length(c) - 1
    ccall((:ignorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), c, c, m, γ)
    c
end

function ignorm(normalizedC::Vector{Cdouble}, γ)
    c = similar(normalizedC)
    m = length(c) - 1
    ccall((:ignorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble), normalizedC, c, m, γ)
    c
end

function freqt!(dst_ceps::Vector{Cdouble}, src_ceps::Vector{Cdouble}, α)
    src_order = length(src_ceps) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:freqt, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, src_order, dst_ceps, dst_order, α)
    dst_ceps
end

function freqt(ceps::Vector{Cdouble}, order, α)
    dst_ceps = Array(Cdouble, order + 1)
    freqt!(dst_ceps, ceps, α)
end

function frqtr!(dst_ceps::Vector{Cdouble}, src_ceps::Vector{Cdouble}, α)
    src_order = length(src_ceps) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:frqtr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, src_order, dst_ceps, dst_order, α)
    dst_ceps
end

function frqtr(c::Vector{Cdouble}, order, α)
    dst_ceps = Array(Cdouble, order + 1)
    frqtr!(dst_ceps, c, α)
end

function mgc2mgc!(dst_ceps::Vector{Cdouble}, dst_α, dst_γ,
                  src_ceps::Vector{Cdouble}, src_α, src_γ)
    src_order = length(src_ceps) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:mgc2mgc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Cdouble, Ptr{Cdouble}, Cint, Cdouble,
           Cdouble),
          src_ceps, src_order, src_α, src_γ,
          dst_ceps, dst_order, dst_α, dst_γ)
    dst_ceps
end

function mgc2mgc(src_ceps::Vector{Cdouble}, src_α, src_γ,
                 dst_order, dst_α, dst_γ)
    dst_ceps = Array(Cdouble, dst_order + 1)
    src_order = length(src_ceps) - 1
    mgc2mgc!(dst_ceps, dst_α, dst_γ, src_ceps, src_α, src_γ)
end

function mgc2sp!(sp::Vector{Complex{Cdouble}},
                 mgc::Vector{Cdouble}, α, γ)
    fftlen = (length(sp)-1)<<1
    assert_fftlen(fftlen)

    order = length(mgc) - 1
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

function mgc2sp(mgc::Vector{Cdouble}, α, γ, fftlen)
    assert_fftlen(fftlen)
    order = length(mgc) - 1
    sp = Array(Complex{Cdouble}, fftlen>>1 + 1)
    mgc2sp!(sp, mgc, α, γ)
end

function mgclsp2sp!(sp::Vector{Cdouble}, lsp::Vector{Cdouble}, α, γ;
                    gain::Bool=true)
    fftlen = (length(sp) - 1)<<1
    assert_fftlen(fftlen)
    order = gain ? length(lsp)-1 : length(lsp)
    ccall((:mgclsp2sp, libSPTK), Void, (Cdouble, Cdouble, Ptr{Cdouble}, Cint,
                                        Ptr{Cdouble}, Cint, Cint),
          α, γ, lsp, order, sp, length(sp), convert(Cint, gain))
    sp
end

function mgclsp2sp(lsp::Vector{Cdouble}, α, γ, fftlen; kargs...)
    assert_fftlen(fftlen)
    sp = Array(Cdouble, fftlen>>1 + 1)
    mgclsp2sp!(sp, lsp, α, γ; kargs...)
end
