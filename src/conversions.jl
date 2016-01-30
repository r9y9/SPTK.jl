# Conversions

## LPC, LSP, PARCOR conversions

function lpc2c!(dst_ceps::StridedVector{Cdouble},
                src_lpc::StridedVector{Cdouble})
    src_order = length(src_lpc) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:lpc2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          src_lpc, src_order, dst_ceps, dst_order)
    dst_ceps
end

function lpc2c(src_lpc::StridedVector{Cdouble}, dst_order=length(src_lpc)-1)
    src_order = length(src_lpc) - 1
    dst_ceps = Array{Cdouble}(dst_order+1)
    lpc2c!(dst_ceps, src_lpc)
end

function lpc2lsp!(lsp::StridedVector{Cdouble}, lpc::StridedVector{Cdouble};
                  numsp::Int=512,
                  maxiter::Int=4,
                  eps::Float64=1e-6,
                  loggain::Bool=false,
                  otype::Int=0,
                  fs=nothing)
    if length(lsp) != length(lpc)
        throw(DimensionMismatch("inconsistent dimentions"))
    end
    fill!(lsp, zero(eltype(lsp)))
    order = length(lsp) - 1
    ccall((:lpc2lsp, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble),
          lpc, pointer(lsp, 2), order, numsp, maxiter, eps)

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

function lpc2lsp(lpc::StridedVector{Cdouble}; kargs...)
    lsp = zeros(lpc)
    lpc2lsp!(lsp, lpc; kargs...)
end

function lpc2par!(par::StridedVector{Cdouble}, lpc::StridedVector{Cdouble})
    if length(par) != length(lpc)
        throw(DimensionMismatch("inconsistent dimentions"))
    end
    ccall((:lpc2par, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
          lpc, par, length(lpc)-1)
    par
end

function lpc2par(lpc::StridedVector{Cdouble})
    par = similar(lpc)
    lpc2par!(par, lpc)
end

function par2lpc!(lpc::StridedVector{Cdouble}, par::StridedVector{Cdouble})
    if length(lpc) != length(par)
        throw(DimensionMismatch("inconsistent dimentions"))
    end
    # NOTE: Since par2lpc destroys input PARCOR coefficients,
    # so passing its copy to the ccall to keep PARCOR coeff. unchanged
    ccall((:par2lpc, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
          copy(par), lpc, length(par)-1)
    lpc
end

function par2lpc(par::StridedVector{Cdouble})
    lpc = similar(par)
    par2lpc!(lpc, par)
    return lpc
end

# assume lsp has loggain at lsp[1]
function lsp2sp!(sp::StridedVector{Cdouble}, lsp::StridedVector{Cdouble})
    fftlen = (length(sp) - 1)<<1
    assert_fftlen(fftlen)
    ccall((:lsp2sp, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint),
          lsp, length(lsp)-1, sp, length(sp), 1)
    sp
end

function lsp2sp(lsp::StridedVector{Cdouble}, fftlen=256)
    assert_fftlen(fftlen)
    sp = Array{Cdouble}(fftlen>>1+1)
    lsp2sp!(sp, lsp)
end

## Mel-generalized cepstrum conversions

function mc2b!(b::StridedVector{Cdouble}, mc::StridedVector{Cdouble}, α=0.35)
    if length(b) != length(mc)
        throw(DimensionMismatch("inconstent dimensions"))
    end
    order = length(b) - 1
    ccall((:mc2b, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          mc, b, order, α)
    b
end

function mc2b(mc::StridedVector{Cdouble}, α=0.35)
    b = similar(mc)
    mc2b!(b, mc, α)
end

function b2mc!(mc::StridedVector{Cdouble}, b::StridedVector{Cdouble}, α=0.35)
    if length(mc) != length(b)
        throw(DimensionMismatch("inconstent dimensions"))
    end
    order = length(b) - 1
    ccall((:b2mc, libSPTK), Void, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          b, mc, order, α)
    mc
end

function b2mc(b::StridedVector{Cdouble}, α=0.35)
    mc = similar(b)
    b2mc!(mc, b, α)
end

function b2c!(dst_ceps::StridedVector{Cdouble}, src_b::StridedVector{Cdouble},
              α)
    src_order = length(src_b) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:b2c, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          src_b, src_order, dst_ceps, dst_order, α)
    dst_ceps
end

function b2c(src_b::StridedVector{Cdouble}, dst_order=length(src_b)-1, α=0.35)
    dst_ceps = Array{Cdouble}(dst_order + 1)
    b2c!(dst_ceps, src_b, α)
end

function c2acr!(r::StridedVector{Cdouble}, c::StridedVector{Cdouble},
                fftlen=256)
    assert_fftlen(fftlen)
    dst_order = length(r) - 1
    ccall((:c2acr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint),
          c, length(c) - 1, r, dst_order, fftlen)
    r
end

function c2acr(c::StridedVector{Cdouble}, order=length(c)-1, fftlen=256)
    r = Array{Cdouble}(order + 1)
    c2acr!(r, c, fftlen)
end

function c2ir!(h::StridedVector{Cdouble}, c::StridedVector{Cdouble})
    order = length(c) # NOT length(c) - 1
    ccall((:c2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, order, h, length(h))
    h
end

function c2ir(c::StridedVector{Cdouble}, len=256)
    h = Array{Cdouble}(len)
    c2ir!(h, c)
end

function ic2ir!(c::StridedVector{Cdouble}, h::StridedVector{Cdouble})
    ccall((:ic2ir, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          h, length(h), c, length(c))
    c
end

function ic2ir(h::StridedVector{Cdouble}, order=25)
    c = Array{Cdouble}(order+1)
    ic2ir!(c, h)
end

function c2ndps!(ndps::StridedVector{Cdouble}, c::StridedVector{Cdouble})
    fftlen = (length(ndps) - 1)<<1
    assert_fftlen(fftlen)
    buf = Array{Cdouble}(fftlen)
    m = length(c)-1
    ccall((:c2ndps, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          c, m, buf, fftlen)
    for i = 1:length(ndps)
        @inbounds ndps[i] = buf[i]
    end

    ndps
end

function c2ndps(c::StridedVector{Cdouble}, fftlen=256)
    assert_fftlen(fftlen)
    ndps = Array{Cdouble}(fftlen>>1 + 1)
    c2ndps!(ndps, c)
end

function ndps2c!(dst_ceps::StridedVector{Cdouble}, ndps::StridedVector{Cdouble})
    order = length(dst_ceps) - 1
    fftlen = (length(ndps) - 1)<<1 # assuming the length of npds is fftsize/2+1
    assert_fftlen(fftlen)
    ccall((:ndps2c, libSPTK), Void, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint),
          ndps, fftlen, dst_ceps, order)
    dst_ceps
end

function ndps2c(ndps::StridedVector{Cdouble}, order=25)
    dst_ceps = Array{Cdouble}(order + 1)
    ndps2c!(dst_ceps, ndps)
end

function gc2gc!(dst_ceps::StridedVector{Cdouble}, dst_γ,
                src_ceps::StridedVector{Cdouble}, src_γ)
    assert_gamma(dst_γ)
    assert_gamma(src_γ)
    dst_order = length(dst_ceps) - 1
    src_order = length(src_ceps) - 1
    ccall((:gc2gc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, src_order, src_γ, dst_ceps, dst_order, dst_γ)
    dst_ceps
end

function gc2gc(src_ceps::StridedVector{Cdouble}, src_γ=0.0,
               dst_order=length(src_ceps)-1, dst_γ=0.0)
    src_order = length(src_ceps) - 1
    dst_ceps = Array{Cdouble}(dst_order + 1)
    gc2gc!(dst_ceps, dst_γ, src_ceps, src_γ)
end

function gc2gc(src_ceps::StridedVector{Cdouble}, src_γ::AbstractFloat,
        dst_γ::AbstractFloat)
    gc2gc(src_ceps, src_γ, length(src_ceps)-1, dst_γ)
end

function gnorm!(dst_ceps::StridedVector{Cdouble},
                src_ceps::StridedVector{Cdouble}, γ=0.0)
    assert_gamma(γ)
    if length(dst_ceps) != length(src_ceps)
        throw(DimensionMismatch("inconsistent dimensions"))
    end
    order = length(dst_ceps) - 1
    ccall((:gnorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, dst_ceps, order, γ)
    dst_ceps
end

function gnorm(src_ceps::StridedVector{Cdouble}, γ=0.0)
    dst_ceps = similar(src_ceps)
    gnorm!(dst_ceps, src_ceps, γ)
end

function ignorm!(dst_ceps::StridedVector{Cdouble},
                 src_ceps::StridedVector{Cdouble}, γ=0.0)
    assert_gamma(γ)
    if length(dst_ceps) != length(src_ceps)
        throw(DimensionMismatch("inconsistent dimensions"))
    end
    order = length(dst_ceps) - 1
    ccall((:ignorm, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, dst_ceps, order, γ)
    dst_ceps
end

function ignorm(src_ceps::StridedVector{Cdouble}, γ=0.0)
    dst_ceps = similar(src_ceps)
    ignorm!(dst_ceps, src_ceps, γ)
end

function freqt!(dst_ceps::StridedVector{Cdouble},
                src_ceps::StridedVector{Cdouble}, α=0.0)
    src_order = length(src_ceps) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:freqt, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, src_order, dst_ceps, dst_order, α)
    dst_ceps
end

function freqt(ceps::StridedVector{Cdouble}, order=25, α=0.0)
    dst_ceps = Array{Cdouble}(order + 1)
    freqt!(dst_ceps, ceps, α)
end

function freqt(ceps::StridedVector{Cdouble}, α::AbstractFloat)
    freqt(ceps, length(ceps)-1, α)
end

function frqtr!(dst_ceps::StridedVector{Cdouble},
                src_ceps::StridedVector{Cdouble}, α=0.0)
    src_order = length(src_ceps) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:frqtr, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          src_ceps, src_order, dst_ceps, dst_order, α)
    dst_ceps
end

function frqtr(ceps::StridedVector{Cdouble}, order=25, α=0.0)
    dst_ceps = Array{Cdouble}(order + 1)
    frqtr!(dst_ceps, ceps, α)
end

function frqtr(ceps::StridedVector{Cdouble}, α::AbstractFloat)
    frqtr(ceps, length(ceps)-1, α)
end

function mgc2mgc!(dst_ceps::StridedVector{Cdouble}, dst_α, dst_γ,
                  src_ceps::StridedVector{Cdouble}, src_α, src_γ)
    assert_gamma(dst_γ)
    assert_gamma(src_γ)
    src_order = length(src_ceps) - 1
    dst_order = length(dst_ceps) - 1
    ccall((:mgc2mgc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Cdouble, Cdouble, Ptr{Cdouble}, Cint, Cdouble,
           Cdouble),
          src_ceps, src_order, src_α, src_γ,
          dst_ceps, dst_order, dst_α, dst_γ)
    dst_ceps
end

function mgc2mgc(src_ceps::StridedVector{Cdouble}, src_α=0.0, src_γ=0.0,
                 dst_order=length(src_ceps)-1, dst_α=0.0, dst_γ=0.0)
    dst_ceps = Array{Cdouble}(dst_order + 1)
    src_order = length(src_ceps) - 1
    mgc2mgc!(dst_ceps, dst_α, dst_γ, src_ceps, src_α, src_γ)
end

function mgc2mgc(src_ceps::StridedVector{Cdouble}, src_α::AbstractFloat,
        src_γ::AbstractFloat,
        dst_α::AbstractFloat, dst_γ::AbstractFloat)
    mgc2mgc(src_ceps, src_α, src_γ, length(src_ceps)-1, dst_α, dst_γ)
end

function mgc2sp!(sp::Vector{Complex{Cdouble}},
                 mgc::StridedVector{Cdouble}, α=0.0, γ=0.0)
    assert_gamma(γ)
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

function mgc2sp(mgc::StridedVector{Cdouble}, α=0.0, γ=0.0, fftlen=256)
    assert_fftlen(fftlen)
    order = length(mgc) - 1
    sp = Array{Complex{Cdouble}}(fftlen>>1 + 1)
    mgc2sp!(sp, mgc, α, γ)
end

function mgclsp2sp!(sp::StridedVector{Cdouble},
                    lsp::StridedVector{Cdouble}, α=0.0, γ=0.0;
                    gain::Bool=true)
    assert_gamma(γ)
    fftlen = (length(sp) - 1)<<1
    assert_fftlen(fftlen)
    order = gain ? length(lsp)-1 : length(lsp)
    ccall((:mgclsp2sp, libSPTK), Void, (Cdouble, Cdouble, Ptr{Cdouble}, Cint,
                                        Ptr{Cdouble}, Cint, Cint),
          α, γ, lsp, order, sp, length(sp), convert(Cint, gain))
    sp
end

function mgclsp2sp(lsp::StridedVector{Cdouble}, α=0.0, γ=0.0, fftlen=256;
                   kargs...)
    assert_fftlen(fftlen)
    sp = Array{Cdouble}(fftlen>>1 + 1)
    mgclsp2sp!(sp, lsp, α, γ; kargs...)
end
