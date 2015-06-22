# Speech waveform synthesis filters

function poledf(x::Cdouble, a::Vector{Cdouble}, delay::Vector{Cdouble})
    ccall((:poledf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, a, length(a)-1, delay)
end

poledf_delay(order) = zeros(order)

function lmadf(x::Cdouble, b::Vector{Cdouble}, pd, delay::Vector{Cdouble})
    assert_pade(pd)
    ccall((:lmadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, pd, delay)
end

lmadf_delay(order, pd) = zeros(2pd*(order+1))

function lspdf(x::Cdouble, f::Vector{Cdouble}, delay::Vector{Cdouble})
    order = length(f) - 1
    if iseven(order)
        lspdf_even(x, f, delay)
    else
        lspdf_odd(x, f, delay)
    end
end

lspdf_delay(order) = zeros(2order + 1)

function lspdf_even(x::Cdouble, f::Vector{Cdouble}, delay::Vector{Cdouble})
    ccall((:lspdf_even, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, f, length(f) - 1, delay)
end

function lspdf_odd(x::Cdouble, f::Vector{Cdouble}, delay::Vector{Cdouble})
    ccall((:lspdf_odd, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, f, length(f) - 1, delay)
end

function ltcdf(x::Cdouble, k::Vector{Cdouble}, delay::Vector{Cdouble})
    if length(k) != length(delay)
        throw(ArgumentError("invalid delay length"))
    end
    ccall((:ltcdf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, k, length(k) - 1, delay)
end

ltcdf_delay(order) = zeros(order + 1)

# NOTE: stage = -1/γ
function glsadf(x::Cdouble, c::Vector{Cdouble}, stage, delay::Vector{Cdouble})
    stage >= 1 || throw(ArgumentError("stage >= 1 (-1 <= γ < 0)"))
    ccall((:glsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
          x, c, length(c) - 1, stage, delay)
end

glsadf_delay(order, stage) = zeros(order * (stage + 1) + 1)

function mlsadf(x::Cdouble, b::Vector{Cdouble}, α, pd, delay::Vector{Cdouble})
    assert_pade(pd)
    ccall((:mlsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, α, pd, delay)
end

# see mlsadf.c in original SPTK for this magic allocation
mlsadf_delay(order, pd) = zeros(3*(pd+1) + pd*(order+2))

function mglsadf(x::Cdouble, b::Vector{Cdouble}, α, stage,
                 delay::Vector{Cdouble})
    stage >= 1 || throw(ArgumentError("stage >= 1 (-1 <= γ < 0)"))
    ccall((:mglsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, α, stage, delay)
end

# see mglsadf.c in original SPTK for this magic allocation
mglsadf_delay(order, stage) = zeros((order+1)*stage)
