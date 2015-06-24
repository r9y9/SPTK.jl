# Speech waveform synthesis filters

poledf_delay_length(order) = order
poledf_delay(order) = zeros(poledf_delay_length(order))

function poledf(x::Cdouble, a::Vector{Cdouble}, delay::Vector{Cdouble})
    order = length(a) - 1
    if length(delay) != poledf_delay_length(order)
        throw(DimensionMismatch("inconsistent delay length"))
    end
    ccall((:poledf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, a, length(a)-1, delay)
end

lmadf_delay_length(order, pd) = 2pd*(order + 1)
lmadf_delay(order, pd) = zeros(lmadf_delay_length(order, pd))

function lmadf(x::Cdouble, b::Vector{Cdouble}, pd, delay::Vector{Cdouble})
    assert_pade(pd)
    order = length(b) - 1
    if length(delay) != lmadf_delay_length(order, pd)
        throw(DimensionMismatch("inconsistent delay length"))
    end
    ccall((:lmadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, pd, delay)
end

lspdf_delay_length(order) = 2order + 1
lspdf_delay(order) = zeros(lspdf_delay_length(order))

function lspdf(x::Cdouble, f::Vector{Cdouble}, delay::Vector{Cdouble})
    order = length(f) - 1
    if iseven(order)
        lspdf_even(x, f, delay)
    else
        lspdf_odd(x, f, delay)
    end
end

function lspdf_even(x::Cdouble, f::Vector{Cdouble}, delay::Vector{Cdouble})
    order = length(f) - 1
    if length(delay) != lspdf_delay_length(order)
        throw(DimensionMismatch("inconsistent delay length"))
    end
    ccall((:lspdf_even, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, f, length(f) - 1, delay)
end

function lspdf_odd(x::Cdouble, f::Vector{Cdouble}, delay::Vector{Cdouble})
    order = length(f) - 1
    if length(delay) != lspdf_delay_length(order)
        throw(DimensionMismatch("inconsistent delay length"))
    end
    ccall((:lspdf_odd, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, f, length(f) - 1, delay)
end

ltcdf_delay_length(order) = order + 1
ltcdf_delay(order) = zeros(ltcdf_delay_length(order))

function ltcdf(x::Cdouble, k::Vector{Cdouble}, delay::Vector{Cdouble})
    order = length(k) - 1
    if length(delay) != ltcdf_delay_length(order)
        throw(DimensionMismatch("invalid delay length"))
    end
    ccall((:ltcdf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
          x, k, length(k) - 1, delay)
end

glsadf_delay_length(order, stage) = order * (stage + 1) + 1
glsadf_delay(order, stage) = zeros(glsadf_delay_length(order, stage))

# NOTE: stage = -1/γ
function glsadf(x::Cdouble, c::Vector{Cdouble}, stage, delay::Vector{Cdouble})
    stage >= 1 || throw(ArgumentError("stage >= 1 (-1 <= γ < 0)"))
    order = length(c) - 1
    if length(delay) != glsadf_delay_length(order, stage)
        throw(DimensionMismatch("invalid delay length"))
    end
    ccall((:glsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
          x, c, length(c) - 1, stage, delay)
end

# see mlsadf.c in original SPTK for this magic allocation
mlsadf_delay_length(order, pd) = 3*(pd+1) + pd*(order+2)
mlsadf_delay(order, pd) = zeros(mlsadf_delay_length(order, pd))

function mlsadf(x::Cdouble, b::Vector{Cdouble}, α, pd, delay::Vector{Cdouble})
    assert_pade(pd)
    order = length(b) - 1
    if length(delay) != mlsadf_delay_length(order, pd)
        throw(DimensionMismatch("invalid delay length"))
    end
    ccall((:mlsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, α, pd, delay)
end

# see mglsadf.c in original SPTK for this magic allocation
mglsadf_delay_length(order, stage) = (order+1)*stage
mglsadf_delay(order, stage) = zeros(mglsadf_delay_length(order, stage))

function mglsadf(x::Cdouble, b::Vector{Cdouble}, α, stage,
                 delay::Vector{Cdouble})
    stage >= 1 || throw(ArgumentError("stage >= 1 (-1 <= γ < 0)"))
    order = length(b) - 1
    if length(delay) != mglsadf_delay_length(order, stage)
        throw(DimensionMismatch("invalid delay length"))
    end
    ccall((:mglsadf, libSPTK), Cdouble,
          (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cint, Ptr{Cdouble}),
          x, b, length(b)-1, α, stage, delay)
end
