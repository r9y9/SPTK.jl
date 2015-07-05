# Adaptive mel-cepstrum analysis

function acep!(c::StridedVector{Cdouble}, x::Cdouble;
               λ::Float64=0.98,
               step::Float64=0.1,
               τ::Float64=0.9,
               pd::Int=4,
               eps::Float64=1.0e-6)
    assert_pade(pd)
    order = length(c) - 1
    prederr = ccall((:acep, libSPTK), Cdouble,
                    (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cdouble, Cdouble,
                     Cint, Cdouble),
                    x, c, order, λ, step, τ, pd, eps)
    prederr
end

# γ = -1/stage (-1 <= γ < 0)
function agcep!(c::StridedVector{Cdouble}, x::Cdouble, stage=1;
                λ::Float64=0.98,
                step::Float64=0.1,
                τ::Float64=0.9,
                eps::Float64=1.0e-6)
    stage >= 1 || throw(ArgumentError("stage >= 1 (-1 <= γ < 0)"))
    order = length(c) - 1
    prederr = ccall((:agcep, libSPTK), Cdouble,
                    (Cdouble, Ptr{Cdouble}, Cint, Cint, Cdouble, Cdouble,
                     Cdouble, Cdouble),
                    x, c, order, stage, λ, step, τ, eps)
    prederr
end

function amcep!(b::StridedVector{Cdouble}, x::Cdouble, α=0.41;
                λ::Float64=0.98,
                step::Float64=0.1,
                τ::Float64=0.9,
                pd::Int=4,
                eps::Float64=1.0e-6)
    assert_pade(pd)
    order = length(b) - 1
    prederr = ccall((:amcep, libSPTK), Cdouble,
                    (Cdouble, Ptr{Cdouble}, Cint, Cdouble, Cdouble, Cdouble,
                     Cdouble, Cint, Cdouble),
                    x, b, order, α, λ, step, τ, pd, eps)
    prederr
end
