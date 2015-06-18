# Mel-generalized cepstrum analysis

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

function fftcep(logsp::Vector{Cdouble}, order;
                itr::Int=0, accelerationf::Float64=0.0)
    c = zeros(order+1)

    ccall((:fftcep, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Cdouble),
          logsp, length(logsp), c, length(c), itr, accelerationf)
    c
end

function lpc(x::Vector{Cdouble}, order;
             f::Float64=1e-6)
    a = Array(Cdouble, order+1)
    ccall((:lpc, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
          x, length(x), a, order, f)
    a
end
