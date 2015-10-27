# Mel-generalized cepstrum analysis

function mcep!(mc::StridedVector{Cdouble}, windowed::StridedVector{Cdouble},
               α=0.35;
               miniter::Int=2,
               maxiter::Int=30,
               threshold::Float64=0.001,
               etype::Int=0,
               eps::Float64=0.0,
               min_det::Float64=1.0e-6,
               itype::Int=0)
    if itype ∉ 0:4
        throw(ArgumentError("unsupported itype: $itype, must be ∈ 0:4"))
    end
    if etype ∉ 0:2
        throw(ArgumentError("unsupported etype: $etype, must be ∈ 0:2"))
    end
    if etype == 0 && eps != 0.0
        throw(ArgumentError("eps cannot be specified for etype = 0"))
    end
    if etype ∈ 1:2 && eps < 0.0
        throw(ArgumentError("eps: $eps, must be >= 0"))
    end
    if min_det < 0.0
        throw(ArgumentError("min_det must be positive: min_det = $min_det"))
    end

    order = length(mc) - 1
    ret = ccall((:mcep, libSPTK), Cint,
                (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint,
                 Cdouble, Cint, Cint, Cdouble, Cint, Cdouble, Cdouble, Cint),
                windowed, length(windowed), mc, order, α,
                miniter, maxiter, threshold, etype, eps, min_det, itype)
    @assert ret ∈ -1:0 || ret ∈ 3:4
    if ret == 3
        error("failed to compute mcep; error occured in theq")
    elseif ret == 4
        error("zero(s) are found in periodogram, use eps option to floor")
    end

    mc
end

function mcep(windowed::StridedVector{Cdouble}, order=25, α=0.35; kargs...)
    mc = Array{Cdouble}(order+1)
    mcep!(mc, windowed, α; kargs...)
end

function gcep!(gc::StridedVector{Cdouble}, windowed::StridedVector{Cdouble},
               γ=0.0;
               miniter::Int=2,
               maxiter::Int=30,
               threshold::Float64=0.001,
               etype::Int=0,
               eps::Float64=0.0,
               min_det::Float64=1.0e-6,
               itype::Int=0,
               norm::Bool=false)
    assert_gamma(γ)
    if itype ∉ 0:4
        throw(ArgumentError("unsupported itype: $itype, must be ∈ 0:4"))
    end
    if etype ∉ 0:2
        throw(ArgumentError("unsupported etype: $etype, must be ∈ 0:2"))
    end
    if etype == 0 && eps != 0.0
        throw(ArgumentError("eps cannot be specified for etype = 0"))
    end
    if etype ∈ 1:2 && eps < 0.0
        throw(ArgumentError("eps: $eps, must be >= 0"))
    end
    if min_det < 0.0
        throw(ArgumentError("min_det must be positive: min_det = $min_det"))
    end

    order = length(gc) - 1
    # Note that the output is gain-normalized
    ret = ccall((:gcep, libSPTK), Cint,
                (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint,
                 Cdouble, Cint, Cint, Cdouble, Cint, Cdouble, Cdouble, Cint),
                windowed, length(windowed), gc, order, γ,
                miniter, maxiter, threshold, etype, eps, min_det, itype)
    @assert ret ∈ -1:0 || ret == 3
    if ret == 3
        error("failed to compute gcep; error occured in theq")
    end

    !norm && ignorm!(gc, gc, γ)

    gc
end

function gcep(windowed::StridedVector{Cdouble}, order=25, γ=0.0; kargs...)
    gc = Array{Cdouble}(order + 1)
    gcep!(gc, windowed, γ; kargs...)
end

function mgcep!(mgc::StridedVector{Cdouble}, windowed::StridedVector{Cdouble},
                α=0.35,
                γ=0.0;
                num_recursions::Int=length(windowed)-1,
                miniter::Int=2,
                maxiter::Int=30,
                threshold::Float64=0.001,
                etype::Int=0,
                eps::Float64=0.0,
                min_det::Float64=1.0e-6,
                itype::Int=0,
                otype::Int=0)
    assert_gamma(γ)
    if itype ∉ 0:4
        throw(ArgumentError("unsupported itype: $itype, must be ∈ 0:4"))
    end
    if etype ∉ 0:2
        throw(ArgumentError("unsupported etype: $etype, must be ∈ 0:2"))
    end
    if etype == 0 && eps != 0.0
        throw(ArgumentError("eps cannot be specified for etype = 0"))
    end
    if etype ∈ 1:2 && eps < 0.0
        throw(ArgumentError("eps: $eps, must be >= 0"))
    end
    if min_det < 0.0
        throw(ArgumentError("min_det must be positive: min_det = $min_det"))
    end
    if otype ∉ 0:5
        throw(ArgumentError("unsupported otype: $otype, must be ∈ 0:5"))
    end

    order = length(mgc) - 1
    ret = ccall((:mgcep, libSPTK), Cint,
                (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble, Cdouble,
                 Cint, Cint, Cint, Cdouble, Cint, Cdouble, Cdouble, Cint),
                windowed, length(windowed), mgc, order, α, γ, num_recursions,
                miniter, maxiter, threshold, etype, eps, min_det, itype)
    @assert ret ∈ -1:0 || ret == 3
    if ret == 3
        error("failed to compute mgcep; error occured in theq")
    end

    if otype ∈ 0:2 || otype == 4
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

    if otype ∈ 4:5
        for i=2:length(mgc)
            mgc[i] *= γ
        end
    end

    mgc
end

function mgcep(windowed::StridedVector{Cdouble}, order=25, α=0.35, γ=0.0;
               kargs...)
    mgc = Array{Cdouble}(order + 1)
    mgcep!(mgc, windowed, α, γ; kargs...)
end

function uels!(c::StridedVector{Cdouble}, windowed::StridedVector{Cdouble};
               miniter::Int=2,
               maxiter::Int=30,
               threshold::Float64=0.001,
               etype::Int=0,
               eps::Float64=0.0,
               itype::Int=0)
    if itype ∉ 0:4
        throw(ArgumentError("unsupported itype: $itype, must be ∈ 0:4"))
    end
    if etype ∉ 0:2
        throw(ArgumentError("unsupported etype: $etype, must be ∈ 0:2"))
    end
    if etype == 0 && eps != 0.0
        throw(ArgumentError("eps cannot be specified for etype = 0"))
    end
    if etype ∈ 1:2 && eps < 0.0
        throw(ArgumentError("eps: $eps, must be >= 0"))
    end

    order = length(c) - 1
    ret = ccall((:uels, libSPTK), Cint,
                (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble, Cint,
                 Cdouble, Cint),
                windowed, length(windowed), c, order,
                miniter, maxiter, threshold, etype, eps, itype)
    @assert ret ∈ -1:0 || ret == 3
    if ret == 3
        error("zero(s) are found in periodogram, use eps option to floor")
    end

    c
end

function uels(windowed::StridedVector{Cdouble}, order=25; kargs...)
    c = Array{Cdouble}(order + 1)
    uels!(c, windowed; kargs...)
end

function fftcep!(c::StridedVector{Cdouble}, logsp::StridedVector{Cdouble};
                 num_iter::Int=0,
                 acceleration_factor::Float64=0.0)
    order = length(c) - 1
    ccall((:fftcep, libSPTK), Void,
          (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cint, Cdouble),
          logsp, length(logsp), c, order, num_iter, acceleration_factor)
    c
end

function fftcep(logsp::StridedVector{Cdouble}, order=25; kargs...)
    c = Array{Cdouble}(order + 1)
    fftcep!(c, logsp; kargs...)
end

function lpc!(a::StridedVector{Cdouble}, x::StridedVector{Cdouble};
              min_det::Float64=1.0e-6)
    if min_det < 0.0
        throw(ArgumentError("min_det must be positive: min_det = $min_det"))
    end

    order = length(a) - 1
    ret = ccall((:lpc, libSPTK), Cint,
                (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Cdouble),
                x, length(x), a, order, min_det)
    @assert ret ∈ -2:0
    if ret == -2
        warn("failed to compute `stable` LPC. Please try again with different parameters")
    elseif ret == -1
        error("failed to compute LPC. Please try again with different parameters")
    end

    a
end

function lpc(x::StridedVector{Cdouble}, order=25; kargs...)
    a = Array{Cdouble}(order+1)
    lpc!(a, x; kargs...)
end
