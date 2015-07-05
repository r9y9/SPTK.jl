# Library routines

function agexp(r, x, y)
    ccall((:agexp, libSPTK), Cdouble, (Cdouble, Cdouble, Cdouble), r, x, y)
end

gexp(r, x) = ccall((:gexp, libSPTK), Cdouble, (Cdouble, Cdouble), r, x)
glog(r, x) = ccall((:glog, libSPTK), Cdouble, (Cdouble, Cdouble), r, x)

mseq() = ccall((:mseq, libSPTK), Cint, ())

# solve (T + H)a = b
# NOTE: the order of arguments is different in theq in SPTK
# I would put the solution vector `a` as the first argument instead of third.
function theq!(a::StridedVector{Cdouble}, t::StridedVector{Cdouble},
               h::StridedVector{Cdouble}, b::StridedVector{Cdouble};
               min_det::Float64=1.0e-6)
    n = length(a)
    if length(b) != n || length(t) != n || length(h) != 2n-1
        throw(DimensionMismatch("inconsistent vector length"))
    end
    if min_det < 0.0
        throw(ArgumentError("min_det must be positive: min_det = $min_det"))
    end

    ret = ccall((:theq, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint,
                 Cdouble), t, h, a, b, n, min_det)
    if ret != 0
        error("failed to solve Toeplitz plus Hankel matrix system (T + H)a = b")
    end

    a
end

function theq(t::StridedVector{Cdouble}, h::StridedVector{Cdouble},
              b::StridedVector{Cdouble};
              kargs...)
    a = similar(t)
    theq!(a, t, h, b, kargs...)
end

# solve Ta = b
# NOTE: the order of arguments is different in toeplitz in SPTK
# I would put the solution vector `a` as the first argument instead of second.
function toeplitz!(a::StridedVector{Cdouble}, t::StridedVector{Cdouble},
                   b::StridedVector{Cdouble};
                   eps::Float64=1.0e-6)
    if length(a) != length(t) || length(a) != length(b)
        throw(DimensionMismatch("inconsistent vector length"))
    end
    if eps < 0.0
        throw(ArgumentError("eps must be positive: eps = $eps"))
    end

    ret = ccall((:toeplitz, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
                t, a, b, length(t), eps)
    if ret != 0
        error("failed to solve a symmetric Toeplitz set of linear equations")
    end

    a
end

function toeplitz(t::StridedVector{Cdouble}, b::StridedVector{Cdouble};
                  kargs...)
    a = similar(t)
    toeplitz!(a, t, b, kargs...)
end
