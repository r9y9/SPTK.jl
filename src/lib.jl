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
function theq!(a::Vector{Cdouble}, t::Vector{Cdouble}, h::Vector{Cdouble},
               b::Vector{Cdouble}, eps=1.0e-6)
    n = length(a)
    if length(b) != n || length(t) != n || length(h) != 2n-1
        throw(DimensionMismatch("inconsistent vector length"))
    end

    ret = ccall((:theq, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint,
                 Cdouble), t, h, a, b, n, eps)
    if ret != 0
        error("failed to solve Toeplitz plus Hankel matrix system (T + H)a = b")
    end

    a
end

function theq(t::Vector{Cdouble}, h::Vector{Cdouble}, b::Vector{Cdouble};
              eps::Float64=1.0e-6)
    a = similar(t)
    theq!(a, t, h, b, eps)
end

function toeplitz!(t::Vector{Cdouble}, a::Vector{Cdouble}, b::Vector{Cdouble},
                   eps)
    if length(c) != length(a) || length(c) != length(b)
        error("input vectors should have same length")
    end

    ret = ccall((:toeplitz, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
                t, a, b, length(t), eps)
    if ret != 0
        error("failed to solve a symmetric Toeplitz set of linear equations")
    end

    a
end

function toeplitz(t::Vector{Cdouble}, b::Vector{Cdouble};
                  eps::Float64=1.0e-6)
    a = similar(t)
    toeplitz!(t, a, b, eps)
end
