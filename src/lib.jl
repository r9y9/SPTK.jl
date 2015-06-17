# Library routines

function agexp(r, x, y)
    ccall((:agexp, libSPTK), Cdouble, (Cdouble, Cdouble, Cdouble), r, x, y)
end

# solve Ca = b
function cholesky!(c::Vector{Cdouble}, a::Vector{Cdouble}, b::Vector{Cdouble}, eps)
    if length(c) != length(a) || length(c) != length(b)
        error("input vectors should have same length")
    end

    ret = ccall((:cholesky, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cdouble),
                c, a, b, length(c), eps)
    if ret != 0
        error("failed to compute choleskey decomposition")
    end

    a
end

function cholesky!(c::Vector{Cdouble}, b::Vector{Cdouble};
                   eps::Float64=1.0e-6)
    a = similar(c)
    cholesky!(c, a, b, eps)
end

gexp(r, x) = ccall((:gexp, libSPTK), Cdouble, (Cdouble, Cdouble), r, x)
glog(r, x) = ccall((:glog, libSPTK), Cdouble, (Cdouble, Cdouble), r, x)

mseq() = ccall((:mseq, libSPTK), Cint, ())

# solve (T + H)a = b
function theq!(t::Vector{Cdouble}, h::Vector{Cdouble}, a::Vector{Cdouble},
               b::Vector{Cdouble}, eps)
    if length(t) != length(h) || length(t) != length(a) || length(t) != length(b)
        error("input vectors should have same length")
    end

    ret = ccall((:theq, libSPTK), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint,
                 Cdouble), t, h, a, b, length(t), eps)
    if ret != 0
        error("failed to solve Toeplitz plus Hankel matrix system (T + H)a = b")
    end

    a
end

function theq(t::Vector{Cdouble}, h::Vector{Cdouble}, b::Vector{Cdouble},
              eps::Float64=1.0e-6)
    a = similar(t)
    theq!(t, h, a, b, eps)
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

function toeplitz!(t::Vector{Cdouble}, b::Vector{Cdouble};
                   eps::Float64=1.0e-6)
    a = similar(t)
    toeplitz!(t, a, b, eps)
end
