# Window functions

struct Cwindow
    w::Cint

    function Cwindow(w)
        if w ∉ 0:5
            throw(ArgumentError("invalid window index $w: must be ∈ 0:5"))
        end
        new(w)
    end
end

function _window!(wtype::Cwindow, x::StridedVector{Cdouble}; normalize::Int=1)
    # normalize: Int
    #    0 : don't normalize
    #    1 : normalize by power
    #    2 : normalize by magnitude
    if normalize ∉ 0:2
        throw(ArgumentError("invalid normalize flag $normalize, must be ∈ 0:2"))
    end
    g = ccall((:window, libSPTK), Cdouble,
              (Cwindow, Ptr{Cdouble}, Cint, Cint), wtype, x, length(x),
              normalize)
    x
end

for (f, wtype) in [(:blackman, 0),
                   (:hamming, 1),
                   (:hanning, 2),
                   (:bartlett, 3),
                   (:trapezoid, 4),
                   (:rectangular, 5)]
    @eval begin
        function $f(n::Integer; kargs...)
            y = ones(Float64, n)
            _window!(Cwindow($wtype), y; kargs...)
        end
    end
end
