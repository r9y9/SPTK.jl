# Window functions

immutable WindowType
    w::Int
end

function _window!(wtype::WindowType, x::Vector{Cdouble}; normalize::Int=0)
    # normalize: Int
    #    0 : don't normalize
    #    1 : normalize by power
    #    2 : normalize by magnitude
    if normalize ∉ 0:2
        throw(ArgumentError("invalid normalize flag $normalize, must be ∈ 0:2"))
    end
    g = ccall((:window, libSPTK), Cdouble,
              (WindowType, Ptr{Cdouble}, Cint, Cint), wtype, x, length(x),
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
            _window!(WindowType($wtype), y; kargs...)
        end
    end
end
