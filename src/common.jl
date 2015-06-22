function assert_gamma(γ)
    if !(-1 <= γ <= 0.0)
        throw(ArgumentError("unsupported γ: must be -1 <= γ <= 0)"))
    end
end

function assert_pade(pade)
    if pade ∉ 4:5
        throw(ArgumentError("4 or 5 pade approximation is supported"))
    end
end

function assert_fftlen(fftlen)
    if !ispow2(fftlen)
        throw(ArgumentError("fftlen must be power of 2"))
    end
end
