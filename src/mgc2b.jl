# mgc2b converts mel generalized cesptrum to MGLSADF filter coefficients.
function mgc2b(mgc::Vector{Float64}, alpha::Float64, gamma::Float64)
    b = mc2b(mgc, alpha)

    # when gamma = 0, mel-generalized cespstrum corresponds to mel cepstrum
    if gamma == 0.0
        return b
    end

    b = gnorm(b, gamma) # TODO replace with inplace version

    # scale by gamma
    b[1] = log(b[1])
    b[2:end] *= gamma

    return b
end
