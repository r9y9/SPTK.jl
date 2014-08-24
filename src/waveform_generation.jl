# MLSADF represents a Mel Log Spectrum Approximation (MLSA) digtal filter
type MLSADF
    order::Int              # order of mel-cepstrum
    pd::Int                 # order of pade approximation
    delay::Vector{Float64}  # filter delay
end

# see mlsadf.c in the original SPTK for this magic allocation
mlsadf_delay(order::Int, pd::Int) = zeros(3*(pd+1)+pd*(order+2))

MLSADF(order::Int; pd::Int=5) = MLSADF(order, pd, mlsadf_delay(order, pd))

# filter! modifies MLSADF delay.
function filter!(m::MLSADF, x::Real, b::Vector{Float64}, alpha::Float64=0.41)
    order = length(b) - 1
    order == m.order ||
        throw(DimensionMismatch("Order of mel-cepstrum may be wrong."))

    return mlsadf(float64(x), b, alpha, m.pd, m.delay)
end


# MLSASynthesizer represents a waveform synthesizer based on MLSA digital
# filter.
type MLSASynthesizer
    f::MLSADF
end

MLSASynthesizer(order::Int; pd::Int=5) = MLSASynthesizer(MLSADF(order, pd=pd))

# synthesis_one_frame! generates speech waveform for one frame speech signal
# given a excitation signal and successive two mel-cepstrum.
function synthesis_one_frame!(s::MLSASynthesizer, excite::Vector{Float64},
                              previous_mcep::Vector{Float64},
                              current_mcep::Vector{Float64},
                              alpha::Float64=0.41)

    previous_coef = mc2b(previous_mcep, alpha)
    current_coef = mc2b(current_mcep, alpha)

    slope = (current_coef - previous_coef) / float(length(excite))

    part_of_speech = zeros(length(excite))
    interpolated_coef = copy(previous_coef)

    for i=1:endof(excite)
        scaled_excitation = excite[i] * exp(interpolated_coef[1])
        part_of_speech[i] = filter!(s.f, scaled_excitation,
                                    interpolated_coef, alpha)
        interpolated_coef += slope
    end

    return part_of_speech
end

# synthesis! generates a speech waveform given a excitation signal and
# a sequence of mel-cepstrum.
function synthesis!(s::MLSASynthesizer, excite::Vector{Float64},
                    mcep_sequence::Matrix{Float64},
                    alpha::Float64=0.41, hopsize::Int=80)
    synthesized = zeros(length(excite))

    previous_mcep = mcep_sequence[1]
    for i=1:size(mcep_sequence, 1)
        if i > 1
            previous_mcep = mcep_sequence[i-1]
        end
        current_mcep = mcep_sequence[i]

        s, e = (i-1)*hopsize+1, i*hopsize
        if e >= length(excite)
            break
        end

        part_of_speech = synthesis_one_frame!(s, excite[s:e],
                                              previous_mcep,
                                              current_mcep, alpha)
        synthesized[s:e] = part_of_speech
    end

    return synthesized
end
