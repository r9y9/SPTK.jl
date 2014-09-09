using SPTK
using Base.Test

# Just check no segumentation fault happen
function check_no_segfault()
    # Create dummy variables
    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)

    # (Mel-) Cepstrum related analysis
    mcep(dummy_input, 20, 0.41)
    mgcep(dummy_input, 20, 0.41, -1/4)
    uels(dummy_input, 20)
    fftcep(dummy_input, 20)
    mfcc(dummy_input, 12)

    # Conversions
    mc2b(dummy_ceps, 0.41)
    b2mc(dummy_ceps, 0.41)
    c2ir(dummy_ceps, 512)
    gc2gc(dummy_ceps, 0.0, 15, -1/4)
    gnorm(dummy_ceps, -1/4)
    ignorm(dummy_ceps, -1/4)
    freqt(dummy_ceps, 22, 0.41)
    mgc2mgc(dummy_ceps, 0.41, 0.0, 22, 0.41, -1/4)

    # F0 estimation
    swipe(dummy_input, 16000)

    # Waveform generation filters
    # mlsadf
    pd::Int = 5
    order::Int = length(dummy_ceps)-1
    # see mlsadf.c in original SPTK for this magic allocation
    delay_mlsadf = zeros(3*(pd+1) + pd*(order+2))
    mlsadf(dummy_input[1], dummy_ceps, 0.41, pd, delay_mlsadf)

    # mlgasdf
    stage = 12
    # see mglsadf.c in original SPTK for this magic allocation
    delay_mglsadf = zeros((order+1)*stage)
    mglsadf(dummy_input[1], dummy_ceps, 0.41, stage, delay_mglsadf)
end

check_no_segfault()
