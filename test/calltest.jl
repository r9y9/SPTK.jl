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

    # Mel-Cepstrum related analysis
    mcep(dummy_input, order=20, alpha=0.41)
    mgcep(dummy_input, order=20, alpha=0.41, gamma=-1/4)
    fftcep(dummy_input, 20)
    mfcc(dummy_input, order=12)

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
end

check_no_segfault()
