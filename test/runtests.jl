using SPTK
using Base.Test

function check_no_segfault()
    # Create dummy variables
    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)

    # (Mel-) Cepstrum related analysis
    mcep(dummy_input, 20, 0.41)
    gcep(dummy_input, 20, -1/4)
    mgcep(dummy_input, 20, 0.41, -1/4)
    uels(dummy_input, 20)
    fftcep(dummy_input, 20)
    mfcc(dummy_input, 12)

    # Conversions
    mc2b(dummy_ceps, 0.41)
    b2mc(dummy_ceps, 0.41)
    c2ir(dummy_ceps, 512)
    c = c2ndps(dummy_ceps, 512)
    @test length(c) == 256 + 1
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
    d = mlsadf_delay(order, pd)
    mlsadf(dummy_input[1], dummy_ceps, 0.41, pd, d)

    # mlgasdf
    stage = 12
    d = mglsadf_delay(order, stage)
    mglsadf(dummy_input[1], dummy_ceps, 0.41, stage, d)
end

function testmfcc()
  srand(20121)
  dummy = rand(1024)

  # By default, c0 is not contained
  cc = mfcc(dummy, 12)
  @test length(cc) == 12

  # with c0
  cc = mfcc(dummy, 12; czero=true)
  @test length(cc) == 13

  # with c0 + power
  cc = mfcc(dummy, 12; czero=true, power=true)
  @test length(cc) == 14
end

# MGCep analysis when gamma = 0, MGCep analysis is corresponds to MCep analysis
function mgcep_as_special_case_of_mcep()
  srand(20121)
  dummy = rand(1024)

  mc = mcep(dummy, 20, 0.41)
  mgc = mgcep(dummy, 20, 0.41, 0.0)

  # Order + 0-th
  @test length(mc) == 21
  @test length(mgc) == 21

  @test_approx_eq_eps mc mgc 1.0e-4
end

check_no_segfault()
testmfcc()
mgcep_as_special_case_of_mcep()
