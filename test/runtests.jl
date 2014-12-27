using SPTK
using Base.Test

function test_mgcep()
    println("testing: Mel-generalized cesptrum analysis")
    srand(98765)
    dummy_input = rand(1024)

    c = mcep(dummy_input, 20, 0.41)
    @test length(c) == 21
    c = gcep(dummy_input, 20, -1/4)
    @test length(c) == 21
    c = mgcep(dummy_input, 20, 0.41, -1/4)
    @test length(c) == 21
    c = uels(dummy_input, 20)
    @test length(c) == 21
    c = fftcep(dummy_input, 20)
    @test length(c) == 21

    # mgcep when gamma = 0, the result of mgcep is corresponds to mcep
    mc = mcep(dummy_input, 20, 0.41)
    mgc = mgcep(dummy_input, 20, 0.41, 0.0)
    @test_approx_eq_eps mc mgc 1.0e-4
end

# TODO
function fail_tests()
    srand(98765)
    dummy_input = rand(1024)

    # segfaults
    mc = mcep(dummy_input, 40, 0.41)
    mgc = mgcep(dummy_input, 20, 0.41, 0.0)
end

function test_mfcc()
    println("testing: MFCC")
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

function test_mgcep_conversion()
    println("testing: Mel-genearlized cepstrum conversions")

    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)

    b = mc2b(dummy_ceps, 0.41)
    @test length(b) == length(dummy_ceps)
    c = b2mc(dummy_ceps, 0.41)
    @test length(c) == length(dummy_ceps)
    c = c2ir(dummy_ceps, 512)
    @test length(c) == 512
    ndps = c2ndps(dummy_ceps, 512)
    @test length(ndps) == div(512, 2) + 1
    c = ndps2c(ndps, 20)
    @test length(c) == 21
    c = gc2gc(dummy_ceps, 0.0, 15, -1/4)
    @test length(c) == 16
    c = gnorm(dummy_ceps, -1/4)
    @test length(c) == length(dummy_ceps)
    c = ignorm(dummy_ceps, -1/4)
    @test length(c) == length(dummy_ceps)
    c = freqt(dummy_ceps, 22, 0.41)
    @test length(c) == 23
    c = mgc2mgc(dummy_ceps, 0.41, 0.0, 22, 0.41, -1/4)
    @test length(c) == 23
end

function test_f0()
    println("testing: F0 estimation")

    srand(98765)
    dummy_input = rand(1024)

    f0 = swipe(dummy_input, 16000, 80)
    @test length(f0) == div(length(dummy_input), 80) + 1
    @test !any(isnan(f0))
end

function test_waveform_generation()
    println("testing: Waveform generation filters")

    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)

    # mlsadf
    pd::Int = 5
    order::Int = length(dummy_ceps)-1
    d = mlsadf_delay(order, pd)
    for x in dummy_input
        y = mlsadf(x, dummy_ceps, 0.41, pd, d)
        @test !isnan(y)
    end

    # mlgasdf
    stage = 12
    d = mglsadf_delay(order, stage)
    for x in dummy_input
        y = mglsadf(x, dummy_ceps, 0.41, stage, d)
        @test !isnan(y)
    end
end

test_mgcep()
test_mfcc()
test_mgcep_conversion()
test_f0()
test_waveform_generation()
