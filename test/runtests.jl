using SPTK
using Base.Test

function test_mgcep()
    println("testing: Mel-generalized cesptrum analysis")
    srand(98765)
    dummy_input = rand(1024)
    dummy_input_mat = repmat(dummy_input, 1, 10)

    c = mcep(dummy_input, 20, 0.41)
    @test length(c) == 21
    cmat = mcep(dummy_input_mat, 20, 0.41)
    @test size(cmat) == (21, 10)

    c = gcep(dummy_input, 20, -1/4)
    @test length(c) == 21
    cmat = gcep(dummy_input_mat, 20, -1/4)
    @test size(cmat) == (21, 10)

    c = mgcep(dummy_input, 20, 0.41, -1/4)
    @test length(c) == 21
    cmat = mgcep(dummy_input_mat, 20, 0.41, -1/4)
    @test size(cmat) == (21, 10)

    c = uels(dummy_input, 20)
    @test length(c) == 21
    cmat = uels(dummy_input_mat, 20)
    @test size(cmat) == (21, 10)

    c = fftcep(dummy_input, 20)
    cmat = fftcep(dummy_input_mat, 20)
    @test length(c) == 21
    @test size(cmat) == (21, 10)

    # mgcep when gamma = 0, the result of mgcep is corresponds to mcep
    mc = mcep(dummy_input, 20, 0.41)
    mgc = mgcep(dummy_input, 20, 0.41, 0.0)
    @test_approx_eq_eps mc mgc 1.0e-4

    sp = mgclsp2sp(dummy_input, 0.41, -1/4, 512)
    @test length(sp) == 256+1

    # refs #5
    # change order 20 -> 40
    mcep(dummy_input, 40, 0.41)
    mgcep(dummy_input, 40, 0.41, 0.0)
end

function test_mfcc()
    println("testing: MFCC")
    srand(20121)
    dummy = rand(1024)
    dummy_mat = repmat(dummy, 1, 10)

    # By default, c0 is not contained
    cc = mfcc(dummy, 12)
    @test length(cc) == 12
    ccmat = mfcc(dummy_mat, 12)
    @test size(ccmat) == (12, 10)

    # with c0
    cc = mfcc(dummy, 12; czero=true)
    @test length(cc) == 13
    ccmat = mfcc(dummy_mat, 12; czero=true)
    @test size(ccmat) == (13, 10)

    # with c0 + power
    cc = mfcc(dummy, 12; czero=true, power=true)
    @test length(cc) == 14
    ccmat = mfcc(dummy_mat, 12; czero=true, power=true)
    @test size(ccmat) == (14, 10)
end

function test_mgcep_conversion()
    println("testing: Mel-genearlized cepstrum conversions")

    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)
    dummy_ceps_mat = repmat(dummy_ceps, 1, 10)

    b = mc2b(dummy_ceps, 0.41)
    @test length(b) == length(dummy_ceps)
    bmat = mc2b(dummy_ceps_mat, 0.41)
    @test size(bmat) == (length(dummy_ceps), 10)

    c = b2mc(dummy_ceps, 0.41)
    @test length(c) == length(dummy_ceps)
    cmat = b2mc(dummy_ceps_mat, 0.41)
    @test size(cmat) == (length(dummy_ceps), 10)

    c = c2ir(dummy_ceps, 512)
    @test length(c) == 512
    cmat = c2ir(dummy_ceps_mat, 512)
    @test size(cmat) == (512, 10)

    ndps = c2ndps(dummy_ceps, 512)
    @test length(ndps) == div(512, 2) + 1
    ndpsmat = c2ndps(dummy_ceps_mat, 512)
    @test size(ndpsmat) == (div(512, 2) + 1, 10)

    c = ndps2c(ndps, 20)
    @test length(c) == 21
    cmat = ndps2c(repmat(ndps, 1, 10), 20)
    @test size(cmat) == (21, 10)

    c = gc2gc(dummy_ceps, 0.0, 15, -1/4)
    @test length(c) == 16
    cmat = gc2gc(dummy_ceps_mat, 0.0, 15, -1/4)
    @test size(cmat) == (16, 10)

    c = gnorm(dummy_ceps, -1/4)
    @test length(c) == length(dummy_ceps)
    cmat = gnorm(dummy_ceps_mat, -1/4)
    @test size(cmat) == (length(dummy_ceps), 10)

    c = ignorm(dummy_ceps, -1/4)
    @test length(c) == length(dummy_ceps)
    cmat = ignorm(dummy_ceps_mat, -1/4)
    @test size(cmat) == (length(dummy_ceps), 10)

    c = freqt(dummy_ceps, 22, 0.41)
    @test length(c) == 23
    cmat = freqt(dummy_ceps_mat, 22, 0.41)
    @test size(cmat) == (23, 10)

    c = mgc2mgc(dummy_ceps, 0.41, 0.0, 22, 0.41, -1/4)
    @test length(c) == 23
    cmat = mgc2mgc(dummy_ceps_mat, 0.41, 0.0, 22, 0.41, -1/4)
    @test size(cmat) == (23, 10)
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

function test_windows()
    srand(98765)
    x = rand(1024)
    xmat = rand(1024, 10)

    y = blackman(x)
    @test length(y) == length(x)
    y = hamming(x)
    @test length(y) == length(x)
    y = hanning(x)
    @test length(y) == length(x)
    y = bartlett(x)
    @test length(y) == length(x)
    y = trapezoid(x)
    @test length(y) == length(x)
    y = rectangular(x)
    @test length(y) == length(x)

    y = blackman(xmat)
    @test size(y) == size(xmat)
    y = hamming(xmat)
    @test size(y) == size(xmat)
    y = hanning(xmat)
    @test size(y) == size(xmat)
    y = bartlett(xmat)
    @test size(y) == size(xmat)
    y = trapezoid(xmat)
    @test size(y) == size(xmat)
    y = rectangular(xmat)
    @test size(y) == size(xmat)
end

test_mgcep()
test_mfcc()
test_mgcep_conversion()
test_f0()
test_waveform_generation()
test_windows()
