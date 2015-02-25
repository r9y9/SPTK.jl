using SPTK
using Base.Test

let
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

let
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

let
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

    c = SPTK.frqtr(dummy_ceps, 22, 0.41)
    @test length(c) == 23
    cmat = SPTK.frqtr(dummy_ceps_mat, 22, 0.41)
    @test size(cmat) == (23, 10)

    c = mgc2mgc(dummy_ceps, 0.41, 0.0, 22, 0.41, -1/4)
    @test length(c) == 23
    cmat = mgc2mgc(dummy_ceps_mat, 0.41, 0.0, 22, 0.41, -1/4)
    @test size(cmat) == (23, 10)

    sp = mgc2sp(dummy_ceps, 0.41, 0.0, 1024)
    @test length(sp) == 1024>>1+1
    spmat = mgc2sp(dummy_ceps_mat, 0.41, 0.0, 1024)
    @test size(spmat) == (1024>>1+1, 10)
end

let
    println("testing: LPC")

    srand(98765)
    dummy_input = rand(1024)

    l = lpc(dummy_input, 20)
    @test length(l) == 21
    @test !any(isnan(l))

    println("testing: LPC conversions")

    c = lpc2c(l)
    @test length(c) == length(l)
    @test !any(isnan(c))

    lsp = lpc2lsp(l, 20)
    @test length(lsp) == 21
    @test !any(isnan(lsp))

    par = lpc2par(l)
    @test length(par) == 21
    @test !any(isnan(par))

    sp = lsp2sp(lsp, 1024)
    @test length(sp) == 1024>>1+1
    @test !any(isnan(sp))

    try lspcheck(l); catch @assert false; end
end

let
    println("testing: F0 estimation")

    srand(98765)
    dummy_input = rand(1024)

    f0 = swipe(dummy_input, 16000, 80)
    @test length(f0) == div(length(dummy_input), 80) + 1
    @test !any(isnan(f0))
    @test all(f0 .>= 0)
end

let
    println("testing: Waveform generation filters")

    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)
    dummy_lpc = rand(20)

    pd = 5
    order = length(dummy_ceps)-1

    # poledf
    d = poledf_delay(order)
    for x in dummy_input
        y = poledf(x, dummy_lpc, d)
        @test !isnan(y)
    end

    # lmadf
    d = lmadf_delay(order, pd)
    for x in dummy_input
        y = lmadf(x, dummy_ceps, pd, d)
        @test !isnan(y)
    end

    # mlsadf
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

let
    println("testing: Window functions")

    srand(98765)
    x = rand(1024)

    y = blackman(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    y = hamming(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    y = hanning(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    y = bartlett(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    y = trapezoid(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    y = rectangular(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))
end
