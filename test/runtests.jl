using SPTK
using Base.Test

# Library routines

println("test library routines")

function test_agexp()
    println("-- test_agexp")
    agexp(1, 1, 1)
end

function test_cholesky()
    println("-- test_choleskey")
    srand(98765)
    a = rand(10)
    b = copy(a)
    c = copy(a)
    @test_throws Exception cholesky!(c, a, b, 1.0e-6)
    @test_throws Exception cholesky(c, b)
end

function test_gexp()
    println("-- test_gexp")
    gexp(1, 1)
end

function test_glog()
    println("-- test_glog")
    glog(1, 1)
end

function test_mseq()
    println("-- test_mseq")
    mseq()
    mseq()
    mseq()
end

function test_theq()
    println("-- test_theq")
    srand(98765)
    t = rand(10)
    h = copy(t)
    a = copy(t)
    b = copy(t)
    @test_throws Exception theq!(t, h, a, b, 1.0e-6)
    @test_throws Exception theq(t, h, b)
end

function test_toeplitz()
    println("-- test_toeplitz")
    srand(98765)
    a = rand(10)
    b = copy(a)
    c = copy(a)
    @test_throws Exception toeplitz!(c, a, b, 1.0e-6)
    @test_throws Exception toeplitz(c, b)
end

test_agexp()
test_cholesky()
test_gexp()
test_glog()
test_mseq()
test_theq()
test_toeplitz()

# SPTK APIs

function test_adaptive_mcep()
    println("test adaptive mel-cepstrum analysis")
    srand(98765)
    dummy_input = rand(64)
    order = 20
    c = zeros(order+1)

    println("-- test_acep!")
    for x in dummy_input
        acep!(c, x)
    end
    @assert !any(isnan(c))

    println("-- test_agcep!")
    for stage in 1:10
        fill!(c, zero(eltype(c)))
        for x in dummy_input
            agcep!(c, x, stage)
        end
        @assert !any(isnan(c))
    end

    println("-- test_amcep!")
    for α in[0.35, 0.41, 0.5]0
        fill!(c, zero(eltype(c)))
        delay = zeros(length(c))
        for x in dummy_input
            amcep!(c, x, α)
            phidf!(x, length(c)-1, α, delay)
        end
        @assert !any(isnan(c))
    end
end

test_adaptive_mcep()

function test_mgcep()
    println("test mel-generalized cepstrum analysis")
    srand(98765)
    dummy_input = rand(1024)
    dummy_input_mat = repmat(dummy_input, 1, 10)

    println("-- test_mcep")
    c = mcep(dummy_input, 20, 0.41)
    @test length(c) == 21
    cmat = mcep(dummy_input_mat, 20, 0.41)
    @test size(cmat) == (21, 10)

    println("-- test_gcep")
    c = gcep(dummy_input, 20, -1/4)
    @test length(c) == 21
    cmat = gcep(dummy_input_mat, 20, -1/4)
    @test size(cmat) == (21, 10)

    println("-- test_mgcep")
    c = mgcep(dummy_input, 20, 0.41, -1/4)
    @test length(c) == 21
    cmat = mgcep(dummy_input_mat, 20, 0.41, -1/4)
    @test size(cmat) == (21, 10)

    # TODO fix in windows
    @unix_only begin
        println("-- test_uels")
        c = uels(dummy_input, 20)
        @test length(c) == 21
        cmat = uels(dummy_input_mat, 20)
        @test size(cmat) == (21, 10)
    end

    println("-- test_fftcep")
    c = fftcep(dummy_input, 20)
    cmat = fftcep(dummy_input_mat, 20)
    @test length(c) == 21
    @test size(cmat) == (21, 10)

    println("-- test mcep and mgcep consistency")
    # mgcep when gamma = 0, the result of mgcep is corresponds to mcep
    mc = mcep(dummy_input, 20, 0.41)
    mgc = mgcep(dummy_input, 20, 0.41, 0.0)
    @test_approx_eq_eps mc mgc 1.0e-4

    println("-- test_mgclsp2sp")
    sp = mgclsp2sp(dummy_input, 0.41, -1/4, 512)
    @test length(sp) == 256+1

    # refs #5
    # change order 20 -> 40
    println("-- test changing order (ref: issue #5)")
    mcep(dummy_input, 40, 0.41)
    mgcep(dummy_input, 40, 0.41, 0.0)
end

test_mgcep()

function test_mfcc()
    println("test_mfcc")
    srand(20121)
    dummy = rand(1024)
    dummy_mat = repmat(dummy, 1, 10)

    # By default, c0 is not contained
    println("-- test no c0 is expected")
    cc = mfcc(dummy, 12)
    @test length(cc) == 12
    ccmat = mfcc(dummy_mat, 12)
    @test size(ccmat) == (12, 10)

    # with c0
    println("-- test c0 is expected")
    cc = mfcc(dummy, 12; czero=true)
    @test length(cc) == 13
    ccmat = mfcc(dummy_mat, 12; czero=true)
    @test size(ccmat) == (13, 10)

    # with c0 + power
    println("-- test c0 and power are expected")
    cc = mfcc(dummy, 12; czero=true, power=true)
    @test length(cc) == 14
    ccmat = mfcc(dummy_mat, 12; czero=true, power=true)
    @test size(ccmat) == (14, 10)

    # with power
    println("-- test power is expected")
    cc = mfcc(dummy, 12; czero=false, power=true)
    @test length(cc) == 13
    ccmat = mfcc(dummy_mat, 12; czero=false, power=true)
    @test size(ccmat) == (13, 10)
end

test_mfcc()

function test_mgcep_conversions()
    println("test mgcep conversions")

    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)
    dummy_ceps_mat = repmat(dummy_ceps, 1, 10)

    println("-- test_mc2b")
    b = mc2b(dummy_ceps, 0.41)
    @test length(b) == length(dummy_ceps)
    bmat = mc2b(dummy_ceps_mat, 0.41)
    @test size(bmat) == (length(dummy_ceps), 10)

    println("-- test_b2mc")
    c = b2mc(dummy_ceps, 0.41)
    @test length(c) == length(dummy_ceps)
    cmat = b2mc(dummy_ceps_mat, 0.41)
    @test size(cmat) == (length(dummy_ceps), 10)

    println("-- test_c2ir")
    c = c2ir(dummy_ceps, 512)
    @test length(c) == 512
    cmat = c2ir(dummy_ceps_mat, 512)
    @test size(cmat) == (512, 10)

    println("-- test_c2ndps")
    ndps = c2ndps(dummy_ceps, 512)
    @test length(ndps) == div(512, 2) + 1
    ndpsmat = c2ndps(dummy_ceps_mat, 512)
    @test size(ndpsmat) == (div(512, 2) + 1, 10)

    println("-- test_ndps2c")
    c = ndps2c(ndps, 20)
    @test length(c) == 21
    cmat = ndps2c(repmat(ndps, 1, 10), 20)
    @test size(cmat) == (21, 10)

    println("-- test_gc2gc")
    c = gc2gc(dummy_ceps, 0.0, 15, -1/4)
    @test length(c) == 16
    cmat = gc2gc(dummy_ceps_mat, 0.0, 15, -1/4)
    @test size(cmat) == (16, 10)

    println("-- test_gnorm")
    c = gnorm(dummy_ceps, -1/4)
    @test length(c) == length(dummy_ceps)
    cmat = gnorm(dummy_ceps_mat, -1/4)
    @test size(cmat) == (length(dummy_ceps), 10)

    println("-- test_ignorm")
    c = ignorm(dummy_ceps, -1/4)
    @test length(c) == length(dummy_ceps)
    cmat = ignorm(dummy_ceps_mat, -1/4)
    @test size(cmat) == (length(dummy_ceps), 10)

    println("-- test_freqt")
    c = freqt(dummy_ceps, 22, 0.41)
    @test length(c) == 23
    cmat = freqt(dummy_ceps_mat, 22, 0.41)
    @test size(cmat) == (23, 10)

    println("-- test_frqtr")
    c = SPTK.frqtr(dummy_ceps, 22, 0.41)
    @test length(c) == 23
    cmat = SPTK.frqtr(dummy_ceps_mat, 22, 0.41)
    @test size(cmat) == (23, 10)

    println("-- test_mgc2mgc")
    c = mgc2mgc(dummy_ceps, 0.41, 0.0, 22, 0.41, -1/4)
    @test length(c) == 23
    cmat = mgc2mgc(dummy_ceps_mat, 0.41, 0.0, 22, 0.41, -1/4)
    @test size(cmat) == (23, 10)

    println("-- test_mgc2sp")
    sp = mgc2sp(dummy_ceps, 0.41, 0.0, 1024)
    @test length(sp) == 1024>>1+1
    spmat = mgc2sp(dummy_ceps_mat, 0.41, 0.0, 1024)
    @test size(spmat) == (1024>>1+1, 10)
end

test_mgcep_conversions()

function test_lpc()
    println("test lpc")
    srand(98765)
    dummy_input = rand(1024)

    println("-- test_lpc")
    l = lpc(dummy_input, 20)
    @test length(l) == 21
    @test !any(isnan(l))

    println("-- test_lpc2c")
    c = lpc2c(l)
    @test length(c) == length(l)
    @test !any(isnan(c))

    println("-- test_lpc2lsp")
    lsp = lpc2lsp(l, 20)
    @test length(lsp) == 21
    @test !any(isnan(lsp))

    println("-- test_lpc2par")
    par = lpc2par(l)
    @test length(par) == 21
    @test !any(isnan(par))

    println("-- test_lsp2sp")
    sp = lsp2sp(lsp, 1024)
    @test length(sp) == 1024>>1+1
    @test !any(isnan(sp))

    println("-- test_lspcheck")
    try lspcheck(l); catch @assert false; end
end

test_lpc()

function test_f0()
    println("test f0 estimatino")
    srand(98765)
    dummy_input = rand(1024)

    println("-- test_swipe")
    f0 = swipe(dummy_input, 16000, 80)
    @test length(f0) == div(length(dummy_input), 80) + 1
    @test !any(isnan(f0))
    @test all(f0 .>= 0)
end

test_f0()

function test_waveform_generation_filters()
    println("test waveform generation filters")

    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs(fft(dummy_input))
    dummy_logsp = log(dummy_sp)
    dummy_ceps = rand(20)
    dummy_lpc = rand(20)

    pd = 5
    order = length(dummy_ceps)-1

    println("-- test_poledf")
    # poledf
    d = poledf_delay(order)
    for x in dummy_input
        y = poledf(x, dummy_lpc, d)
        @test !isnan(y)
    end

    println("-- test_lmadf")
    # lmadf
    d = lmadf_delay(order, pd)
    for x in dummy_input
        y = lmadf(x, dummy_ceps, pd, d)
        @test !isnan(y)
    end

    println("-- test_mlsadf")
    # mlsadf
    d = mlsadf_delay(order, pd)
    for x in dummy_input
        y = mlsadf(x, dummy_ceps, 0.41, pd, d)
        @test !isnan(y)
    end

    println("-- test_mglsadf")
    # mlgasdf
    stage = 12
    d = mglsadf_delay(order, stage)
    for x in dummy_input
        y = mglsadf(x, dummy_ceps, 0.41, stage, d)
        @test !isnan(y)
    end
end

test_waveform_generation_filters()

function test_window_functions()
    println("test windows functions")
    srand(98765)
    x = rand(1024)

    println("-- test_blackman")
    y = blackman(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    println("-- test_hamming")
    y = hamming(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    println("-- test_hanning")
    y = hanning(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    println("-- test_bartlett")
    y = bartlett(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    println("-- test_trapezoid")
    y = trapezoid(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))

    println("-- test_rectangular")
    y = rectangular(length(x))
    @test length(y) == length(x)
    @test !any(isnan(y))
end

test_window_functions()
