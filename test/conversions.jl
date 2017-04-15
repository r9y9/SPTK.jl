# Conversions

import SPTK: lpc2c!

function test_conversion_base(f::Function, f!::Function,
                              src_order, dst_order, args...;
                              kargs...)
    # expected dst length
    ka = Dict{Symbol,Any}(kargs)
    dst_length = haskey(ka, :dst_length) ? ka[:dst_length] : dst_order + 1

    srand(98765)
    src = rand(src_order + 1)
    dst_inplace = zeros(dst_length)
    @assert applicable(f, src, dst_order, args...)
    @assert applicable(f!, dst_inplace, src, args...)

    dst = f(src, dst_order, args...)
    f!(dst_inplace, src, args...)
    @test all(isfinite.(dst))
    @test dst ≈ dst_inplace
end

function test_transform_base(f::Function, f!::Function, order, args...)
    srand(98765)
    src = rand(order + 1)
    dst_inplace = zeros(order + 1)
    @assert applicable(f, src, args...)
    @assert applicable(f!, dst_inplace, src, args...)

    dst = f(src, args...)
    f!(dst_inplace, src, args...)
    @test all(isfinite.(dst))
    @test dst ≈ dst_inplace
end

function test_lpc2lsp(order)
    srand(98765)
    dummy_lpc = lpc(rand(512), order)
    lsp_inplace = zeros(order + 1)

    lsp = lpc2lsp(dummy_lpc)
    lpc2lsp!(lsp_inplace, dummy_lpc)
    @test all(isfinite.(lsp))
    @test lsp ≈ lsp_inplace
end

function test_ic2ir_invertibility(order, len)
    srand(98765)
    dummy_ceps = rand(order+1)
    ir = c2ir(dummy_ceps, len)
    c = ic2ir(ir, order)
    @test c ≈ dummy_ceps
end

function test_gc2gc(src_order, dst_order, src_γ, dst_γ)
    srand(98765)
    src_ceps = rand(src_order + 1)
    dst_ceps_inplace = zeros(dst_order + 1)

    dst_ceps = gc2gc(src_ceps, src_γ, dst_order, dst_γ)
    gc2gc!(dst_ceps_inplace, dst_γ, src_ceps, src_γ)
    @test all(isfinite.(dst_ceps))
    @test dst_ceps ≈ dst_ceps_inplace
end

function test_mgc2mgc(dst_order, dst_α, dst_γ, src_order, src_α, src_γ)
    srand(98765)
    src_ceps = rand(src_order + 1)
    dst_ceps_inplace = zeros(dst_order + 1)

    dst_ceps = mgc2mgc(src_ceps, src_α, src_γ, dst_order, dst_α, dst_γ)
    mgc2mgc!(dst_ceps_inplace, dst_α, dst_γ, src_ceps, src_α, src_γ)
    @test all(isfinite.(dst_ceps))
    @test dst_ceps ≈ dst_ceps_inplace
end

function test_mgc2sp(order, α, γ, fftlen)
    srand(98765)
    dummy_ceps = rand(order + 1)
    sp = mgc2sp(dummy_ceps, α, γ, fftlen)
    @test length(sp) == fftlen>>1 + 1
    @test all(isfinite.(sp))
end

function test_mgclsp2sp(order, α, γ, fftlen)
    srand(98765)
    dummy_lsp = rand(order + 1)
    sp = mgclsp2sp(dummy_lsp, α, γ, fftlen)
    @test length(sp) == fftlen>>1 + 1
    @test all(isfinite.(sp))
end

println("-- test_lpc2c")
for src_order in [15, 20, 25, 30]
    for dst_order in [15, 20, 25, 30]
        test_conversion_base(lpc2c, lpc2c!, src_order, dst_order)
    end
end

println("-- test_lpc2lsp")
for order in [15, 20, 25, 30, 40, 50]
    println(" where order = $order")
    test_transform_base(lpc2lsp, lpc2lsp!, order)
    test_lpc2lsp(order)
end

let
    srand(98765)
    dummy_lpc = lpc(rand(512), 21)
    @test_throws ErrorException lpc2lsp(dummy_lpc, otype=2)
    @test_throws ErrorException lpc2lsp(dummy_lpc, otype=3)
    lsp1 = lpc2lsp(dummy_lpc, otype=2, fs=16000)
    lsp2 = lpc2lsp(dummy_lpc, otype=3, fs=16)
    @test lsp1 ≈ lsp2
    lsp3 = lpc2lsp(dummy_lpc, otype=3, fs=16, loggain=true)
    @test first(lsp3) ≈ log(first(lsp2))
end

let
    srand(98765)
    dummy_lpc = lpc(rand(512), 21)
    dummy_lsp_inplace = rand(50)
    @test_throws DimensionMismatch lpc2lsp!(dummy_lsp_inplace, dummy_lpc)
end

println("-- test_lpc2par")
for order in [15, 20, 25, 30]
    println(" where order = $order")
    test_transform_base(lpc2par, lpc2par!, order)
end

let
    @test_throws DimensionMismatch lpc2par!(ones(10), ones(9))
end

println("-- test_par2lpc")
for order in [15, 20, 25, 30]
    println(" where order = $order")
    test_transform_base(par2lpc, par2lpc!, order)
end

let
    @test_throws DimensionMismatch lpc2par!(ones(10), ones(9))
end

println("-- test_lsp2sp")
for order in [15, 20, 25, 30]
    for fftlen in [256, 512, 1024]
        println(" where order = $order, fftlen = $fftlen")
        test_conversion_base(lsp2sp, lsp2sp!, order, fftlen;
                             dst_length=fftlen>>1 + 1)
    end
end

let
    try lsp2sp(ones(20), 512) catch @test false end
    try lsp2sp!(ones(257), ones(20)) catch @test false end
    # invalid fftlen
    @test_throws ArgumentError lsp2sp(ones(20), 513)
    @test_throws ArgumentError lsp2sp!(ones(256), ones(20))
end

println("-- test_mc2b")
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        println(" where order = $order, α = $α")
        test_transform_base(mc2b, mc2b!, order, α)
    end
end

let
    @test_throws DimensionMismatch mc2b!(ones(10), ones(9), 0.41)
end

println("-- test_b2mc")
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        println(" where order = $order, α = $α")
        test_transform_base(b2mc, b2mc!, order, α)
    end
end

let
    @test_throws DimensionMismatch b2mc!(ones(10), ones(9), 0.41)
end

println("-- test_b2c")
for src_order in [15, 20, 25, 30]
    for dst_order in [15, 20, 25, 30]
        for α in [0.35, 0.41, 0.5]
            println(" where src_order = $src_order, dst_order = $dst_order, α = $α")
            test_conversion_base(b2c, b2c!, src_order, dst_order, α)
        end
    end
end

println("-- test_c2acr")
for src_order in [15, 20, 25, 30]
    for dst_order in [15, 20, 25, 30]
        for fftlen in [256, 512, 1024]
            println(" where src_order = $src_order, dst_order = $dst_order, fftlen = $fftlen")
            test_conversion_base(c2acr, c2acr!, src_order, dst_order, fftlen)
        end
    end
end

let
    # invalid fftlen
    @test_throws ArgumentError c2acr(ones(20), 19, 513)
end

println("-- test_c2ir")
for order in [15, 20, 25, 30]
    for len in [256, 512, 1024]
        println(" where order = $order, len = $len")
        test_conversion_base(c2ir, c2ir!, order, len; dst_length=len)
    end
end

println("-- test_ic2ir")
for len in [256, 512, 1024]
    for order in [15, 20, 25, 30]
        println(" where len = $len, order = $order")
        test_conversion_base(ic2ir, ic2ir!, len, order)
    end
end

println("-- test_ic2ir invertibility")
for order in [15, 20, 25, 30]
    for len in [256, 512, 1024]
        println(" where order = $order, len = $len")
        test_ic2ir_invertibility(order, len)
    end
end

println("-- test_gc2gc")
for src_order in [15, 20, 25, 30]
    for dst_order in [15, 20, 25, 30]
        for src_γ in [-1.0, -0.5, 0.0]
            for dst_γ in [-1.0, -0.5, 0.0]
                println(" where dst_order = $dst_order, dst_γ = $dst_γ, src_order = $src_order, src_γ = $src_γ")
                test_gc2gc(src_order, dst_order, src_γ, dst_γ)
            end
        end
    end
end

let
    @test_throws ArgumentError gc2gc!(ones(20), -0.1, ones(20), 0.1)
    @test_throws ArgumentError gc2gc!(ones(20), 0.1, ones(20), -0.1)
end

let
    srand(98765)
    order = 20
    src = rand(order + 1)
    @test gc2gc(src, 0.0, order, -0.1) ≈ gc2gc(src, 0.0, -0.1)
end

println("-- test_gnorm")
for order in [15, 20, 25, 30]
    for γ in [-1.0, -0.5, 0.0]
        println(" where order = $order, γ = $γ")
        test_transform_base(gnorm, gnorm!, order, γ)
    end
end

println("-- test_ignorm")
for order in [15, 20, 25, 30]
    for γ in [-1.0, -0.5, 0.0]
        println(" where order = $order, γ = $γ")
        test_transform_base(ignorm, ignorm!, order, γ)
    end
end

let
    c1 = ones(10)
    c2 = ones(11)
    @test_throws DimensionMismatch gnorm!(c1, c2, 0.0)
    @test_throws DimensionMismatch ignorm!(c1, c2, 0.0)
    @test_throws ArgumentError gnorm!(c1, c2, 0.1)
    @test_throws ArgumentError ignorm!(c1, c2, 0.1)
end

println("-- test_freqt and fqtr")
for src_order in [15, 20, 25, 30]
    for dst_order in [15, 20, 25, 30]
        for α in [0.35, 0.41, 0.5]
            println(" where dst_order = $dst_order, src_order = $src_order, α = $α")
            test_conversion_base(freqt, freqt!, src_order, dst_order, α)
            test_conversion_base(frqtr, frqtr!, src_order, dst_order, α)
        end
    end
end

let
    srand(98765)
    order = 20
    src = rand(order + 1)
    @test freqt(src, order, 0.41) ≈ freqt(src, 0.41)
    @test frqtr(src, order, 0.41) ≈ frqtr(src, 0.41)
end

println("-- test_mgc2mgc")
for dst_order in [15, 20, 25, 30]
    for dst_α in [-0.35, 0.41, 0.5]
        for src_γ in [-1.0, -0.5, 0.0]
            for src_order in [15, 20, 25, 30]
                for src_α in [0.35, 0.41, 0.5]
                    for dst_γ in [-1.0, -0.5, 0.0]
                        println(" where dst_order = $dst_order, dst_α = $dst_α, dst_γ = $dst_γ, src_order = $src_order, src_α = $src_α, dst_γ = $dst_γ")
                        test_mgc2mgc(dst_order, dst_α, dst_γ, src_order, src_α, src_γ)
                    end
                end
            end
        end
    end
end

let
    c1 = ones(20)
    c2 = copy(c1)
    @test_throws ArgumentError mgc2mgc!(c2, 0.0, -0.1, c1, 0.0, 1.0)
    @test_throws ArgumentError mgc2mgc!(c2, 0.0, 1.0, c1, 0.0, -0.1)
end

let
    srand(98765)
    order = 20
    src = rand(order + 1)
    @test mgc2mgc(src, 0.0, 0.0, order, 0.41, -0.1) ≈ mgc2mgc(src, 0.0, 0.0, 0.41, -0.1)
end

println("-- test_mgc2sp")
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        for γ in [-1.0, -0.5, 0.0]
            for fftlen in [256, 512, 1024]
                println(" where order = $order, α = $α, γ = $γ, fftlen = $fftlen")
                test_mgc2sp(order, α, γ, fftlen)
            end
        end
    end
end

function mgc2sp_exceptions()
    try mgc2sp(zeros(20), 0.0, 0.0, 512) catch @test false end
    let sp = zeros(Complex{Cdouble}, 513)
        try mgc2sp!(sp, zeros(20), 0.0, 0.0) catch @test false end
    end
    # invalid fftlen
    @test_throws ArgumentError mgc2sp(zeros(20), 0.0, 0.0, 513)
    let sp = zeros(Complex{Cdouble}, 512)
        @test_throws ArgumentError mgc2sp!(sp, zeros(20), 0.0, 0.0)
    end

    # invalid γ
    @test_throws ArgumentError mgc2sp(zeros(20), 0.0, 1.0, 513)
end

mgc2sp_exceptions()

println("-- test_mgclsp2sp")

# TODO: Inf/-Inf will happens when γ = 0
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        for γ in [-1.0, -0.5] # , 0.0]
            for fftlen in [256, 512, 1024]
                println(" where order = $order, α = $α, γ = $γ, fftlen = $fftlen")
                test_mgclsp2sp(order, α, γ, fftlen)
            end
        end
    end
end

let
    try mgclsp2sp(zeros(20), 0.0, 0.0, 512) catch @test false end
    try mgclsp2sp!(zeros(513), zeros(20), 0.0, 0.0) catch @test false end
    # invalid γ
    @test_throws ArgumentError mgclsp2sp!(zeros(513), zeros(20), 0.0, 1.0)
    # invalid fftlen
    @test_throws ArgumentError mgclsp2sp(zeros(20), 0.0, 0.0, 513)
    @test_throws ArgumentError mgclsp2sp!(zeros(512), zeros(20), 0.0, 0.0)
end

println("-- test_c2ndps")
for order in [15, 20, 25, 30]
    for fftlen in [256, 512, 1024]
        println(" where order = $order, fftlen = $fftlen")
        test_conversion_base(c2ndps, c2ndps!, order, fftlen,
                             dst_length=fftlen>>1 + 1)
    end
end

let
    try c2ndps(zeros(20), 512) catch @test false end
    try c2ndps!(zeros(257), zeros(20)) catch @test false end
    # invalid fftlen
    @test_throws ArgumentError c2ndps(zeros(20), 513)
    @test_throws ArgumentError c2ndps!(zeros(256), zeros(20))
end

println("-- test_ndps2c")
for order in [15, 20, 25, 30]
    for fftlen in [256, 512, 1024]
        println(" where order = $order, fftlen = $fftlen")
        test_conversion_base(ndps2c, ndps2c!, fftlen>>1, order)
    end
end

let
    try ndps2c(zeros(257), 20) catch @test false end
    # invalid fftlen
    @test_throws ArgumentError ndps2c(zeros(256), 20)
end
