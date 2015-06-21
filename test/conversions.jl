# Conversions

function test_mc2b(order, α)
    srand(98765)
    dummy_ceps = rand(order+1)
    b = mc2b(dummy_ceps, α)
    @test length(b) == length(dummy_ceps)
    @test !any(isnan(b))
end

function test_b2mc(order, α)
    srand(98765)
    dummy_b = rand(order+1)
    c = b2mc(dummy_b, α)
    @test length(c) == length(dummy_b)
    @test !any(isnan(c))
end

function test_c2acr(order, desired_order, fftlen)
    srand(98765)
    dummy_ceps = rand(order+1)
    r = c2acr(dummy_ceps, desired_order, fftlen)
    @test length(r) == desired_order + 1
    @test !any(isnan(r))
end

function test_c2ir(order, len)
    srand(98765)
    dummy_ceps = rand(order+1)
    ir = c2ir(dummy_ceps, len)
    @test length(ir) == len
    @test !any(isnan(ir))
end

function test_ic2ir(len, order)
    srand(98765)
    dummy_ir = rand(len)
    c = ic2ir(dummy_ir, order)
    @test length(c) == order + 1
    @test !any(isnan(c))
end

function test_ic2ir_invertibility(order, len)
    srand(98765)
    dummy_ceps = rand(order+1)
    ir = c2ir(dummy_ceps, len)
    c = ic2ir(ir, order)
    @test_approx_eq c dummy_ceps
end

function test_gc2gc(src_order, dst_order, src_γ, dst_γ)
    srand(98765)
    src_ceps = rand(src_order + 1)
    dst_ceps_inplace = zeros(dst_order + 1)

    dst_ceps = gc2gc(src_ceps, src_γ, dst_order, dst_γ)
    gc2gc!(dst_ceps_inplace, dst_γ, src_ceps, src_γ)
    @test !any(isnan(dst_ceps))
    @test_approx_eq dst_ceps dst_ceps_inplace
end

function test_gnorm(order, γ)
    srand(98765)
    dummy_ceps = rand(order + 1)
    c = gnorm(dummy_ceps, γ)
    @test length(c) == length(dummy_ceps)
    @test !any(isnan(c))
    gnorm!(dummy_ceps, γ)
    @test_approx_eq dummy_ceps c
end

function test_ignorm(order, γ)
    srand(98765)
    dummy_ceps = rand(order + 1)
    c = ignorm(dummy_ceps, γ)
    @test length(c) == length(dummy_ceps)
    @test !any(isnan(c))
    ignorm!(dummy_ceps, γ)
    @test_approx_eq dummy_ceps c
end

function test_freqt(src_order, dst_order, α, f::Symbol=:freqt)
    f = eval(f)
    inplacef = eval(symbol(string(f, :!)))

    srand(98765)
    src_ceps = rand(src_order + 1)
    dst_ceps_inplace = zeros(dst_order + 1)

    dst_ceps = f(src_ceps, dst_order, α)
    inplacef(dst_ceps_inplace, src_ceps, α)

    @test !any(isnan(dst_ceps))
    @test_approx_eq dst_ceps dst_ceps_inplace
end

function test_mgc2mgc(dst_order, dst_α, dst_γ, src_order, src_α, src_γ)
    srand(98765)
    src_ceps = rand(src_order + 1)
    dst_ceps = zeros(dst_order + 1)
    mgc2mgc!(dst_ceps, dst_α, dst_γ, src_ceps, src_α, src_γ)
    @test !any(isnan(dst_ceps))
end

function test_mgc2sp(order, α, γ, fftlen)
    srand(98765)
    dummy_ceps = rand(order + 1)
    sp = mgc2sp(dummy_ceps, α, γ, fftlen)
    @test length(sp) == fftlen>>1 + 1
    @test !any(isnan(sp))
end

function test_c2ndps(order, fftlen)
    srand(98765)
    dummy_ceps = rand(order + 1)
    ndps = c2ndps(dummy_ceps, fftlen)
    @test length(ndps) == fftlen>>1 + 1
    @test !any(isnan(ndps))
end

function test_ndps2c(order, fftlen)
    srand(98765)
    dummy_ndps = rand(fftlen>>1 + 1)
    c = ndps2c(dummy_ndps, order)
    @test length(c) == order + 1
    @test !any(isnan(c))
end

println("-- test_mc2b")
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        println(" where order = $order, α = $α")
        test_mc2b(order, α)
    end
end

let
    @test_throws DimensionMismatch mc2b!(ones(10), ones(9), 0.41)
end

println("-- test_b2mc")
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        println(" where order = $order, α = $α")
        test_b2mc(order, α)
    end
end

let
    @test_throws DimensionMismatch b2mc!(ones(10), ones(9), 0.41)
end

println("-- test_c2acr")
for order in [15, 20, 25, 30]
    for desired_order in [15, 20, 25, 30]
        for fftlen in [256, 512, 1024]
            println(" where order = $order, desired_order = $desired_order, fftlen = $fftlen")
            test_c2acr(order, desired_order, fftlen)
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
        test_c2ir(order, len)
    end
end

println("-- test_ic2ir")
for len in [256, 512, 1024]
    for order in [15, 20, 25, 30]
        println(" where len = $len, order = $order")
        test_ic2ir(len, order)
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

println("-- test_gnorm")
for order in [15, 20, 25, 30]
    for γ in [-1.0, -0.5, 0.0]
        println(" where order = $order, γ = $γ")
        test_gnorm(order, γ)
    end
end

println("-- test_ignorm")
for order in [15, 20, 25, 30]
    for γ in [-1.0, -0.5, 0.0]
        println(" where order = $order, γ = $γ")
        test_ignorm(order, γ)
    end
end

println("-- test_freqt and fqtr")
for src_order in [15, 20, 25, 30]
    for dst_order in [15, 20, 25, 30]
        for α in [0.35, 0.41, 0.5]
            println(" where dst_order = $dst_order, src_order = $src_order, α = $α")
            test_freqt(src_order, dst_order, α, :freqt)
            test_freqt(src_order, dst_order, α, :frqtr)
        end
    end
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

let
    try mgc2sp(zeros(20), 0.0, 0.0, 512) catch @test false end
    # invalid fftlen
    @test_throws ArgumentError mgc2sp(zeros(20), 0.0, 0.0, 513)
end

println("-- test_c2ndps")
for order in [15, 20, 25, 30]
    for fftlen in [256, 512, 1024]
        println(" where order = $order, fftlen = $fftlen")
        test_c2ndps(order, fftlen)
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
        test_ndps2c(order, fftlen)
    end
end

let
    try ndps2c(zeros(257), 20) catch @test false end
    # invalid fftlen
    @test_throws ArgumentError ndps2c(zeros(256), 20)
end
