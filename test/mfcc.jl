import SPTK: mfcc

function test_mfcc(len, order)
    srand(98765)
    dummy = rand(len)

    # By default, c0 is not contained
    cc = mfcc(dummy, order, czero=true, power=true)
    # @test length(cc) == order
    if any(isnan.(cc))
        @show cc
    end
    @test all(isfinite.(cc))
end

function test_mfcc_options()
    srand(98765)
    dummy = rand(512)

    try mfcc(dummy, 19, numfilterbunks=20) catch @test false end
    @test_throws ArgumentError mfcc(dummy, 20, numfilterbunks=20)

    # with c0
    println("-- test c0 is expected")
    cc = mfcc(dummy, 12; czero=true)
    @test length(cc) == 13

    # with c0 + power
    println("-- test c0 and power are expected")
    cc = mfcc(dummy, 12; czero=true, power=true)
    @test length(cc) == 14

    # with power
    println("-- test power is expected")
    cc = mfcc(dummy, 12; czero=false, power=true)
    @test length(cc) == 13
end

# TODO: Update binary dependency for windows
if !is_windows()
println("-- test_mfcc")
for len in [256, 512, 1024, 2048]
    for order in [12, 14, 16, 18]
        println(" where len = $len, order = $order")
        test_mfcc(len, order)
    end
end

test_mfcc_options()
end
