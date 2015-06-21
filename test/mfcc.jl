function test_mfcc()
    println("test_mfcc")
    srand(20121)
    dummy = rand(256)
    dummy_mat = repmat(dummy, 1, 2)

    # By default, c0 is not contained
    println("-- test no c0 is expected")
    for order in [5, 12, 14, 16, 18]
        cc = mfcc(dummy, order)
        @test length(cc) == order
        ccmat = mfcc(dummy_mat, order)
        @test size(ccmat) == (order, 2)
    end

    try mfcc(dummy, 19, numfilterbunks=20) catch @test false end
    @test_throws ArgumentError mfcc(dummy, 20, numfilterbunks=20)

    # with c0
    println("-- test c0 is expected")
    cc = mfcc(dummy, 12; czero=true)
    @test length(cc) == 13
    ccmat = mfcc(dummy_mat, 12; czero=true)
    @test size(ccmat) == (13, 2)

    # with c0 + power
    println("-- test c0 and power are expected")
    cc = mfcc(dummy, 12; czero=true, power=true)
    @test length(cc) == 14
    ccmat = mfcc(dummy_mat, 12; czero=true, power=true)
    @test size(ccmat) == (14, 2)

    # with power
    println("-- test power is expected")
    cc = mfcc(dummy, 12; czero=false, power=true)
    @test length(cc) == 13
    ccmat = mfcc(dummy_mat, 12; czero=false, power=true)
    @test size(ccmat) == (13, 2)
end

test_mfcc()
