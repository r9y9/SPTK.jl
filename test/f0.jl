function test_f0()
    println("test f0 estimation")
    srand(98765)
    dummy_input = rand(1024)

    println("-- test_swipe")
    f0 = swipe(dummy_input, 16000, 80)
    @test length(f0) == div(length(dummy_input), 80) + 1
    @test all(isfinite(f0))
    @test all(f0 .>= 0)

    @test_throws ArgumentError swipe(dummy_input, 16000, 80, otype=-1)
    @test_throws ArgumentError swipe(dummy_input, 16000, 80, otype=-3)
end

test_f0()
