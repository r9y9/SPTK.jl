function test_window_functions()
    try SPTK.Cwindow(0) catch @test false end
    try SPTK.Cwindow(5) catch @test false end
    @test_throws ArgumentError SPTK.Cwindow(-1)
    @test_throws ArgumentError SPTK.Cwindow(6)

    println("test windows functions")
    srand(98765)
    x = rand(1024)

    for f in [blackman, hamming, hanning, bartlett, trapezoid, rectangular]
        println("-- test_$f")
        y = f(length(x))
        @test length(y) == length(x)
        @test all(isfinite.(y))
    end

    # invalid normalize flag
    @test_throws ArgumentError blackman(512, normalize=-1)
    @test_throws ArgumentError blackman(512, normalize=3)
end

test_window_functions()
