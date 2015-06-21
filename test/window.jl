function test_window_functions()
    try SPTK.Cwindow(0) catch @test false end
    try SPTK.Cwindow(5) catch @test false end
    @test_throws ArgumentError SPTK.Cwindow(-1)
    @test_throws ArgumentError SPTK.Cwindow(6)

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

    # invalid normalize flag
    @test_throws ArgumentError blackman(512, normalize=-1)
    @test_throws ArgumentError blackman(512, normalize=3)
end

test_window_functions()
