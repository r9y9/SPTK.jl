function test_f0_otypes()
    srand(98765)
    dummy_input = rand(Float64, 4096)

    for (otype_int, otype_str) in zip(0:2, SPTK.supported_otypes)
        f01 = swipe(dummy_input, 16000, 80, otype=otype_int)
        f02 = swipe(dummy_input, 16000, 80, otype=otype_str)
        @test f01 ≈ f02
    end

    srand(98765)
    dummy_input = rand(Float32, 4096)

    for (otype_int, otype_str) in zip(0:2, SPTK.supported_otypes)
        f01 = rapt(dummy_input, 16000, 80, otype=otype_int)
        f02 = rapt(dummy_input, 16000, 80, otype=otype_str)
        @test f01 ≈ f02
    end
end

function test_swipe()
    srand(98765)
    dummy_input = rand(1024)

    println("-- test_swipe")
    for otype in [0, 1, 2]
        for fs in [16000]
            for hopsize in [40, 80, 160, 320]
                f0 = swipe(dummy_input, fs, hopsize, otype=otype)
                @test all(isfinite.(f0))
                otype == 0 && @test all(f0 .>= 0)
            end
        end
    end

    fs = 16000

    # invalid otype
    @test_throws ArgumentError swipe(dummy_input, fs, 80, otype=-1)
    @test_throws ArgumentError swipe(dummy_input, fs, 80, otype=-3)

    # invalid min/max
    @test_throws ArgumentError swipe(dummy_input, fs, 80, min=60.0, max=60.0)
    swipe(dummy_input, 16000, 80, min=60.0, max=7999.0)
    @test_throws ArgumentError swipe(dummy_input, fs, 80, min=60.0, max=8000.0)
end

function test_rapt()
    srand(98765)
    dummy_input = rand(Float32, 1024)

    println("-- test_rapt")
    for otype in [0, 1, 2]
        for fs in [16000]
            for hopsize in [40, 80, 160, 320]
                println(" where otype = $otype, fs = $fs and hopsize = $hopsize")
                f0 = rapt(dummy_input, fs, hopsize, otype=otype)
                @test all(isfinite.(f0))
                otype == 0 && @test all(f0 .>= 0)
            end
        end
    end

    fs = 16000

    # invalid otype
    @test_throws ArgumentError rapt(dummy_input, fs, 80, otype=-1)
    @test_throws ArgumentError rapt(dummy_input, fs, 80, otype=-3)

    # invalid min/max
    @test_throws ArgumentError rapt(dummy_input, fs, 80, min=60.0, max=60.0)
    @test_throws ArgumentError rapt(dummy_input, fs, 80, min=60.0, max=8000.0)

    # valid frame period (corner case)
    # TODO: should be 1599 instead of 1600
    rapt(rand(Float32, 10000), fs, 1599)

    # invalid frame_period
    @test_throws ArgumentError rapt(dummy_input, fs, 1601)

    # invalid input length (too small)
    @test_throws ErrorException rapt(dummy_input[1:100], fs, 80)
end

println("test f0 estimation")
test_f0_otypes()
test_swipe()
test_rapt()
