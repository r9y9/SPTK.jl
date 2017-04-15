function test_excite()
    println("test excitation generation")
    srand(98765)
    dummy_input = rand(1024)

    println("-- test_excite")
    for hopsize in [40, 80, 160, 320]
        for g in [true, false]
            ex = excite(dummy_input, hopsize, gaussian=g)
            @test all(isfinite.(ex))
        end
    end
end

test_excite()
