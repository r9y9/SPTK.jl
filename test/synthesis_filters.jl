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
    d = poledf_delay(order)
    for x in dummy_input
        @test !isnan(poledf(x, dummy_lpc, d))
    end

    println("-- test_lmadf")
    d = lmadf_delay(order, pd)
    for x in dummy_input
        @test !isnan(lmadf(x, dummy_ceps, pd, d))
    end

    println("-- test_lspdf")
    dummy_lsp1 = rand(20) # odd order
    d1 = lspdf_delay(19)
    dummy_lsp2 = rand(21) # even order
    d2 = lspdf_delay(20)
    for x in dummy_input
        @test !isnan(lspdf(x, dummy_lsp1, d1))
        @test !isnan(lspdf(x, dummy_lsp2, d2))
    end

    println("-- test_ltcdf")
    d = ltcdf_delay(order)
    for x in dummy_input
        @test !isnan(ltcdf(x, dummy_ceps, d))
    end

    println("-- test_glsadf")
    for stage in [1, 2, 4, 8]
        d = glsadf_delay(order, stage)
        for x in dummy_input
            @test !isnan(glsadf(x, dummy_ceps, stage, d))
        end
    end

    println("-- test_mlsadf")
    d = mlsadf_delay(order, pd)
    for x in dummy_input
        @test !isnan(mlsadf(x, dummy_ceps, 0.41, pd, d))
    end

    println("-- test_mglsadf")
    for stage in [1, 2, 4, 8, 12]
        d = mglsadf_delay(order, stage)
        for x in dummy_input
            @test !isnan(mglsadf(x, dummy_ceps, 0.41, stage, d))
        end
    end
end

test_waveform_generation_filters()
