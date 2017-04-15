function test_filt_base!(f::Function, order, delay, args...)
    srand(98765)
    dummy_input = rand(1024)
    dummy_sp = abs.(fft(dummy_input))
    dummy_logsp = log.(dummy_sp)
    dummy_mgceps = rand(order + 1)

    @assert applicable(f, first(dummy_input), dummy_mgceps, args..., delay)
    for x in dummy_input
        @test isfinite.(f(x, dummy_mgceps, args..., delay))
    end
end

println("test waveform generation filters")

println("-- test_poledf")
for order in [20,  25, 30]
    println(" where order = $order")
    delay = poledf_delay(order)
    test_filt_base!(poledf, order, delay)
end

let
    @test_throws DimensionMismatch poledf(0.0, ones(10), ones(1))
end

println("-- test_lmadf")
for order in [20,  25, 30]
    for pd in 4:5
        println(" where order = $order, pd = $pd")
        delay = lmadf_delay(order, pd)
        test_filt_base!(lmadf, order, delay, pd)
    end
end

let
    # invalid delay length
    @test_throws DimensionMismatch lmadf(0.0, ones(10), 5, ones(1))
    # invalid pade
    @test_throws ArgumentError lmadf(0.0, ones(10), 3, lmadf_delay(9, 3))
end

println("-- test_lspdf")
for order in [20, 21, 25, 26, 30, 31]
    println(" where order = $order")
    delay = lspdf_delay(order)
    test_filt_base!(lspdf, order, delay)
end

let
    @test_throws DimensionMismatch lspdf(0.0, ones(10), ones(1))
    @test_throws DimensionMismatch lspdf(0.0, ones(9), ones(1))
end

println("-- test_ltcdf")
for order in [20, 21, 25, 26, 30, 31]
    println(" where order = $order")
    delay = ltcdf_delay(order)
    test_filt_base!(ltcdf, order, delay)
end

let
    @test_throws DimensionMismatch ltcdf(0.0, ones(10), ones(1))
end

println("-- test_glsadf")
for order in [20, 21, 25, 26, 30, 31]
    for stage in 1:10
        println(" where order = $order, stage = $stage")
        delay = glsadf_delay(order, stage)
        test_filt_base!(glsadf, order, delay, stage)
    end
end

let
    # invalid delay length
    @test_throws DimensionMismatch glsadf(0.0, ones(10), 1, ones(1))
    # invalid stage
    @test_throws ArgumentError glsadf(0.0, ones(10), 0, glsadf_delay(9, 0))
end

println("-- test_mlsadf")
for order in [20,  25, 30]
    for α in [0.0, 0.35, 0.41, 0.5]
        for pd in 4:5
            println(" where order = $order, α = $α, pd = $pd")
            delay = mlsadf_delay(order, pd)
            test_filt_base!(mlsadf, order, delay, α, pd)
        end
    end
end

let
    # invalid delay length
    @test_throws DimensionMismatch mlsadf(0.0, ones(10), 0.41, 5, ones(1))
    # invalid pade
    @test_throws ArgumentError mlsadf(0.0, ones(10), 0.41, 3, mlsadf_delay(9, 3))
end

println("-- test_mglsadf")
for order in [20,  25, 30]
    for α in [0.0, 0.35, 0.41, 0.5]
        for stage in 1:10
            println(" where order = $order, α = $α, stage = $stage")
            delay = mglsadf_delay(order, stage)
            test_filt_base!(mglsadf, order, delay, α, stage)
        end
    end
end

let
    # invalid delay length
    @test_throws DimensionMismatch mglsadf(0.0, ones(10), 0.41, 1, ones(1))
    # invalid stage
    @test_throws ArgumentError mglsadf(0.0, ones(10), 0.41, 0, glsadf_delay(9, 0))
end

let
    @test applicable(SPTK.poledf_delay_length, 10)
    @test applicable(SPTK.lmadf_delay_length, 10, 5)
    @test applicable(SPTK.lspdf_delay_length, 10)
    @test applicable(SPTK.ltcdf_delay_length, 10)
    @test applicable(SPTK.glsadf_delay_length, 10, 1)
    @test applicable(SPTK.mlsadf_delay_length, 10, 5)
    @test applicable(SPTK.mglsadf_delay_length, 10, 1)
end
