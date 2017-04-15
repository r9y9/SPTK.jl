println("test mel-generalized cepstrum analysis")

function test_mgcep_base(f::Function, order, args...; kargs...)
    srand(98765)
    dummy_input = rand(256)
    @assert applicable(f, dummy_input, order, args...)
    c = f(dummy_input, order, args...)
    @test length(c) == order + 1
    @test all(isfinite.(c))
end

function test_mcep_and_mgcep_consistency(order)
    srand(98765)
    dummy_input = rand(256)
    # mgcep when gamma = 0, the result of mgcep is corresponds to mcep
    mc = mcep(dummy_input, order, 0.41)
    mgc = mgcep(dummy_input, order, 0.41, 0.0)
    @test isapprox(mc, mgc, atol=1.0e-3)
end

function test_mcep_exceptions()
    dummy_input = ones(256)
    @test_throws ArgumentError mcep(dummy_input, itype=-1)
    @test_throws ArgumentError mcep(dummy_input, itype=5)
    @test_throws ArgumentError mcep(dummy_input, eps=-1.0)
    @test_throws ArgumentError mcep(dummy_input, etype=-1)
    @test_throws ArgumentError mcep(dummy_input, etype=-3)
    @test_throws ArgumentError mcep(dummy_input, etype=1, eps=-1.0)
    @test_throws ArgumentError mcep(dummy_input, etype=2, eps=-1.0)
    @test_throws ArgumentError mcep(dummy_input, min_det=-1.0)

    # should have error in computing log periodogram
    @test_throws ErrorException mcep(ones(256), 40, 0.41)
end

function test_gcep_exceptions()
    dummy_input = ones(256)

    # invalid γ
    @test_throws ArgumentError gcep(dummy_input, 40, 0.1)
    @test_throws ArgumentError gcep(dummy_input, 40, -2.1)

    # invalid optinal paramters
    @test_throws ArgumentError gcep(dummy_input, itype=-1)
    @test_throws ArgumentError gcep(dummy_input, itype=5)
    @test_throws ArgumentError gcep(dummy_input, eps=-1.0)
    @test_throws ArgumentError gcep(dummy_input, etype=-1)
    @test_throws ArgumentError gcep(dummy_input, etype=-3)
    @test_throws ArgumentError gcep(dummy_input, etype=1, eps=-1.0)
    @test_throws ArgumentError gcep(dummy_input, etype=2, eps=-1.0)
    @test_throws ArgumentError gcep(dummy_input, min_det=-1.0)

    # should have error in theq
    @test_throws ErrorException gcep(ones(256), 40, 0.0)
end

function test_mgcep_exceptions()
    dummy_input = ones(256)

    # invalid γ
    @test_throws ArgumentError mgcep(dummy_input, 40, 0.41, 0.1)
    @test_throws ArgumentError mgcep(dummy_input, 40, 0.41, -2.0)

    # invalid optinal paramters
    @test_throws ArgumentError mgcep(dummy_input, itype=-1)
    @test_throws ArgumentError mgcep(dummy_input, itype=5)
    @test_throws ArgumentError mgcep(dummy_input, eps=-1.0)
    @test_throws ArgumentError mgcep(dummy_input, etype=-1)
    @test_throws ArgumentError mgcep(dummy_input, etype=-3)
    @test_throws ArgumentError mgcep(dummy_input, etype=1, eps=-1.0)
    @test_throws ArgumentError mgcep(dummy_input, etype=2, eps=-1.0)
    @test_throws ArgumentError mgcep(dummy_input, min_det=-1.0)
    @test_throws ArgumentError mgcep(dummy_input, otype=-1)
    @test_throws ArgumentError mgcep(dummy_input, otype=6)

    # should have error in theq
    @test_throws ErrorException mgcep(ones(256))
end

function test_uels_exceptions()
    dummy_input = ones(256)

    # invalid optinal paramters
    @test_throws ArgumentError uels(dummy_input, itype=-1)
    @test_throws ArgumentError uels(dummy_input, itype=5)
    @test_throws ArgumentError uels(dummy_input, eps=-1.0)
    @test_throws ArgumentError uels(dummy_input, etype=-1)
    @test_throws ArgumentError uels(dummy_input, etype=-3)
    @test_throws ArgumentError uels(dummy_input, etype=1, eps=-1.0)
    @test_throws ArgumentError uels(dummy_input, etype=2, eps=-1.0)

    # should have error in computing log periodogram
    @test_throws ErrorException uels(ones(256), 40)
end

function test_lpc_exceptions()
    @test_throws ArgumentError lpc(ones(256), 40, min_det=-1.0)
    @test_throws ErrorException lpc(zeros(256), 40)
end

println("-- test_mcep")
for order in [20, 22, 24]
    for α in [0.35, 0.41, 0.5]
        println(" where order = $order, α = $α")
        test_mgcep_base(mcep, order, α)
    end
end

test_mcep_exceptions()

println("-- test_gcep")
for order in [20, 22, 24]
    for γ in [-1.0, -0.5, 0.0]
        println(" where order = $order, γ = $γ")
        test_mgcep_base(gcep, order, γ)
    end
end

test_gcep_exceptions()

println("-- test_mgcep")
for order in [15, 20, 25, 30]
    for α in [0.35, 0.41, 0.5]
        for γ in [-1.0, -0.5, 0.0]
            for otype in 1:5
                println(" where order = $order, α = $α, γ = $γ, otype=$otype")
                test_mgcep_base(mgcep, order, α, γ, otype=otype)
            end
        end
    end
end

test_mgcep_exceptions()

println("-- test_fftcep")
for order in [20, 22, 24]
    println(" where order = $order")
    test_mgcep_base(fftcep, order)
end

println("-- test_uels")
for order in [20, 22, 24]
    println(" where order = $order")
    test_mgcep_base(uels, order)
end

test_uels_exceptions()

println("-- test_lpc")
for order in [20, 22, 24]
    println(" where order = $order")
    test_mgcep_base(lpc, order)
end


test_lpc_exceptions()

println("-- test mcep and mgcep consistency")
for order in [20, 22, 24]
    println(" where order = $order")
    test_mcep_and_mgcep_consistency(order)
end
