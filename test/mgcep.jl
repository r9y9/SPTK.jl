println("test mel-generalized cepstrum analysis")

function test_mgcep()
    srand(98765)
    dummy_input = rand(256)
    dummy_input_mat = repmat(dummy_input, 1, 2)

    println("-- test_mcep")
    # Since some SPTK functions have static variables inside and
    # this may return different results even if the same input is given
    c1 = mcep(dummy_input, 20, 0.41)
    c2 = mcep(dummy_input, 20, 0.41)
    @test_approx_eq c1 c2

    for order in [20, 22, 24]
        for α in [0.35, 0.41, 0.5]
            c = mcep(dummy_input, order, α)
            @test length(c) == order+1
            cmat = mcep(dummy_input_mat, order, α)
            @test size(cmat) == (order+1, 2)
        end
    end

    # invalid optinal paramters
    @test_throws ArgumentError mcep(dummy_input, itype=-1)
    @test_throws ArgumentError mcep(dummy_input, itype=5)
    @test_throws ArgumentError mcep(dummy_input, eps=-1.0)
    @test_throws ArgumentError mcep(dummy_input, etype=-1)
    @test_throws ArgumentError mcep(dummy_input, etype=-3)
    @test_throws ArgumentError mcep(dummy_input, etype=1, eps=-1.0)
    @test_throws ArgumentError mcep(dummy_input, etype=2, eps=-1.0)
    @test_throws ArgumentError mcep(dummy_input, min_det=-1.0)

    println("-- test_gcep")
    for order in [20, 22, 24]
        for γ in [-1.0, -0.5, 0.0]
            c = gcep(dummy_input, order, γ)
            @test length(c) == order+1
            cmat = gcep(dummy_input_mat, order, γ)
            @test size(cmat) == (order+1, 2)
        end
    end

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

    println("-- test_mgcep")
    for order in [20, 22, 24]
        for α in [0.35, 0.41, 0.5]
            for γ in [-1.0, -0.5, 0.0]
                c = mgcep(dummy_input, order, α, γ)
                @test length(c) == order+1
                cmat = mgcep(dummy_input_mat, order, α, γ)
                @test size(cmat) == (order+1, 2)
            end
        end
    end

    for otype in 1:5
        try mgcep(dummy_input) catch @test false end
    end

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
    @test_throws ArgumentError mgcep(dummy_input, otype=-6)

    println("-- test_uels")
    for order in [20, 22, 24]
        c = uels(dummy_input, order)
        @test length(c) == order+1
        cmat = uels(dummy_input_mat, order)
        @test size(cmat) == (order+1, 2)
    end

    # invalid optinal paramters
    @test_throws ArgumentError uels(dummy_input, itype=-1)
    @test_throws ArgumentError uels(dummy_input, itype=5)
    @test_throws ArgumentError uels(dummy_input, eps=-1.0)
    @test_throws ArgumentError uels(dummy_input, etype=-1)
    @test_throws ArgumentError uels(dummy_input, etype=-3)
    @test_throws ArgumentError uels(dummy_input, etype=1, eps=-1.0)
    @test_throws ArgumentError uels(dummy_input, etype=2, eps=-1.0)

    println("-- test_fftcep")
    for order in [20, 22, 24]
        c = fftcep(dummy_input, order)
        cmat = fftcep(dummy_input_mat, order)
        @test length(c) == order+1
        @test size(cmat) == (order+1, 2)
    end

    println("-- test mcep and mgcep consistency")
    # mgcep when gamma = 0, the result of mgcep is corresponds to mcep
    mc = mcep(dummy_input, 20, 0.41)
    mgc = mgcep(dummy_input, 20, 0.41, 0.0)
    @test_approx_eq_eps mc mgc 1.0e-4

    println("-- test_mgclsp2sp")
    sp = mgclsp2sp(dummy_input, 0.41, -1/4, 512)
    @test length(sp) == 256+1

    # refs #5
    # change order 20 -> 40
    println("-- test changing order (ref: issue #5)")
    mcep(dummy_input, 40, 0.41)
    mgcep(dummy_input, 40, 0.41, 0.0)
end

test_mgcep()
