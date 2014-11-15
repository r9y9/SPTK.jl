function mlsadf_test()
    srand(98765)
    dummy_input = rand(1024)
    dummy_ceps = rand(20)

    order=length(dummy_ceps)-1
    pd = 5

    m = MLSADF(order, pd=pd)

    for x in dummy_input
        y = filter!(m, x, dummy_ceps, 0.41)
        @test !isnan(y)
    end
end

function synthesis()
    srand(98765)
    dummy_input = rand(4096)
    order = 20
    timelength = 10
    dummy_ceps = rand(order+1, timelength)

    # MLSADF based synthesis
    m = MLSADF(order)
    result = synthesis!(m, dummy_input, dummy_ceps, 0.41, 80)
    @test !any(isnan(result))

    # MGLSADF based synthesis
    stage = 12
    gamma = float(1/stage)
    g = MGLSADF(order, stage)
    result = synthesis!(g, dummy_input, dummy_ceps, 0.41, 80, gamma)
    @test !any(isnan(result))
end

mlsadf_test()
synthesis()
