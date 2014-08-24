using SPTK
using Base.Test

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

# TODO test for real data

mlsadf_test()
