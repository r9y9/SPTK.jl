println("test adaptive mel-cepstrum analysis")

function test_acep(order, pd)
    srand(98765)
    dummy_input = rand(64)

    c = zeros(order+1)
    for x in dummy_input
        acep!(c, x, pd=pd)
        @test !any(isnan(c))
    end
end

println("-- test_acep")
for order in [20, 22, 24]
    for pd = 4:5
        println("  where order = $order, pd = $pd")
        test_acep(order, pd)
    end
end

let
    println("-- test_acep exceptions")
    c = zeros(21)
    @test_throws ArgumentError acep!(c, 0.0, pd=3)
    @test_throws ArgumentError acep!(c, 0.0, pd=6)
end

function test_agcep(order, stage)
    srand(98765)
    dummy_input = rand(64)

    c = zeros(order + 1)
    for x in dummy_input
        agcep!(c, x, stage)
        @test !any(isnan(c))
    end
end

println("-- test_agcep")
for order in [20, 22, 24]
    for stage in 1:10
        println("  where order = $order, stage = $stage")
        test_agcep(order, stage)
    end
end

let
    println("-- test_agcep exceptions")
    c = zeros(21)
    @test_throws ArgumentError agcep!(c, 0.0, 0)
    @test_throws ArgumentError agcep!(c, 0.0, -1)
end

function test_amcep(order, α, pd)
    srand(98765)
    dummy_input = rand(64)

    c = zeros(order + 1)
    for x in dummy_input
        amcep!(c, x, α, pd=pd)
        @test !any(isnan(c))
    end
end

println("-- test_amcep")
for order in [20, 22, 24]
    for α in [0.35, 0.41, 0.5]
        for pd = 4:5
            println("  where order = $order, α = $α, pd = $pd")
            test_amcep(order, α, pd)
        end
    end
end

let
    println("-- test_amcep exceptions")
    c = zeros(21)
    @test_throws ArgumentError amcep!(c, 0.0, 0.35, pd=3)
    @test_throws ArgumentError amcep!(c, 0., 0.35, pd=6)
end
