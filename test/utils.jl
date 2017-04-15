function test_phidf(order, α)
    srand(98765)
    dummy_input = rand(64)

    delay = zeros(order + 1)
    for x in dummy_input
        SPTK.phidf!(x, order, α, delay)
        @test all(isfinite.(delay))
    end
end

function test_lspcheck(order)
    srand(98765)
    dummy_lsp = rand(1024)

    # TODO valid check
    lspcheck(dummy_lsp)
end

println("-- test_phidf")
for order in [20, 22, 24]
    for α in [0.35, 0.41, 0.5]
        println(" where order = $order, α = $α")
        test_phidf(order, α)
    end
end

let
    println("-- test_phidf exceptions")
    c = zeros(21)
    delay = zeros(length(c)-1)
    @test_throws ArgumentError SPTK.phidf!(0.0, length(c)-1, 0.41, delay)
end

println("-- test_lspcheck")
for order in [20, 22, 24]
    println(" where order = $order")
    test_lspcheck(order)
end
