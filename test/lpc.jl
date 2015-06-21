function test_lpc()
    println("test lpc")
    srand(98765)
    dummy_input = rand(1024)

    println("-- test_lpc")
    for order in [20, 22, 24]
        l = lpc(dummy_input, order)
        @test length(l) == order + 1
        @test !any(isnan(l))
    end

    @test_throws ArgumentError lpc(dummy_input, 40, min_det=-1.0)

    println("-- test_lpc2c")
    l = lpc(dummy_input, 20)
    c = lpc2c(l)
    @test length(c) == length(l)
    @test !any(isnan(c))

    println("-- test_lpc2lsp")
    lsp = lpc2lsp(l, 20)
    @test length(lsp) == 21
    @test !any(isnan(lsp))

    println("-- test_lpc2par")
    par = lpc2par(l)
    @test length(par) == 21
    @test !any(isnan(par))

    println("-- test_lsp2sp")
    sp = lsp2sp(lsp, 1024)
    @test length(sp) == 1024>>1+1
    @test !any(isnan(sp))

    println("-- test_lspcheck")
    try lspcheck(l); catch @test false; end
end

test_lpc()
