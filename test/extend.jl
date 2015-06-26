function test_mat2mat()
    src = zeros(20)
    srcmat = zeros(20, 2)

    # all extended functions should accept function calls like:
    # 1. f(x) where x is a input vector
    # 2. f(x) where x is a input matrix
    for f in SPTK.vec2vec
        @show f
        @test applicable(eval(f), src)
        @test applicable(eval(f), srcmat)
    end
end

println("-- test mat2mat conversion applicable (col-wise)")
test_mat2mat()
