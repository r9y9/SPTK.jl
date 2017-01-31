function test_mat2mat_base(f::Function, M, N)
    srand(98765)
    srcmat = repmat(rand(M), 1, N)
    dstlen = length(f(srcmat[:,1]))
    dstmat = Array{Float64}(dstlen, N)
    for i = 1:N
        dstmat[:,i] = f(srcmat[:,i])
    end
    dstmat2 = f(srcmat)
    @test dstmat2 ≈ dstmat
end

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

    # test_mat2mat_base(mcep, 1024, 2)
    # test_mat2mat_base(mfcc, 512, 2)
    test_mat2mat_base(lpc2c, 20, 2)
    test_mat2mat_base(mc2b, 20, 2)
    test_mat2mat_base(gc2gc, 20, 2)
    test_mat2mat_base(mgc2mgc, 20, 2)
end

println("-- test mat2mat conversion applicable (col-wise)")
test_mat2mat()
