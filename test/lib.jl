# tests for library routines

import SPTK: agexp, gexp, glog, mseq, theq!, theq, toeplitz!, toeplitz

function test_agexp()
    println("-- test_agexp")
    @test agexp(1, 1, 1) ≈ 5.0
    @test agexp(1, 2, 3) ≈ 18.0
    @test agexp(2, 3, 5) ≈ 12.206555615733702
end

function test_gexp()
    println("-- test_gexp")
    @test gexp(1, 1) ≈ 2.0
    @test gexp(2, 4) ≈ 3.0
end

function test_glog()
    println("-- test_glog")
    @test glog(1, 2) ≈ 1.0
    @test glog(2, 3) ≈ 4.0
end

function test_mseq()
    println("-- test_mseq")
    for i = 1:1000 @test !isnan(mseq()) end
end

function fill_toeplitz!{T}(A::AbstractMatrix{T}, t::AbstractVector{T})
    n = length(t)
    for j=1:n, i=1:n
        A[i,j] = i-j+1 >= 1 ? t[i-j+1] : t[j-i+1]
    end
    A
end

function fill_hankel!{T}(A::AbstractMatrix{T}, h::AbstractVector{T})
    n = length(h)>>1 + 1
    for j=1:n, i=1:n
        A[i,j] = h[i+j-1]
    end
    A
end

# (T + H)a = b
function test_theq()
    println("-- test_theq")
    srand(98765)

    n = 2
    # toeplitz elements
    t = [1.0, -1.0]
    T = zeros(2,2)
    fill_toeplitz!(T, t)

    # hankel elements
    h = zeros(2n-1)
    h[1] = 1.0
    h[3] = 1.0
    H = zeros(2,2)
    fill_hankel!(H, h)

    b = ones(n)
    a = zeros(n)

    # solve (T + H)a = b
    theq!(a, t, h, b)
    @test a ≈ ones(n)
    @test a ≈ theq(t, h, b)
    @test a ≈ (T + H) \ b

    # fail to solve toeplitz plus hankel matrix system
    for m in [3, 5, 10]
        a = ones(m)
        t = ones(a)
        h = ones(2length(a)-1)
        b = ones(a)
        @test_throws ErrorException theq!(a, t, h, b)
    end

    a = ones(5)
    t = ones(a)
    h = ones(2length(a)-1)
    b = ones(a)

    @test_throws ArgumentError theq!(a, t, h, b, min_det=-1.0)
    @test_throws DimensionMismatch theq!(a, t, h, ones(length(b)-1))
    @test_throws DimensionMismatch theq!(a, t, ones(length(h)-1), b)
    @test_throws DimensionMismatch theq!(a, t, ones(length(h)-1), b)
    @test_throws DimensionMismatch theq!(a, ones(length(t)-1), h, b)
    @test_throws DimensionMismatch theq!(ones(length(a)-1), t, h, b)
end

# Ta = b
function test_toeplitz()
    println("-- test_toeplitz")
    srand(98765)

    n = 2
    # toeplitz elements
    t = [2.0, -1.0]
    T = zeros(2,2)
    fill_toeplitz!(T, t)
    b = [1.0, 1.0]
    b = ones(n)
    a = zeros(n)

    # solve Ta = b
    toeplitz!(a, t, b)
    @test a ≈ ones(n)
    @test a ≈ toeplitz(t, b)
    @test a ≈ T \ b

    # fail to solve toeplitz set of linear equations
    for m in [3, 5, 10]
        a = ones(m)
        t = ones(a)
        b = ones(a)
        @test_throws ErrorException toeplitz!(a, t, b)
    end

    a = ones(5)
    t = ones(a)
    b = ones(a)
    @test_throws ArgumentError toeplitz!(a, t, b, eps=-1.0)
    @test_throws DimensionMismatch toeplitz!(a, t, ones(length(b)-1))
    @test_throws DimensionMismatch toeplitz!(a, ones(length(t)-1), b)
    @test_throws DimensionMismatch toeplitz!(ones(length(a)-1), t, b)
end

println("test library routines")
test_agexp()
test_gexp()
test_glog()
test_mseq()
test_theq()
test_toeplitz()
