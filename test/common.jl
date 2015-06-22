import SPTK: assert_gamma, assert_pade, assert_fftlen

let
    try assert_gamma(-1) catch @test false end
    try assert_gamma(0) catch @test false end
    @test_throws ArgumentError assert_gamma(-2)
end

let
    try assert_pade(4) catch @test false end
    try assert_pade(5) catch @test false end
    @test_throws ArgumentError assert_pade(3)
end

let
    try assert_fftlen(256) catch @test false end
    try assert_fftlen(512) catch @test false end
    @test_throws ArgumentError assert_fftlen(255)
end
