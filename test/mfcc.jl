using SPTK
using Base.Test

function testmfcc()
  srand(20121)
  dummy = rand(1024)

  # By default, c0 is not contained
  cc = mfcc(dummy, 12)
  @test length(cc) == 12

  # with c0
  cc = mfcc(dummy, 12, czero=true)
  @test length(cc) == 13

  # with c0 + power
  cc = mfcc(dummy, 12, czero=true, power=true)
  @test length(cc) == 14
end

testmfcc()
