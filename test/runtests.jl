using SPTK
using Base.Test

# MGCep analysis when gamma = 0, MGCep analysis is corresponds to MCep analysis
function mgcep_as_special_case_of_mcep()
  srand(20121)
  dummy = rand(1024)

  mc = mcep(dummy, order=20, alpha=0.41)
  mgc = mgcep(dummy, order=20, alpha=0.41, gamma=0.0)

  # Order + 0-th
  @test length(mc) == 21
  @test length(mgc) == 21

  @test_approx_eq_eps mc mgc 1.0e-4
end

function testmfcc()
  srand(20121)
  dummy = rand(1024)

  # By default, c0 is not contained
  cc = mfcc(dummy, order=12)
  @test length(cc) == 12

  # with c0
  cc = mfcc(dummy, order=12, czero=true)
  @test length(cc) == 13

  # with c0 + power
  cc = mfcc(dummy, order=12, czero=true, power=true)
  @test length(cc) == 14
end

mgcep_as_special_case_of_mcep()
testmfcc()
