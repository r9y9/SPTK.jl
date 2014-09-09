using SPTK
using Base.Test

# MGCep analysis when gamma = 0, MGCep analysis is corresponds to MCep analysis
function mgcep_as_special_case_of_mcep()
  srand(20121)
  dummy = rand(1024)

  mc = mcep(dummy, 20, 0.41)
  mgc = mgcep(dummy, 20, 0.41, 0.0)

  # Order + 0-th
  @test length(mc) == 21
  @test length(mgc) == 21

  @test_approx_eq_eps mc mgc 1.0e-4
end

mgcep_as_special_case_of_mcep()
