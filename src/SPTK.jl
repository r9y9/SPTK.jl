module SPTK

# Julia wrapper for Speech signal Processing Toolkit (SPTK).
# 
# The original SPTK: http://sp-tk.sourceforge.net/

export b2mc, c2ir, gc2gc, gnorm, ignorm, freqt, mc2b, mcep, mgc2mgc, mgcep,
       swipe, mfcc

deps = joinpath(Pkg.dir("SPTK"), "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("SPTK not properly installed. Please run Pkg.build(\"SPTK\")")
end

include("bridge.jl")

end # module SPTK
