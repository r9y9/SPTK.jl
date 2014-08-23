module SPTK

# Julia wrapper for Speech signal Processing Toolkit (SPTK).
#
# The original SPTK: http://sp-tk.sourceforge.net/

# Mel-Cepstrum related analysis
export mcep, mgcep, fftcep, mfcc

# Conversions
export mc2b, b2mc, c2ir, gc2gc, gnorm, ignorm, freqt, mgc2mgc

# F0 estimation
export swipe

deps = joinpath(Pkg.dir("SPTK"), "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("SPTK not properly installed. Please run Pkg.build(\"SPTK\")")
end

include("bridge.jl")

end # module SPTK
