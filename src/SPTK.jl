module SPTK

# Julia wrapper for Speech signal Processing Toolkit (SPTK).
#
# The original SPTK: http://sp-tk.sourceforge.net/
# Note that the SPTK.jl is based on the slightly modified SPTK
# https://github.com/r9y9/SPTK

export
    # (Mel-) Cepstrum related analysis
    mcep,
    gcep,
    mgcep,
    uels,
    fftcep,
    mfcc,

    # Conversions
    mc2b,
    b2mc,
    c2ir,
    c2ndps,
    ndps2c,
    gc2gc,
    lpc,
    lpc2lsp,
    lpc2par,
    lsp2sp,
    gnorm,
    ignorm,
    freqt,
    mgc2mgc,
    mgc2sp,
    mgclsp2sp,

    # F0 estimation
    swipe,

    # Waveform generation filters
    mlsadf,
    mlsadf_delay,
    mglsadf,
    mglsadf_delay,

    # window functions
    blackman,
    hamming,
    hanning,
    bartlett,
    trapezoid,
    rectangular,

    # lib
    theq

deps = joinpath(Pkg.dir("SPTK"), "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("SPTK not properly installed. Please run Pkg.build(\"SPTK\")")
end

for fname in ["bridge", "extend"]
    include(string(fname, ".jl"))
end

end # module SPTK
