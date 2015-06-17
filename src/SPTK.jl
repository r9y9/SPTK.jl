module SPTK

# Julia wrapper for Speech signal Processing Toolkit (SPTK).
#
# The original SPTK: http://sp-tk.sourceforge.net/
# Note that the SPTK.jl is based on the slightly modified SPTK
# https://github.com/r9y9/SPTK

export
    # library routines
    agexp,
    cholesky!,
    cholesky,
    gexp,
    glog,
    mseq,
    theq!,
    theq,
    toeplitz!,
    toeplitz,

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
    lpc2c,
    lpc2lsp,
    lpc2par,
    lsp2sp,
    lspcheck,
    gnorm,
    ignorm,
    freqt,
    mgc2mgc,
    mgc2sp,
    mgclsp2sp,

    # F0 estimation
    swipe,

    # Waveform generation filters
    poledf,
    poledf_delay,
    lmadf,
    lmadf_delay,
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
    rectangular

deps = joinpath(Pkg.dir("SPTK"), "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("SPTK not properly installed. Please run Pkg.build(\"SPTK\")")
end

for fname in ["lib", "sptk", "extend"]
    include(string(fname, ".jl"))
end

end # module SPTK
