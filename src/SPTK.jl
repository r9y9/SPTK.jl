module SPTK

# Julia wrapper for Speech signal Processing Toolkit (SPTK).
#
# The original SPTK: http://sp-tk.sourceforge.net/
# Note that the SPTK.jl is based on the slightly modified SPTK
# https://github.com/r9y9/SPTK
#
# NOTE: Function interfaces might be different betweeen C and Julia,
# as the following convensions.
#
# Conventions for wrapping C interface
#
# 1. Avoid really short names for variables (e.g. a, b, c, aa, bb, dd) .
#
#    Variable names should be informative. If the C functions have such names,
#    use self-descriptive names for Julia interfaces, unless they have clear
#    meanings in context.
#
# 2. Avoid too many function arguments.
#
#    Less is better. If the C functions have too many function arguments, use
#    keyword arguments with proper default values for optional ones in Julia.
#
# 3. Handle errors in Julia
#
#    Since C functions might `exit` (unfortunately) inside thier functions for
#    unexpected inputs, it should be check if the inputs are supported or not
#    in Julia before `ccall`.

export
    # Library routines
    agexp,
    gexp,
    glog,
    mseq,
    theq!,
    theq,
    toeplitz!,
    toeplitz,

    # Adaptive mel-cepstrum analysis
    acep!,
    agcep!,
    amcep!,

    phidf!,

    # Mel-generalized cepstrum analysis
    mcep,
    gcep,
    mgcep,
    uels,
    fftcep,
    lpc,

    # MFCC
    mfcc,

    # Conversions
    mc2b,
    b2mc,
    c2ir,
    ic2ir,
    c2acr,
    c2ndps,
    ndps2c,
    gc2gc,
    lpc2c,
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
    poledf,
    poledf_delay,
    lmadf,
    lmadf_delay,
    lspdf,
    lspdf_delay,
    ltcdf,
    ltcdf_delay,
    glsadf,
    glsadf_delay,
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

    # Utils
    lspcheck

deps = joinpath(Pkg.dir("SPTK"), "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("SPTK not properly installed. Please run Pkg.build(\"SPTK\")")
end

for fname in [
              "lib",
              "adaptive",
              "conversions",
              "f0",
              "mfcc",
              "mgcep",
              "synthesis_filters",
              "utils",
              "window",
              "extend"
    ]
    include(string(fname, ".jl"))
end

end # module SPTK
