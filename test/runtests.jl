using SPTK
using FFTW
using Base.Test

import SPTK:
     # Adaptive mel-generalized analysis
     acep!,
     agcep!,
     amcep!,

     # Mel-generalized cepstrum analysis
     mcep!,
     mcep,
     gcep!,
     gcep,
     mgcep!,
     mgcep,
     uels!,
     uels,
     fftcep!,
     fftcep,
     lpc!,
     lpc,

     # LPC, LSP and PARCOR conversions
     lpc2c!,
     lpc2c,
     lpc2lsp!,
     lpc2lsp,
     lpc2par!,
     lpc2par,
     par2lpc!,
     par2lpc,
     lsp2sp!,
     lsp2sp,

     # Mel-generalized cepstrum conversions
     mc2b!,
     mc2b,
     b2mc!,
     b2mc,
     b2c!,
     b2c,
     c2acr!,
     c2acr,
     c2ir!,
     c2ir,
     ic2ir!,
     ic2ir,
     c2ndps!,
     c2ndps,
     ndps2c!,
     ndps2c,
     gc2gc!,
     gc2gc,
     gnorm!,
     gnorm,
     ignorm!,
     ignorm,
     freqt!,
     freqt,
     frqtr!,
     frqtr,
     mgc2mgc!,
     mgc2mgc,
     mgc2sp!,
     mgc2sp,
     mgclsp2sp!,
     mgclsp2sp,

     # F0 estimation
     swipe,
     rapt,

     # Excitation generation
     excite,

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

for fname in [
              "common",
              "lib",
              "adaptive",
              "conversions",
              "excite",
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
