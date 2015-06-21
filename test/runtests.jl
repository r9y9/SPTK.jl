using SPTK
using Base.Test

for fname in [
              "lib",
              "adaptive",
              "conversions",
              "f0",
              "mfcc",
              "mgcep",
              "lpc",
              "synthesis_filters",
              "window",
    ]
    include(string(fname, ".jl"))
end
