using SPTK
using Base.Test

for fname in [
              "lib",
              "adaptive",
              "conversions",
              "f0",
              "mfcc",
              "mgcep",
              "synthesis_filters",
              "utils",
              "window"
    ]
    include(string(fname, ".jl"))
end
