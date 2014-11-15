using SPTK
using Base.Test

for fname in ["call",
              "mgcep",
              "mfcc",
              "waveform_generation"]
    include(string(fname, ".jl"))
end
